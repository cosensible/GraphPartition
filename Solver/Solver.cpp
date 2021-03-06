#include "Solver.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>
#include <cmath>
#include <numeric>
#include <queue>

#include "CsvReader.h"
#include "MpSolver.h"

using namespace std;

namespace szx {

#pragma region Solver::Cli
	int Solver::Cli::run(int argc, char * argv[]) {
		Log(LogSwitch::Szx::Cli) << "parse command line arguments." << endl;
		Set<String> switchSet;
		Map<String, char*> optionMap({ // use string as key to compare string contents instead of pointers.
			{ InstancePathOption(), nullptr },
			{ SolutionPathOption(), nullptr },
			{ RandSeedOption(), nullptr },
			{ TimeoutOption(), nullptr },
			{ MaxIterOption(), nullptr },
			{ JobNumOption(), nullptr },
			{ RunIdOption(), nullptr },
			{ EnvironmentPathOption(), nullptr },
			{ ConfigPathOption(), nullptr },
			{ LogPathOption(), nullptr }
			});

		for (int i = 1; i < argc; ++i) { // skip executable name.
			auto mapIter = optionMap.find(argv[i]);
			if (mapIter != optionMap.end()) { // option argument.
				mapIter->second = argv[++i];
			}
			else { // switch argument.
				switchSet.insert(argv[i]);
			}
		}

		Log(LogSwitch::Szx::Cli) << "execute commands." << endl;
		if (switchSet.find(HelpSwitch()) != switchSet.end()) {
			cout << HelpInfo() << endl;
		}

		if (switchSet.find(AuthorNameSwitch()) != switchSet.end()) {
			cout << AuthorName() << endl;
		}

		Solver::Environment env;
		env.load(optionMap);
		if (env.instPath.empty() || env.slnPath.empty()) { return -1; }

		Solver::Configuration cfg;
		cfg.load(env.cfgPath);

		Log(LogSwitch::Szx::Input) << "load instance " << env.instPath << " (seed=" << env.randSeed << ")." << endl;
		Problem::Input input;
		if (!input.load(env.instPath)) { return -1; }

		Solver solver(input, env, cfg);
		solver.solve();

		pb::Submission submission;
		submission.set_instance(env.friendlyInstName());
		submission.set_partnum(solver.aux.partnum);
		submission.set_obj(solver.output.obj);
		submission.set_imbalance(solver.aux.imbalance);
		submission.set_duration(to_string(solver.timer.elapsedSeconds()) + "s");

		solver.output.save(env.slnPath, submission);
#if SZX_DEBUG
		solver.output.save(env.solutionPathWithTime(), submission);
		solver.record();
#endif // SZX_DEBUG

		return 0;
	}
#pragma endregion Solver::Cli

#pragma region Solver::Environment
	void Solver::Environment::load(const Map<String, char*> &optionMap) {
		char *str;

		str = optionMap.at(Cli::EnvironmentPathOption());
		if (str != nullptr) { loadWithoutCalibrate(str); }

		str = optionMap.at(Cli::InstancePathOption());
		if (str != nullptr) { instPath = str; }

		str = optionMap.at(Cli::SolutionPathOption());
		if (str != nullptr) { slnPath = str; }

		str = optionMap.at(Cli::RandSeedOption());
		if (str != nullptr) { randSeed = atoi(str); }

		str = optionMap.at(Cli::TimeoutOption());
		if (str != nullptr) { msTimeout = static_cast<Duration>(atof(str) * Timer::MillisecondsPerSecond); }

		str = optionMap.at(Cli::MaxIterOption());
		if (str != nullptr) { maxIter = atoi(str); }

		str = optionMap.at(Cli::JobNumOption());
		if (str != nullptr) { jobNum = atoi(str); }

		str = optionMap.at(Cli::RunIdOption());
		if (str != nullptr) { rid = str; }

		str = optionMap.at(Cli::ConfigPathOption());
		if (str != nullptr) { cfgPath = str; }

		str = optionMap.at(Cli::LogPathOption());
		if (str != nullptr) { logPath = str; }

		calibrate();
	}

	void Solver::Environment::loadWithoutCalibrate(const String &filePath) {
		// EXTEND[szx][8]: load environment from file.
		// EXTEND[szx][8]: check file existence first.
	}

	void Solver::Environment::load(const String &filePath) {
		loadWithoutCalibrate(filePath);
		calibrate();
	}

	void Solver::Environment::save(const String &filePath) const {
		// EXTEND[szx][8]: save environment to file.
	}

	void Solver::Environment::calibrate() {
		// adjust thread number.
		int threadNum = thread::hardware_concurrency();
		if ((jobNum <= 0) || (jobNum > threadNum)) { jobNum = threadNum; }

		// adjust timeout.
		msTimeout -= Environment::SaveSolutionTimeInMillisecond;
	}
#pragma endregion Solver::Environment

#pragma region Solver::Configuration
	void Solver::Configuration::load(const String &filePath) {
		// EXTEND[szx][5]: load configuration from file.
		// EXTEND[szx][8]: check file existence first.
	}

	void Solver::Configuration::save(const String &filePath) const {
		// EXTEND[szx][5]: save configuration to file.
	}
#pragma endregion Solver::Configuration

#pragma region Solver

	void Solver::record() const {
#if SZX_DEBUG
		ostringstream log;
		System::MemoryUsage mu = System::peakMemoryUsage();

		// record basic information.
		log << env.friendlyLocalTime() << "," << env.instPath << ","
			<< aux.partnum << "," << output.obj << "," << aux.imbalance << ",";

		log << timer.elapsedSeconds() << "," << env.randSeed << ","
			<< mu.physicalMemory << "," << mu.virtualMemory << endl;

		// append all text atomically.
		static mutex logFileMutex;
		lock_guard<mutex> logFileGuard(logFileMutex);

		ofstream logFile(env.logPath, ios::app);
		logFile.seekp(0, ios::end);
		if (logFile.tellp() <= 0) {
			logFile << "Time,Instance,Partnum,Obj,Imbalance,Duration,RandSeed,PhysMem,VirtMem" << endl;
		}
		logFile << log.str();
		logFile.close();
#endif // SZX_DEBUG
	}

	void Solver::solve() {
		cout << "begin coarsen graph" << endl;
		coarsenGraph();

		cout << "begin initial partition" << endl;
		GraphPartition gp(graphList.back(), aux.partnum, graphList.size() - 1);
		initialPartition(gp);
		cout << "obj of initial sol: " << getObj(gp) << endl;

		cout << "begin optimize initial sol" << endl;
		its(gp);
		cout << "obj of optimized initial sol: " << getObj(gp) << endl;

		cout << "begin uncoarsen graph" << endl;
		uncoarsen(gp);
		output.obj = getObj(gp);
		cout << "Obj = " << output.obj << endl;

		auto &nodepart(*output.mutable_nodepart());
		nodepart.Resize(input.nodes_size(), 0);
		for (int i = 0; i < nodepart.size(); ++i) {
			nodepart[i] = gp.nodesPart[i];
		}
		aux.imbalance = getImbalance(gp);

		cout << "partnum = " << aux.partnum << ", imbalance = " << aux.imbalance << endl;
	}

	template<class ForwardIt, class T>
	ForwardIt binarySearch(ForwardIt first, ForwardIt last, const T& value)
	{
		first = std::lower_bound(first, last, value);
		if (!(first == last) && !(value < *first)) { return first; }
		return last;
	}

	void Solver::coarsenGraph() {
		const auto &inputNodes(*input.mutable_nodes());
		const auto &inputEdges(*input.mutable_edges());
		// protobuf 格式转内部数据结构
		List<int> curNodes(inputNodes.size());
		map<pair<int, int>, int> curEdges; // 或者自定义哈希用 unordered_map
		for (int i = 0; i < inputNodes.size(); ++i) {
			curNodes[i] = inputNodes[i].wgt();
		}
		for (auto it = inputEdges.begin(); it != inputEdges.end(); ++it) {
			curEdges[{it->beg(), it->end()}] = it->wgt();
		}

		while (curNodes.size() > 100) {
			// 构造当前图的邻接表
			auto pGraph = std::make_shared<GraphAdjList>(curNodes, curEdges);
			graphList.push_back(pGraph);
			// 用 HEM 方法求最大匹配
			nodeMap.push_back(List<int>(curNodes.size())); // 添加该层节点到下一层节点的映射
			unordered_set<int> unmatchSet; // 该层尚未匹配的节点集合
			for (int i = 0; i < curNodes.size(); ++i) { unmatchSet.insert(i); }
			int newNodeNum = 0; // 新节点个数
			while ( !unmatchSet.empty()) {
				// 随机挑选未匹配节点 nid
				auto it = unmatchSet.begin();
				for (int index = rand.pick(unmatchSet.size()); index > 0; --index, ++it) {}
				int nid = *it;
				unmatchSet.erase(it);

				auto &adjList(pGraph->nodes[nid].adjList);
				List<AdjNode> candAdjs; // 对应边权大的邻接点作为候选节点
				for (auto it = adjList.begin(); it != adjList.end(); ++it) {
					if (unmatchSet.count(it->adjId) == 0) { continue; } // 邻接点已经匹配
					if (candAdjs.empty() || candAdjs[0].eWgt == it->eWgt) { candAdjs.push_back(*it); }
					else { break; } //碰到对应边权更小的邻接点
				}
				if (!candAdjs.empty()) {
					int index = rand.pick(candAdjs.size());
					int candid = candAdjs[index].adjId;
					unmatchSet.erase(candid);
					nodeMap.back()[nid] = nodeMap.back()[candid] = newNodeNum++;
				}
				else { nodeMap.back()[nid] = newNodeNum++; }
			}

			// 根据最大匹配（节点映射关系）压缩图
			List<int> newNodes(newNodeNum);
			map<pair<int, int>, int> newEdges;
			for (int i = 0; i < curNodes.size(); ++i) {
				newNodes[nodeMap.back()[i]] += curNodes[i];
			}
			for (auto e = curEdges.begin(); e != curEdges.end(); ++e) {
				int beg = nodeMap.back()[e->first.first], end = nodeMap.back()[e->first.second];
				if (beg == end) { continue; }
				newEdges[{beg, end}] += e->second;
			}
			swap(curNodes, newNodes);
			swap(curEdges, newEdges);
		}
		auto pGraph = std::make_shared<GraphAdjList>(curNodes, curEdges);
		graphList.push_back(pGraph);
	}

	void Solver::coarsenGraphImpv() {
		const auto &inputNodes(*input.mutable_nodes());
		const auto &inputEdges(*input.mutable_edges());
		// protobuf 格式转内部数据结构
		List<int> curNodes(inputNodes.size());
		map<pair<int, int>, int> curEdges; // 或者自定义哈希用 unordered_map
		for (int i = 0; i < inputNodes.size(); ++i) {
			curNodes[i] = inputNodes[i].wgt();
		}
		for (auto it = inputEdges.begin(); it != inputEdges.end(); ++it) {
			curEdges[{it->beg(), it->end()}] = it->wgt();
		}

		while (curNodes.size() > 15*aux.partnum) {
			// 构造当前图的邻接表
			auto pGraph = std::make_shared<GraphAdjList>(curNodes, curEdges);
			graphList.push_back(pGraph);
			// 用 HEM* 方法求最大匹配
			nodeMap.push_back(List<int>(curNodes.size())); // 添加该层节点到下一层节点的映射
			unordered_set<int> unmatchSet; // 该层尚未匹配的节点集合
			for (int i = 0; i < curNodes.size(); ++i) { unmatchSet.insert(i); }
			int newNodeNum = 0; // 新节点个数
			while (!unmatchSet.empty()) {
				// 随机挑选未匹配节点 nid
				auto it = unmatchSet.begin();
				for (int index = rand.pick(unmatchSet.size()); index > 0; --index, ++it) {}
				int nid = *it;
				unmatchSet.erase(it);
				
				auto &adjList(pGraph->nodes[nid].adjList);
				List<AdjNode> candAdjs; // 对应边权大的邻接点作为候选节点
				for (auto it = adjList.begin(); it != adjList.end(); ++it) {
					if (unmatchSet.count(it->adjId) == 0) { continue; } // 邻接点已经匹配
					if (candAdjs.empty() || candAdjs[0].eWgt == it->eWgt) { candAdjs.push_back(*it); }
					else { break; } //碰到对应边权更小的邻接点
				}
				if (candAdjs.empty()) { // 无匹配节点
					nodeMap.back()[nid] = newNodeNum++;
					continue;
				}
				auto &nodes(pGraph->nodes);
				int candId = -1, maxWgtSum = -1;
				for (int i = 0; i < candAdjs.size(); ++i) {
					int candWgtSum = 0, cid = candAdjs[i].adjId;
					for (auto &vw : nodes[nid].adjs) {
						if (vw.first == cid) { continue; }
						if (nodes[cid].adjs.count(vw.first) != 0) {
							candWgtSum += nodes[cid].adjs[vw.first];
						}
					}
					if (candWgtSum > maxWgtSum) {
						candId = cid;
						maxWgtSum = candWgtSum;
					}
				}
				unmatchSet.erase(candId);
				nodeMap.back()[nid] = nodeMap.back()[candId] = newNodeNum++;
			}

			// 根据最大匹配（节点映射关系）压缩图
			List<int> newNodes(newNodeNum);
			map<pair<int, int>, int> newEdges;
			for (int i = 0; i < curNodes.size(); ++i) {
				newNodes[nodeMap.back()[i]] += curNodes[i];
			}
			for (auto e = curEdges.begin(); e != curEdges.end(); ++e) {
				int beg = nodeMap.back()[e->first.first], end = nodeMap.back()[e->first.second];
				if (beg == end) { continue; }
				newEdges[{beg, end}] += e->second;
			}
			swap(curNodes, newNodes);
			swap(curEdges, newEdges);
		}
		auto pGraph = std::make_shared<GraphAdjList>(curNodes, curEdges);
		graphList.push_back(pGraph);
	}

	void Solver::uncoarsen(GraphPartition &gp) {
		for (int level = graphList.size() - 2; level >= 0; --level) {
			GraphPartition curgp(graphList[level], aux.partnum, level);
			for (int v = 0; v < curgp.nodeNum; ++v) {
				curgp.nodesPart[v] = gp.nodesPart[nodeMap[level][v]];
			}
			its(curgp);
			swap(curgp, gp);
			cout << "level=" << level << ", obj=" << getObj(gp) << endl;
		}
	}

	void Solver::initialPartition(GraphPartition &gp) {
		List<int> idx(gp.nodeNum);
		for (int i = 0; i < idx.size(); ++i) { idx[i] = i; }

		for (int i = idx.size() - 1; i >= 0; --i) {
			int index = rand.pick(i + 1);
			if (i != index) {
				swap(idx[i], idx[index]);
			}
		}
		for (int i = 0; i < idx.size(); ++i) {
			gp.nodesPart[idx[i]] = i % aux.partnum;
		}
	}


	int Solver::findPart(TabuStruct &tss) {
		int maxPart = 0, maxWgt = tss.partWgts[0];
		int randPart = rand.pick(aux.partnum);
		for (int i = 1; i < aux.partnum; ++i) {
			if (maxWgt < tss.partWgts[i]) {
				maxPart = i;
				maxWgt = tss.partWgts[i];
			}
		}
		while (randPart == maxPart) {
			randPart = rand.pick(aux.partnum);
		}
		return randPart;
	}

	int Solver::findSecondPart(TabuStruct &tss, int firstPart) {
		int maxPart = 0, maxWgt = tss.partWgts[0];
		int randPart = rand.pick(aux.partnum);
		for (int i = 1; i < aux.partnum; ++i) {
			if (maxWgt < tss.partWgts[i]) {
				maxPart = i;
				maxWgt = tss.partWgts[i];
			}
		}
		while (randPart == maxPart || randPart == firstPart) {
			randPart = rand.pick(aux.partnum);
		}
		return randPart;
	}

	List<int> Solver::selectSingleMove(int iter, TabuStruct &tss) {
		int randPart = findPart(tss);
		List<int> cand;
		auto &buckets = tss.bktStruct.bktArrList[randPart].buckets;
		int g = tss.bktStruct.bktArrList[randPart].maxGain;
		while (g >= 0) {
			for (auto it = buckets[g].begin(); it != buckets[g].end(); ++it) {
				if (tss.partWgts[tss.vpmap[*it]] > tss.partWgts[randPart]) {
					// 对目标函数有改进或未被禁忌
					if (tss.curObj - (g - tss.maxIndex) < aux.bestObj || iter >= tss.tabuTable[*it][randPart]) {
						cand.push_back(*it);
					}
				}
			}
			if (!cand.empty()) { break; }
			--g;
		}

		if (cand.empty()) return { -1,randPart,-1 };
		//if (tss.curObj - (g - tss.maxIndex) < aux.bestObj) {
		//	return { cand[rand.pick(cand.size())],randPart,g };
		//}

		int oldest = cand[0];
		for (int i = 1; i < cand.size(); ++i) {
			if (tss.lastMvIter[cand[i]] < tss.lastMvIter[oldest]) { oldest = cand[i]; }
			else if (tss.lastMvIter[cand[i]] == tss.lastMvIter[oldest]) {
				//选移动后差最小的
				int minDiff = abs(tss.partWgts[tss.vpmap[oldest]] - 2 * tss.G->nodes[oldest].vWgt - tss.partWgts[randPart]);
				int diff = abs(tss.partWgts[tss.vpmap[cand[i]]] - 2 * tss.G->nodes[cand[i]].vWgt - tss.partWgts[randPart]);
				if (diff < minDiff) { oldest = cand[i]; }
			}
		}
		return { oldest,randPart,g };
	}

	List<int> Solver::selectSecondMove(int iter, TabuStruct &tss, int firstPart) {
		int secondPart = findSecondPart(tss, firstPart);
		List<int> cand;
		auto &buckets = tss.bktStruct.bktArrList[secondPart].buckets;
		int g = tss.bktStruct.bktArrList[secondPart].maxGain;
		while (g >= 0) {
			for (auto it = buckets[g].begin(); it != buckets[g].end(); ++it) {
				if (tss.vpmap[*it] != firstPart) {
					if (tss.curObj - (g - tss.maxIndex) < aux.bestObj || iter >= tss.tabuTable[*it][secondPart]) {
						cand.push_back(*it);
					}
				}
			}
			if (!cand.empty()) { break; }
			--g;
		}

		if (cand.empty()) return { -1,secondPart,-1 };

		int oldest = cand[0];
		for (int i = 1; i < cand.size(); ++i) {
			if (tss.lastMvIter[cand[i]] < tss.lastMvIter[oldest]) { oldest = cand[i]; }
			else if (tss.lastMvIter[cand[i]] == tss.lastMvIter[oldest]) {
				int minDiff = abs(tss.partWgts[tss.vpmap[oldest]] - 2 * tss.G->nodes[oldest].vWgt - tss.partWgts[secondPart]);
				int diff = abs(tss.partWgts[tss.vpmap[cand[i]]] - 2 * tss.G->nodes[cand[i]].vWgt - tss.partWgts[secondPart]);
				if (diff < minDiff) { oldest = cand[i]; }
			}
		}
		return { oldest,secondPart,g };
	}

	// 将 node 节点从原来的分区移动到 target 分区
	void Solver::execMove(int iter, TabuStruct &tss, List<int> &mv, bool enTabu) {
		if (mv[0] < 0) return;
		int node = mv[0], target = mv[1], gain = mv[2], src = tss.vpmap[node];
		if (gain >= 0) { tss.curObj -= (gain - tss.maxIndex); }
		tss.vpmap[node] = target;
		tss.curParts[src].erase(node);
		tss.curParts[target].insert(node);
		tss.partWgts[src] -= tss.G->nodes[node].vWgt;
		tss.partWgts[target] += tss.G->nodes[node].vWgt;
		tss.maxPartWgt = *std::max_element(tss.partWgts.begin(), tss.partWgts.end());

		if (enTabu) {
			tss.tabuTable[node][src] = iter + tss.borNodesOfPart[src].size() * 0.3 + rand.pick(3);
			tss.lastMvIter[node] = iter;
		}

		// 更新 borNodesOfPart
		tss.borNodesOfPart[target].erase(node);
		for (auto &vw : tss.G->nodes[node].adjs) {
			if (tss.vpmap[vw.first] != target) {
				tss.borNodesOfPart[target].insert(vw.first);
			}
		}
		tss.borNodesOfPart[src].clear();
		for (int v : tss.curParts[src]) {
			for (auto &vw : tss.G->nodes[node].adjs) {
				if (tss.vpmap[vw.first] != src) {
					tss.borNodesOfPart[src].insert(vw.first);
				}
			}
		}

		// 更新 node 节点在桶结构中的位置
		for (int k = 0; k < tss.partNum; ++k) {
			if (k == target) {
				tss.bktStruct.remove(node, target);
			}
			else if (k == src && tss.borNodesOfPart[src].find(node) != tss.borNodesOfPart[src].end()) {
				int gainIndex = tss.getGainIndex(node, src);
				tss.bktStruct.insert(gainIndex, node, src);
			}
			else {
				auto &gainPtr = tss.bktStruct.vptrList[node][k];
				if (gainPtr.first < 0) { continue; } // 不是分区 k 的边界点
				int gainIndex = tss.getGainIndex(node, k);
				tss.bktStruct.insert(gainIndex, node, k);
			}
		}
		// 更新 node 邻接点在桶结构中的位置
		for (auto &vw : tss.G->nodes[node].adjs) {
			for (int k = 0; k < tss.curParts.size(); ++k) {
				if (tss.vpmap[vw.first] == k) { continue; }
				if (tss.borNodesOfPart[k].find(vw.first) == tss.borNodesOfPart[k].end()) {
					tss.bktStruct.remove(vw.first, k);
				}
				else {
					int gainIndex = tss.getGainIndex(vw.first, k);
					tss.bktStruct.insert(gainIndex, vw.first, k);
				}
			}
		}
	}

	void Solver::perturbation(int iter, TabuStruct &tss, int level) {
		int strength = 0.02 * tss.nodeNum;
		int randPart = -1, randNode = -1;
		List<int> mv;
		for (int t = 0; t < strength; ++t) {
			randPart = findPart(tss);
			randNode = rand.pick(tss.nodeNum);
			if (level != 0) {
				while (tss.vpmap[randNode] == randPart || tss.partWgts[tss.vpmap[randNode]] < tss.partWgts[randPart]) {
					randNode = rand.pick(tss.nodeNum);
				}
				mv = { randNode,randPart,-1 };
				execMove(iter, tss, mv);
			}
			else {
				while (tss.vpmap[randNode] == randPart) {
					if (tss.partWgts[tss.vpmap[randNode]] - 2 * tss.G->nodes[tss.vpmap[randNode]].vWgt >= tss.partWgts[randPart])
						break;
					randNode = rand.pick(tss.nodeNum);
				}
				mv = { randNode,randPart,-1 };
				execMove(iter, tss, mv);
			}
		}
		tss.curObj = getObj(tss.G, tss.vpmap); // 更新 curObj
	}

	void Solver::oldPerturb(int iter, TabuStruct &tss) {
		int ptbWgt = 0.02 * tss.nodeNum;
		for (int t = 0; t < ptbWgt; ++t) {
			List<int> candPart;
			for (int k = 0; k < tss.partNum; ++k) {
				if (tss.partWgts[k] != tss.maxPartWgt) { candPart.push_back(k); }
			}
			int randPart = -1;
			if (candPart.empty()) {
				randPart = rand.pick(tss.partNum);
				for (int k = 0; k < tss.partNum; ++k) {
					if (k != randPart) { candPart.push_back(k); }
				}
			}
			else {
				randPart = candPart[rand.pick(candPart.size())];
				candPart.clear();
				for (int k = 0; k < tss.partNum; ++k) {
					if (tss.partWgts[k] > tss.partWgts[randPart]) { candPart.push_back(k); }
				}
			}
			// 等概率选择所有候选节点中的一个
			int randNode = -1, i = 0;
			for (int p : candPart) {
				for (auto it = tss.curParts[p].begin(); it != tss.curParts[p].end(); ++it, ++i) {
					int r = rand.pick(i + 1);
					if (r == 0) { randNode = *it; }
				}
			}
			//execMove(iter, tss, List<int>{randNode, randPart, -1}, false);
			execMove(iter, tss, List<int>{randNode, randPart, -1});
		}
		//tss.tabuTable.assign(tss.nodeNum, List<int>(tss.partNum, 0)); // 扰动后禁忌表置零
		tss.curObj = getObj(tss.G, tss.vpmap); // 更新 curObj
	}

	int Solver::getBalanceQuality(TabuStruct &tss) {
		int bq = 0;
		int div = input.nodes_size() / aux.partnum;
		int vmp = input.nodes_size() % aux.partnum;
		for (int i = 0; i < aux.partnum; ++i) {
			if (vmp != 0 && tss.partWgts[i] != (div + 1) && tss.partWgts[i] != div)
				bq += abs(div - tss.partWgts[i]);
			else if (vmp == 0)
				bq += abs(div - tss.partWgts[i]);
		}
		return bq;
	}

	void Solver::its(GraphPartition &gp) {
		long maxIter = (2 * gp.level + 1) * gp.nodeNum;
		//List<double> alpha = { 0.05,0.1,0.2,0.3 };
		int smallOptCnt = 20, bigOptCnt = 0.1*gp.nodeNum;

		int maxAdjWgt = 0; // 各节点边权和的最大值, 作为桶的最大编号
		for (const auto &node : gp.p2G->nodes) {
			if (node.adjWgt > maxAdjWgt) { maxAdjWgt = node.adjWgt; }
		}
		aux.bestObj = getObj(gp);
		TabuStruct tss(gp, maxAdjWgt, aux.bestObj);
		tss.initDataStrcut();
		int bestBalance = getBalanceQuality(tss);

		int localOpt1 = 0, localOpt2 = 0, noImpv = 0;
		// 迭代禁忌搜索 N1-N2 token-ring
		for (long iter = 0; iter < maxIter; ++iter) {
			List<int> smv = selectSingleMove(iter, tss);
			if (localOpt1 < smallOptCnt) {
				if (aux.bestObj < tss.curObj) { ++localOpt1; }
				execMove(iter, tss, smv);
			}
			else {
				if (aux.bestObj < tss.curObj) { ++localOpt2; }
				List<int> mv2 = selectSecondMove(iter, tss, smv[1]);
				execMove(iter, tss, smv);
				execMove(iter, tss, mv2);
				if (localOpt2 == smallOptCnt) { localOpt1 = localOpt2 = 0; }
			}

			int bq = getBalanceQuality(tss);
			if (tss.curObj <= aux.bestObj && bq <= bestBalance) {
				aux.bestObj = tss.curObj;
				gp.nodesPart = tss.vpmap;
				noImpv = 0;
			}
			if (tss.curObj == aux.bestObj && bq < bestBalance) {
				bestBalance = bq;
			}
			++noImpv;
			if (noImpv == bigOptCnt) { oldPerturb(iter, tss); }
			//if (iter % bigOptCnt == 0) { perturbation(iter, tss, gp.level); }
		}
	}

	int Solver::getObj(GraphPartition &gp) {
		int obj = 0;
		for (int v = 0; v < gp.nodeNum; ++v) {
			for (auto &vw : gp.p2G->nodes[v].adjs) {
				if (gp.nodesPart[v] != gp.nodesPart[vw.first]) {
					obj += vw.second;
				}
			}
		}
		return obj / 2;
	}

	int Solver::getObj(shared_ptr<GraphAdjList> &p2G, const List<int> &nodesPart) {
		int obj = 0;
		for (int v = 0; v < nodesPart.size(); ++v) {
			for (auto &vw : p2G->nodes[v].adjs) {
				if (nodesPart[v] != nodesPart[vw.first]) {
					obj += vw.second;
				}
			}
		}
		return obj / 2;
	}

	double Solver::getImbalance(GraphPartition &gp) {
		vector<int> partWgts(gp.partNum, 0);
		for (int i = 0; i < gp.nodeNum; ++i) {
			partWgts[gp.nodesPart[i]] += gp.p2G->nodes[i].vWgt;
		}
		int maxPartWgt = *std::max_element(partWgts.begin(), partWgts.end());
		int wgtSum = std::accumulate(partWgts.begin(), partWgts.end(), 0);
		int optWgt = 1 + wgtSum / gp.partNum;
		return 1.0*maxPartWgt / optWgt;
	}

#pragma endregion Solver

}
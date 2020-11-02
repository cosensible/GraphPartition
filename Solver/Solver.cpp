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

	void Solver::solve() {
		cout << "begin coarsen graph" << endl;
		coarsenGraph();
		cout << "begin initial partition" << endl;
		GraphPartition gp(graphList.back(), input.partnum());
		initialPartition(gp);
		cout << "begin optimize initial sol" << endl;
		its(gp);
		cout << "begin uncoarsen graph" << endl;
		uncoarsen(gp);
		cout << "Obj = " << getObj(gp) << endl;
	}

	template<class ForwardIt, class T>
	ForwardIt binarySearch(ForwardIt first, ForwardIt last, const T& value)
	{
		first = std::lower_bound(first, last, value);
		if (!(first == last) && !(value < *first)) { return first; }
		return last;
	}

	void Solver::coarsenGraph() {
		const auto &inputNodes(*input.mutable_graph()->mutable_nodes());
		const auto &inputEdges(*input.mutable_graph()->mutable_edges());
		// protobuf 格式转内部数据结构
		List<int> curNodes(inputNodes.size());
		map<pair<int, int>, int> curEdges;
		for (int i = 0; i < inputNodes.size(); ++i) {
			curNodes[i] = inputNodes[i].weight();
		}
		for (auto it = inputEdges.begin(); it != inputEdges.end(); ++it) {
			curEdges[{it->beg(), it->end()}] = it->weight();
		}

		while (curNodes.size() > 200) {
			// 构造当前图的邻接表
			auto pGraph = std::make_shared<GraphAdjList>(curNodes, curEdges);
			graphList.push_back(pGraph);
			// 用 HEM 方法求最大匹配
			nodeMap.push_back(List<int>(curNodes.size())); // 添加该层节点到下一层节点的映射
			// ====> TODO：用 hash_set
			List<int> unmatchSet(curNodes.size()); // 该层尚未匹配的节点集合
			for (int i = 0; i < unmatchSet.size(); ++i) { unmatchSet[i] = i; }
			int newNodeNum = 0; // 新节点个数
			for (; !unmatchSet.empty(); ++newNodeNum) {
				int index = rand.pick(unmatchSet.size());
				int nid = unmatchSet[index];
				unmatchSet.erase(unmatchSet.begin() + index);

				AdjNode *adj = pGraph->nodes[nid].adj;
				auto pos = unmatchSet.end();
				while (adj) {
					// 邻接点已经匹配，考虑下一个邻接点
					if ((pos = binarySearch(unmatchSet.begin(), unmatchSet.end(), adj->adjId)) == unmatchSet.end()) {
						adj = adj->next;
					}
					else { break; }
				}
				if (adj) {
					unmatchSet.erase(pos);
					nodeMap.back()[nid] = nodeMap.back()[adj->adjId] = newNodeNum;
				}
				else { nodeMap.back()[nid] = newNodeNum; }
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
			GraphPartition curgp(graphList[level], input.partnum());
			for (int v = 0; v < curgp.nodeNum; ++v) {
				curgp.nodesPart[v] = gp.nodesPart[nodeMap[level - 1][v]];
			}
			its(curgp);
			swap(curgp, gp);
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
			gp.nodesPart[idx[i]] = i % input.partnum();
		}
	}

	List<int> Solver::selectSingleMove(int iter, TabuStruct &tss) {
		int randPart = rand.pick(tss.partNum - 1);
		if (randPart >= tss.maxPart) { randPart++; }

		List<int> cand;
		auto &buckets = tss.bktStruct.bktArrList[randPart].buckets;
		int g = tss.bktStruct.bktArrList[randPart].maxGain;
		for (; g >= 0 && cand.empty(); --g) {
			for (auto it = buckets[g].begin(); it != buckets[g].end(); ++it) {
				if (tss.partWgts[tss.vpmap[*it]] > tss.partWgts[randPart]) {
					// 对目标函数有改进或未被禁忌
					if (g > buckets.size() / 2 || iter > tss.tabuTable[*it][randPart]) {
						cand.push_back(*it);
					}
				}
			}
		}
		assert(!cand.empty());
		if (g > buckets.size() / 2) {
			return { cand[rand.pick(cand.size())],randPart,g };
		}

		int minMvNode = cand[0];
		for (int i = 1; i < cand.size(); ++i) {
			if (tss.mvFreq[cand[i]] < tss.mvFreq[minMvNode]) { minMvNode = cand[i]; }
			else if (tss.mvFreq[cand[i]] == tss.mvFreq[minMvNode]) {
				//选移动后差最小的
				int minDiff = abs(tss.partWgts[tss.vpmap[minMvNode]] - 2 * tss.G->nodes[minMvNode].vWgt - tss.partWgts[randPart]);
				int diff = abs(tss.partWgts[tss.vpmap[cand[i]]] - 2 * tss.G->nodes[cand[i]].vWgt - tss.partWgts[randPart]);
				if (diff < minDiff) {
					minMvNode = cand[i];
				}
			}
		}
		return { minMvNode,randPart,g };
	}

	List<int> Solver::selectDoubleMove(int iter, TabuStruct &tss) {
		auto mv = selectSingleMove(iter, tss);

		List<int> leftParts;
		for (int k = 0; k < tss.partNum; ++k) {
			if (k != tss.maxPart && k != mv[1]) { leftParts.push_back(k); }
		}
		int randPart = leftParts[rand.pick(leftParts.size())];

		List<int> cand;
		auto &buckets = tss.bktStruct.bktArrList[randPart].buckets;
		int g = tss.bktStruct.bktArrList[randPart].maxGain;
		for (; g >= 0 && cand.empty(); --g) {
			for (auto it = buckets[g].begin(); it != buckets[g].end(); ++it) {
				if (tss.vpmap[*it] != mv[1]) {
					// 对目标函数有改进或未被禁忌
					if (g > buckets.size() / 2 || iter > tss.tabuTable[*it][randPart]) {
						cand.push_back(*it);
					}
				}
			}
		}
		assert(!cand.empty());
		if (g > buckets.size() / 2) {
			mv.push_back(cand[rand.pick(cand.size())]);
		}
		else {
			int minMvNode = cand[0];
			for (int i = 1; i < cand.size(); ++i) {
				if (tss.mvFreq[cand[i]] < tss.mvFreq[minMvNode]) { minMvNode = cand[i]; }
				else if (tss.mvFreq[cand[i]] == tss.mvFreq[minMvNode]) {
					//选移动后差最小的
					int minDiff = abs(tss.partWgts[tss.vpmap[minMvNode]] - 2 * tss.G->nodes[minMvNode].vWgt - tss.partWgts[randPart]);
					int diff = abs(tss.partWgts[tss.vpmap[cand[i]]] - 2 * tss.G->nodes[cand[i]].vWgt - tss.partWgts[randPart]);
					if (diff < minDiff) {
						minMvNode = cand[i];
					}
				}
			}
			mv.push_back(minMvNode);
		}
		mv.push_back(randPart);
		mv.push_back(g);
		return mv;
	}

	// 将 node 节点从原来的分区移动到 target 分区
	void Solver::execMove(TabuStruct &tss, int node, int target, int gain) {
		int src = tss.vpmap[node];
		tss.curObj -= (gain - tss.maxIndex);
		tss.vpmap[node] = target;
		tss.curParts[src].erase(node);
		tss.curParts[target].insert(node);
		tss.partWgts[src] -= tss.G->nodes[node].vWgt;
		tss.partWgts[target] += tss.G->nodes[node].vWgt;
		tss.maxPart = 0;
		for (int k = 1; k < tss.partNum; ++k) {
			if (tss.partWgts[k] > tss.partWgts[tss.maxPart]) {
				tss.maxPart = k;
			}
		}

		// 更新 borNodesOfPart
		tss.borNodesOfPart[target].erase(node);
		for (auto p = tss.G->nodes[node].adj; p; p = p->next) {
			if (tss.vpmap[p->adjId] != target) {
				tss.borNodesOfPart[target].insert(p->adjId);
			}
		}
		tss.borNodesOfPart[src].clear();
		for (int v : tss.curParts[src]) {
			for (auto p = tss.G->nodes[v].adj; p; p = p->next) {
				if (tss.vpmap[p->adjId] != src) {
					tss.borNodesOfPart[src].insert(p->adjId);
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
		for (auto p = tss.G->nodes[node].adj; p; p = p->next) {
			for (int k = 0; k < tss.curParts.size(); ++k) {
				if (tss.vpmap[p->adjId] == k) { continue; }
				if (tss.borNodesOfPart[k].find(p->adjId) == tss.borNodesOfPart[k].end()) {
					tss.bktStruct.remove(p->adjId, k);
				}
				else {
					int gainIndex = tss.getGainIndex(p->adjId, k);
					tss.bktStruct.insert(gainIndex, p->adjId, k);
				}
			}
		}
	}

	void Solver::perturbation(TabuStruct &tss) {
		int ptbWgt = 0.02 * input.graph().nodes_size();
		for (int t = 0; t < ptbWgt; ++t) {
			int randPart = rand.pick(tss.partNum - 1);
			if (randPart >= tss.maxPart) { randPart++; }
			List<int> highPart;
			for (int k = 0; k < tss.partNum; ++k) {
				if (tss.partWgts[k] > tss.partWgts[randPart]) { highPart.push_back(k); }
			}
			// 等概率选择所有候选节点中的一个
			assert(!highPart.empty());
			int randNode = -1, i = 0;
			for (int p : highPart) {
				for (auto it = tss.curParts[p].begin(); it != tss.curParts[p].end(); ++it, ++i) {
					int r = rand.pick(i + 1);
					if (r == 0) { randNode = *it; }
				}
			}
			execMove(tss, randNode, randPart); // 未更新 curObj
		}
		tss.tabuTable.assign(tss.nodeNum, List<int>(tss.partNum, 0)); // 扰动后禁忌表置零
		tss.curObj = getObj(tss.G, tss.vpmap); // 更新 curObj
	}

	void Solver::its(GraphPartition &gp) {
		long iterCount = 10 * input.graph().nodes_size();
		List<double> alpha = { 0.05,0.1,0.2,0.3 };

		int maxAdjWgt = 0; // 各节点边权和的最大值, 作为桶的最大编号
		for (const auto &node : gp.p2G->nodes) {
			if (node.adjWgt > maxAdjWgt) { maxAdjWgt = node.adjWgt; }
		}
		int bestObj = getObj(gp);
		TabuStruct tss(gp, maxAdjWgt, bestObj);
		tss.initDataStrcut();
		int maxPartWgt = tss.partWgts[tss.maxPart];

		// 迭代禁忌搜索 N1-N2 token-ring
		int noImprove1 = 0.005*input.graph().nodes_size(), noImprove2 = noImprove1;
		int noImprove = noImprove1 + noImprove2;
		for (long iter = 0; iter < iterCount;) {
			for (long step = 0; step < noImprove && iter < iterCount;) {
				for (long step1 = 0; step1 < noImprove1 && iter < iterCount;) {
					// 挑选 Move 并执行
					List<int> mvs = selectSingleMove(iter, tss);
					int node = mvs[0], target = mvs[1], gain = mvs[2], src = tss.vpmap[node];
					execMove(tss, node, target, gain);
					tss.tabuTable[node][src] = iter + tss.borNodesOfPart[src].size() * alpha[rand.pick(alpha.size())] + rand.pick(3);
					tss.mvFreq[node]++;
					++step1, ++step, ++iter;
					// at least as good as bestPart in terms of the optimization objective and partition balance
					if (tss.curObj <= bestObj && tss.partWgts[tss.maxPart] <= maxPartWgt) {
						bestObj = tss.curObj;
						maxPartWgt = tss.partWgts[tss.maxPart];
						gp.nodesPart = tss.vpmap;
						step1 = step = 0;
					}
				}
				if (iter >= iterCount || step >= noImprove) { break; }
				for (long step2 = 0; step2 < noImprove2 && iter < iterCount;) {
					List<int> mvs = selectDoubleMove(iter, tss);
					for (int i = 0; i < mvs.size(); i += 3) {
						int node = mvs[i], target = mvs[i + 1], gain = mvs[i + 2], src = tss.vpmap[node];
						execMove(tss, node, target, gain);
						tss.tabuTable[node][src] = iter + tss.borNodesOfPart[src].size() * alpha[rand.pick(alpha.size())] + rand.pick(3);
						tss.mvFreq[node]++;
					}
					++step2, ++step, ++iter;
					if (tss.curObj <= bestObj && tss.partWgts[tss.maxPart] <= maxPartWgt) {
						bestObj = tss.curObj;
						maxPartWgt = tss.partWgts[tss.maxPart];
						gp.nodesPart = tss.vpmap;
						step2 = step = 0;
					}
				}
			}
			if (iter < iterCount) {
				perturbation(tss);
				if (tss.curObj <= bestObj && tss.partWgts[tss.maxPart] <= maxPartWgt) {
					bestObj = tss.curObj;
					maxPartWgt = tss.partWgts[tss.maxPart];
					gp.nodesPart = tss.vpmap;
				}
			}
		}
	}

	int Solver::getObj(GraphPartition &gp) {
		int obj = 0;
		for (int v = 0; v < gp.nodeNum; ++v) {
			for (auto p = gp.p2G->nodes[v].adj; p; p = p->next) {
				if (gp.nodesPart[v] != gp.nodesPart[p->adjId]) {
					obj += p->eWgt;
				}
			}
		}
		return obj;
	}

	int Solver::getObj(shared_ptr<GraphAdjList> &p2G, const List<int> &nodesPart) {
		int obj = 0;
		for (int v = 0; v < nodesPart.size(); ++v) {
			for (auto p = p2G->nodes[v].adj; p; p = p->next) {
				if (nodesPart[v] != nodesPart[p->adjId]) {
					obj += p->eWgt;
				}
			}
		}
		return obj;
	}

#pragma endregion Solver

}
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
	bool Solver::solve() {
		coarsenGraph();
		cout << "curNodes.size = " << aux.curNodes.size()
			<< ", curEdges.size = " << aux.curEdges.size() << endl;
		int obj = optimizeCosGraph(iniSolParts, 1800); // 获得初始解
		cout << "initial solution obj = " << obj << endl;
		if (obj < 0) return false;

		getHighGainNodes(aux.highGainNodes);
		cout << "get high gain nodes finished. size = " << aux.highGainNodes.size() << endl;
		getReservedNodes(aux.highGainNodes, aux.reservedNode);
		cout << "get reserved  nodes finished. size = " << aux.reservedNode.size() << endl;

		cout << "curNodes.size = " << aux.curNodes.size()
			<< ", curEdges.size = " << aux.curEdges.size() << endl;
		obj = optimizeCosGraph(aux.parts, 54000);
		cout << "latest solution obj = " << obj << endl;


		//for (int k = 0; k < iniSolParts.size(); ++k) {
		//	cout << "part = " << k << endl;
		//	for (int n : iniSolParts[k]) {
		//		cout << n << ",";
		//	}
		//	cout << "\n" << endl;
		//}
		return true;
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
		aux.curNodes.resize(inputNodes.size());
		//gidmap.resize(inputNodes.size());
		for (int i = 0; i < inputNodes.size(); ++i) {
			aux.curNodes[i] = inputNodes[i].weight();
			//gidmap[i] = i;
		}

		cout << "input edges size = " << inputEdges.size() << endl;
		for (auto it = inputEdges.begin(); it != inputEdges.end(); ++it) {
			aux.curEdges[{it->beg(), it->end()}] = it->weight();
		}

		cout << "curNodes.size = " << aux.curNodes.size()
			<< ", curEdges.size = " << aux.curEdges.size() << endl;

		// 初始化大图邻接表
		//std::shared_ptr<GraphAdjList> pGraph = std::make_shared<GraphAdjList>(aux.curNodes, aux.curEdges);
		//graphList.push_back(pGraph);
		//graph.init(aux.curNodes, aux.curEdges);

		while (aux.curNodes.size() > 200) {
			// 构造压缩图邻接表
			auto pGraph = std::make_shared<GraphAdjList>(aux.curNodes, aux.curEdges);
			graphList.push_back(pGraph);
			// 用 HEM 方法求最大匹配
			nodeMap.push_back(List<int>());
			nodeMap.back().resize(aux.curNodes.size());
			List<int> unmatchSet(aux.curNodes.size()); // 尚未匹配的节点集合
			for (int i = 0; i < aux.curNodes.size(); ++i) { unmatchSet[i] = i; }
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

			// 根据局部映射表更新全局映射表
			//for (auto it = gidmap.begin(); it != gidmap.end(); ++it) {
			//	*it = idmap[*it];
			//}

			// 根据最大匹配（节点映射关系）压缩图
			List<int> newNodes(newNodeNum);
			map<pair<int, int>, int> newEdges;
			for (int i = 0; i < aux.curNodes.size(); ++i) {
				newNodes[nodeMap.back()[i]] += aux.curNodes[i];
			}
			for (auto e = aux.curEdges.begin(); e != aux.curEdges.end(); ++e) {
				int beg = nodeMap.back()[e->first.first], end = nodeMap.back()[e->first.second];
				if (beg == end) { continue; }
				newEdges[{beg, end}] += e->second;
			}
			// aux.curNodes, aux.curEdges 最后保留压缩后的节点集与边集
			std::swap(aux.curNodes, newNodes);
			std::swap(aux.curEdges, newEdges);
		}
		auto pGraph = std::make_shared<GraphAdjList>(aux.curNodes, aux.curEdges);
		graphList.push_back(pGraph);
	}


	List<int> Solver::initialPartition() {
		const auto &g = *graphList.back();
		List<int> initPartition(g.nodes.size());
		List<int> ids(g.nodes.size());
		for (int i = 0; i < ids.size(); ++i) { ids[i] = i; }

		for (int i = ids.size() - 1; i >= 0; --i) {
			int index = rand.pick(i + 1);
			if (i != index) {
				swap(ids[i], ids[index]);
			}
		}
		for (int i = 0; i < ids.size(); ++i) {
			initPartition[ids[i]] = i % input.partnum();
		}
		return initPartition;
	}

	List<int> Solver::selectSingleMove(int iter, shared_ptr<GraphAdjList> pG, const List<int> &nodesPart,
		const List<List<int>> &tabuTable, const BucketStruct &bktStruct, const List<int> &mvFreq) {
		List<int> partWgts(input.partnum());
		int maxPart = 0;
		for (int i = 0; i < nodesPart.size(); ++i) {
			partWgts[nodesPart[i]] += pG->nodes[i].vWgt;
			if (partWgts[nodesPart[i]] > partWgts[maxPart]) { maxPart = nodesPart[i]; }
		}
		int randPart = rand.pick(input.partnum() - 1);
		if (randPart >= maxPart) { randPart++; }

		auto &buckets = bktStruct.bktArrList[randPart].buckets;
		List<int> cand;
		int g = bktStruct.bktArrList[randPart].maxGain;
		for (; g >= 0 && cand.empty(); --g) {
			for (auto it = buckets[g].begin(); it != buckets[g].end(); ++it) {
				if (partWgts[nodesPart[*it]] > partWgts[randPart]) {
					// 对目标函数有改进或未被禁忌
					if (g > buckets.size() / 2 || iter > tabuTable[*it][randPart]) {
						cand.push_back(*it);
					}
				}
			}
		}
		assert(!cand.empty());
		if (g > buckets.size() / 2) {
			return { cand[rand.pick(cand.size())],randPart };
		}

		int minMvNode = cand[0];
		for (int i = 1; i < cand.size(); ++i) {
			if (mvFreq[cand[i]] < mvFreq[minMvNode]) { minMvNode = cand[i]; }
			else if (mvFreq[cand[i]] == mvFreq[minMvNode]) {
				//选移动后差最小的
				int minDiff = abs(partWgts[nodesPart[minMvNode]] - 2 * pG->nodes[minMvNode].vWgt - partWgts[randPart]);
				int diff = abs(partWgts[nodesPart[cand[i]]] - 2 * pG->nodes[cand[i]].vWgt - partWgts[randPart]);
				if (diff < minDiff) {
					minMvNode = cand[i];
				}
			}
		}
		return { minMvNode,randPart };
	}

	List<int> Solver::selectDoubleMove(int iter, shared_ptr<GraphAdjList> pG, const List<int> &nodesPart,
		const List<List<int>> &tabuTable, const BucketStruct &bktStruct, const List<int> &mvFreq) {
		auto mv = selectSingleMove(iter, pG, nodesPart, tabuTable, bktStruct, mvFreq);

		List<int> partWgts(input.partnum());
		int maxPart = 0;
		for (int i = 0; i < nodesPart.size(); ++i) {
			partWgts[nodesPart[i]] += pG->nodes[i].vWgt;
			if (partWgts[nodesPart[i]] > partWgts[maxPart]) { maxPart = nodesPart[i]; }
		}
		List<int> leftParts;
		for (int i = 0; i < input.partnum(); ++i) {
			if (i != maxPart && i != mv[1]) { leftParts.push_back(i); }
		}
		int randPart = leftParts[rand.pick(leftParts.size())];

		auto &buckets = bktStruct.bktArrList[randPart].buckets;
		List<int> cand;
		int g = bktStruct.bktArrList[randPart].maxGain;
		for (; g >= 0 && cand.empty(); --g) {
			for (auto it = buckets[g].begin(); it != buckets[g].end(); ++it) {
				if (nodesPart[*it] != mv[1]) {
					// 对目标函数有改进或未被禁忌
					if (g > buckets.size() / 2 || iter > tabuTable[*it][randPart]) {
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
				if (mvFreq[cand[i]] < mvFreq[minMvNode]) { minMvNode = cand[i]; }
				else if (mvFreq[cand[i]] == mvFreq[minMvNode]) {
					//选移动后差最小的
					int minDiff = abs(partWgts[nodesPart[minMvNode]] - 2 * pG->nodes[minMvNode].vWgt - partWgts[randPart]);
					int diff = abs(partWgts[nodesPart[cand[i]]] - 2 * pG->nodes[cand[i]].vWgt - partWgts[randPart]);
					if (diff < minDiff) {
						minMvNode = cand[i];
					}
				}
			}
			mv.push_back(minMvNode);
		}
		mv.push_back(randPart);
		return mv;
	}



	int Solver::getGain(shared_ptr<GraphAdjList> pG, List<int> &nodesPart, int nid, int pid) {
		int oldCutWgt = 0, newCutWgt = 0; // 节点被移动前后贡献的切边权重
		for (auto p = pG->nodes[nid].adj; p; p = p->next) {
			if (nodesPart[p->adjId] != nodesPart[nid]) { oldCutWgt += p->eWgt; }
			if (nodesPart[p->adjId] != pid) { newCutWgt += p->eWgt; }
		}
		return oldCutWgt - newCutWgt;
	}

	// 待优化的划分 nodesPart 及其所在的图 pG
	List<int> Solver::its(shared_ptr<GraphAdjList> pG, List<int> &nodesPart) {
		long iterCount = 10 * input.graph().nodes_size();
		List<double> alpha = { 0.05,0.1,0.2,0.3 };
		List<int> bestPart = nodesPart;

		List<unordered_set<int>> curParts(input.partnum());       // 各分区包含的节点
		List<unordered_set<int>> borNodesOfPart(input.partnum()); // 各分区相连的边界节点
		for (int i = 0; i < nodesPart.size(); ++i) {
			curParts[nodesPart[i]].insert(i); // i 所在分区添加 i 节点
			for (auto pAdj = pG->nodes[i].adj; pAdj; pAdj = pAdj->next) {
				if (nodesPart[pAdj->adjId] != nodesPart[i]) {
					// i 的邻居和 i 不在一个分区, 该邻居为 i 所在分区的边界节点
					borNodesOfPart[nodesPart[i]].insert(pAdj->adjId);
				}
			}
		}

		// 禁忌表, 禁止节点 v 移动到 分区 s
		List<List<int>> tabuTable(nodesPart.size(), List<int>(input.partnum()));
		List<int> mvFreq(nodesPart.size()); // 各节点的移动次数

		int maxAdjWgt = 0; // 各节点边权和的最大值, 作为桶的最大编号
		for (const auto &node : pG->nodes) {
			if (node.adjWgt > maxAdjWgt) maxAdjWgt = node.adjWgt;
		}
		// 桶结构, 包含 k 个桶数组和一个指向桶中链表节点的迭代器的列表
		BucketStruct bktStruct(input.partnum(), nodesPart.size(), maxAdjWgt);

		// 计算各分区的边界点移动到分区中能够获得的 gain 值
		for (int k = 0; k < borNodesOfPart.size(); ++k) {
			for (int v : borNodesOfPart[k]) {
				int gainIndex = getGain(pG, nodesPart, v, k) + maxAdjWgt;
				bktStruct.insert(gainIndex, v, k);
			}
		}


		// 将 node 节点从 src 分区移动到 target 分区
		auto execMove = [&](int node, int src, int target) {
			nodesPart[node] = target;
			curParts[target].insert(node);
			curParts[src].erase(node);

			// 更新 borNodesOfPart
			borNodesOfPart[target].erase(node);
			for (auto p = pG->nodes[node].adj; p; p = p->next) {
				if (nodesPart[p->adjId] != target) {
					borNodesOfPart[target].insert(p->adjId);
				}
			}
			borNodesOfPart[src].clear();
			for (int v : curParts[src]) {
				for (auto p = pG->nodes[v].adj; p; p = p->next) {
					if (nodesPart[p->adjId] != src) {
						borNodesOfPart[src].insert(p->adjId);
					}
				}
			}

			// 更新 node 节点在桶结构中的位置
			for (int k = 0; k < curParts.size(); ++k) {
				if (k == target) {
					bktStruct.remove(node, target);
				}
				else if (k == src && borNodesOfPart[src].find(node) != borNodesOfPart[src].end()) {
					int gainIndex = getGain(pG, nodesPart, node, src) + maxAdjWgt;
					bktStruct.insert(gainIndex, node, src);
				}
				else {
					auto &gainPtr = bktStruct.vptrList[node][k];
					if (gainPtr.first < 0) { continue; } // 不是分区 k 的边界点
					int gainIndex = getGain(pG, nodesPart, node, k) + maxAdjWgt;
					bktStruct.insert(gainIndex, node, k);
				}
			}
			// 更新 node 邻接点在桶结构中的位置
			for (auto p = pG->nodes[node].adj; p; p = p->next) {
				for (int k = 0; k < curParts.size(); ++k) {
					if (nodesPart[p->adjId] == k) { continue; }
					if (borNodesOfPart[k].find(p->adjId) == borNodesOfPart[k].end()) {
						bktStruct.remove(p->adjId, k);
					}
					else {
						int gainIndex = getGain(pG, nodesPart, p->adjId, k) + maxAdjWgt;
						bktStruct.insert(gainIndex, p->adjId, k);
					}
				}
			}
		};

		auto perturbation = [&]() {
			int ptbWgt = 0.02 * input.graph().nodes_size();
			for (int t = 0; t < ptbWgt; ++t) {
				List<int> partWgts(input.partnum());
				int maxPart = 0;
				for (int i = 0; i < nodesPart.size(); ++i) {
					partWgts[nodesPart[i]] += pG->nodes[i].vWgt;
					if (partWgts[nodesPart[i]] > partWgts[maxPart]) { maxPart = nodesPart[i]; }
				}
				int randPart = rand.pick(input.partnum() - 1);
				if (randPart >= maxPart) { randPart++; }

				List<int> highPart;
				for (int i = 0; i < input.partnum(); ++i) {
					if (partWgts[i] > partWgts[randPart]) { highPart.push_back(i); }
				}

				//int srcPart = highPart[rand.pick(highPart.size())];
				//assert(!curParts[srcPart].empty());
				//auto it = curParts[srcPart].begin();
				//int randNode = *(it++);
				//for (int i=1; it != curParts[srcPart].end(); ++it,++i) {
				//	int r = rand.pick(i + 1);
				//	if (r == 0) { randNode = *it; }
				//}

				assert(!highPart.empty());
				int srcPart = highPart[0], randNode = -1, i = 0;
				for (int p : highPart) {
					for (auto it = curParts[p].begin(); it != curParts[p].end(); ++it, ++i) {
						int r = rand.pick(i + 1);
						if (r == 0) { randNode = *it; }
					}
				}

				execMove(randNode, srcPart, randPart);
			}
		};

		// 迭代禁忌搜索 N1-N2 token-ring
		for (long iter = 0; iter < iterCount; ++iter) {
			int noImprove1 = 0.005*input.graph().nodes_size(), noImprove2 = noImprove1;
			int noImprove = noImprove1 + noImprove2;

			for (int step = 0; step < noImprove;) {
				for (int step1 = 0; step1 < noImprove1; ++step1, ++step) {
					List<int> singleMv = selectSingleMove(iter, pG, nodesPart, tabuTable, bktStruct, mvFreq);
					for (int i = 0; i < singleMv.size(); i += 2) {
						int node = singleMv[i], target = singleMv[i + 1], src = nodesPart[node];
						execMove(node, src, target);
						// 防止 node 不久后移回 src
						tabuTable[node][src] = iter + borNodesOfPart[src].size() * alpha[rand.pick(alpha.size())] + rand.pick(3);
						mvFreq[node]++;
					}
					//if improved
					step1 = -1; step = -1;
				}
				if (step >= noImprove) break;
				for (int step2 = 0; step2 < noImprove2; ++step2, ++step) {
					//...
					step2 = -1; step = -1;
				}
			}
			//perturbation();
		}
	}

	// 根据curNodes和curEdges求解压缩后各节点所在分区
	// 更新cosNodePart
	int Solver::optimizeCosGraph(List<List<int>> &parts, double timeoutInSec) {
		int nodeNum = aux.curNodes.size(), edgeNum = aux.curEdges.size();
		int partNum = input.partnum();
		aux.cosNodePart.clear();
		aux.cosNodePart.resize(nodeNum);
		MpSolver::Configuration mpCfg(MpSolver::InternalSolver::GurobiMip, timeoutInSec);
		MpSolver mp(mpCfg);

		// add decision variables.
		Arr2D<MpSolver::DecisionVar> isNodeInPart(nodeNum, partNum);
		for (auto x = isNodeInPart.begin(); x != isNodeInPart.end(); ++x) {
			*x = mp.addVar(MpSolver::VariableType::Bool, 0, 1, 0);
		}
		Arr<MpSolver::DecisionVar> isCutEdge(edgeNum);
		for (auto x = isCutEdge.begin(); x != isCutEdge.end(); ++x) {
			*x = mp.addVar(MpSolver::VariableType::Bool, 0, 1, 0);
		}

		// add constraints.
		for (int v = 0; v < nodeNum; ++v) {
			MpSolver::LinearExpr vInOnePart;
			for (int k = 0; k < partNum; ++k) {
				vInOnePart += isNodeInPart[v][k];
			}
			mp.addConstraint(vInOnePart == 1);
		}
		int vWeightSum = std::accumulate(aux.curNodes.begin(), aux.curNodes.end(), 0);
		double Lmax = (1 + input.imbalance()) * (vWeightSum / partNum + 1);
		for (int k = 0; k < partNum; ++k) {
			MpSolver::LinearExpr balance;
			for (int v = 0; v < nodeNum; ++v) {
				balance += isNodeInPart[v][k] * aux.curNodes[v];
			}
			mp.addConstraint(balance <= Lmax);
		}
		int i = 0;
		for (auto x = aux.curEdges.begin(); x != aux.curEdges.end(); ++x, ++i) {
			for (int k = 0; k < partNum; ++k) {
				mp.addConstraint(isCutEdge[i] >= isNodeInPart[x->first.first][k] - isNodeInPart[x->first.second][k]);
				mp.addConstraint(isCutEdge[i] >= isNodeInPart[x->first.second][k] - isNodeInPart[x->first.first][k]);
			}
		}

		// set objective.
		i = 0;
		MpSolver::LinearExpr cutEdgeWeight;
		for (auto x = aux.curEdges.begin(); x != aux.curEdges.end(); ++x, ++i) {
			cutEdgeWeight += x->second*isCutEdge[i];
		}
		mp.addObjective(cutEdgeWeight, MpSolver::OptimaOrientation::Minimize);


		// return result.
		int obj = -1;
		if (mp.optimize()) {
			obj = static_cast<int>(lround(mp.getObjectiveValue()));
			for (int v = 0; v < nodeNum; ++v) {
				for (int k = 0; k < partNum; ++k) {
					if (mp.isTrue(isNodeInPart[v][k])) {
						aux.cosNodePart[v] = k;
						break;
					}
				}
			}
			parts.clear();
			parts.resize(input.partnum());
			for (int i = 0; i < input.graph().nodes_size(); ++i) {
				parts[aux.cosNodePart[gidmap[i]]].push_back(i);
			}
		}
		return obj;
	}

	int Solver::getHighGainNodes(List<int> &highGainNodes, int bottomGain) {
		highGainNodes.clear();
		List<int> srcNodePart(input.graph().nodes_size()); // 原始节点所在分区
		for (int i = 0; i < srcNodePart.size(); ++i) {
			srcNodePart[i] = aux.cosNodePart[gidmap[i]];
		}
		for (int i = 0; i < srcNodePart.size(); ++i) {
			int oldCutWeight = 0; // 节点贡献的切边权重
			for (auto p = graph.nodes[i].adj; p != nullptr; p = p->next) {
				if (srcNodePart[p->adjId] != srcNodePart[i]) { oldCutWeight += p->eWgt; }
			}
			// 考虑将当前节点移动至其他分区
			for (int k = 0; k < input.partnum(); ++k) {
				if (k == srcNodePart[i]) continue; // 节点原本在k区, 不考虑移动
				int newCutWeight = 0;
				for (auto p = graph.nodes[i].adj; p != nullptr; p = p->next) {
					if (srcNodePart[p->adjId] != k) { newCutWeight += p->eWgt; }
				}
				// 是否需要按gain排序, 以便从gain高的节点开始BFS???
				if (oldCutWeight - newCutWeight >= bottomGain) {
					highGainNodes.push_back(i);
					break;
				}
			}
		}
		return highGainNodes.size();
	}

	bool Solver::isUpperNonZero(const List<int> &reservedNodes, const List<List<int>> &parts, int upperNonZeros) {
		List<int> idMap(input.graph().nodes_size(), -1);
		int newNodeNum = 0;
		for (int v : reservedNodes) { idMap[v] = newNodeNum++; }
		for (int p = 0; p < parts.size() && newNodeNum < idMap.size(); p++, newNodeNum++) {
			for (int v : parts[p]) {
				if (idMap[v] == -1) { idMap[v] = newNodeNum; }
			}
		}

		if (newNodeNum == idMap.size()) {
			for (int i = 0; i < newNodeNum; ++i) {
				idMap[i] = i;
			}
		}

		// 根据节点映射关系压缩图
		const auto &inputNodes(*input.mutable_graph()->mutable_nodes());
		const auto &inputEdges(*input.mutable_graph()->mutable_edges());
		aux.curNodes.clear();
		aux.curNodes.resize(newNodeNum);
		aux.curEdges.clear();
		for (int i = 0; i < inputNodes.size(); ++i) {
			aux.curNodes[idMap[i]] += inputNodes[i].weight();
		}
		for (auto it = inputEdges.begin(); it != inputEdges.end(); ++it) {
			int beg = idMap[it->beg()], end = idMap[it->end()];
			if (beg == end) { continue; }
			aux.curEdges[{beg, end}] += it->weight();
		}
		int nonZeros = input.partnum()*(6 * aux.curEdges.size() + 2 * newNodeNum);
		cout << "curEdges.size=" << aux.curEdges.size() << ", curNodes.size=" << newNodeNum << ", nonzeros=" << nonZeros << endl;
		return nonZeros > upperNonZeros;
	}

	int Solver::getReservedNodes(const List<int> &highGainNodes, List<int> &reservedNodes) {
		reservedNodes = highGainNodes;
		// 随机打乱 reservedNodes
		for (int i = reservedNodes.size() - 1; i >= 0; --i) {
			int index = rand.pick(i + 1);
			if (i != index) {
				swap(reservedNodes[i], reservedNodes[index]);
			}
		}

		List<bool> isVisited(input.graph().nodes_size(), false);
		queue<int> nodeQueue;
		for (int n : reservedNodes) {
			nodeQueue.push(n);
			isVisited[n] = true;
		}
		if (isUpperNonZero(reservedNodes, iniSolParts)) { return reservedNodes.size(); }

		while (!nodeQueue.empty()) { // BFS
			int v = nodeQueue.front();
			nodeQueue.pop();
			for (auto p = graph.nodes[v].adj; p != nullptr; p = p->next) {
				if (!isVisited[p->adjId]) {
					nodeQueue.push(p->adjId);
					isVisited[p->adjId] = true;
					reservedNodes.push_back(p->adjId);
					/*if (isUpperNonZero(reservedNodes, iniSolParts)) { return reservedNodes.size(); }*/
				}
			}
			if (isUpperNonZero(reservedNodes, iniSolParts)) { return reservedNodes.size(); }
		}
		return reservedNodes.size();
	}

#pragma endregion Solver

}
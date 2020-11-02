#ifndef SMART_SZX_INVENTORY_ROUTING_SOLVER_H
#define SMART_SZX_INVENTORY_ROUTING_SOLVER_H

#include "Config.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <sstream>
#include <thread>
#include <initializer_list>
#include <memory>
#include <unordered_set>
#include <list>

#include "Common.h"
#include "Utility.h"
#include "LogSwitch.h"
#include "Problem.h"
#include "MpSolver.h"

namespace szx {

	class Solver {
#pragma region Type
	public:
		struct AdjNode {
			int adjId; // 邻接点ID
			int eWgt;  // 边权
			AdjNode *next;
			AdjNode() :adjId(0), eWgt(0), next(nullptr) {}
			AdjNode(int id, int w, AdjNode *nx = nullptr) :adjId(id), eWgt(w), next(nx) {}
		};

		struct Node {
			int vWgt;     // 节点权重
			int adjWgt;   // 节点的边权和
			AdjNode *adj; // 节点邻居链表

			Node() :vWgt(0), adjWgt(0), adj(nullptr) {}
		};

		struct GraphAdjList {
			List<Node> nodes;

			GraphAdjList() = default;
			GraphAdjList(int nodeNum) { nodes.resize(nodeNum); }
			GraphAdjList(const List<int> &ns, const Map<std::pair<int, int>, int> &es) {
				init(ns, es);
			}

			void init(const List<int> &ns, const Map<std::pair<int, int>, int> &es) {
				nodes.resize(ns.size());
				for (int i = 0; i < nodes.size(); ++i) {
					nodes[i].vWgt = ns[i];
				}
				for (auto e = es.begin(); e != es.end(); ++e) {
					// 根据邻接点ID和对应边权构造邻接点
					AdjNode *newAdj = new AdjNode(e->first.second, e->second);
					nodes[e->first.first].adjWgt += e->second;
					// 将该邻接点插入邻接链表头部
					newAdj->next = nodes[e->first.first].adj;
					nodes[e->first.first].adj = newAdj;
					// 寻找新插入节点按边权升序排序后的位置
					AdjNode *adj = nodes[e->first.first].adj, *pos = adj;
					while (pos->next && pos->next->eWgt > adj->eWgt) {
						pos = pos->next;
					}
					if (pos == adj) { continue; }
					nodes[e->first.first].adj = adj->next;
					adj->next = pos->next;
					pos->next = adj;
				}
			}
			~GraphAdjList() {
				for (auto &n : nodes) {
					while (n.adj) {
						AdjNode *cur = n.adj;
						n.adj = n.adj->next;
						delete cur;
					}
				}
			}
		};

		struct BucketArr
		{
			int maxGain;
			List<std::list<int>> buckets;

			BucketArr() :maxGain(0) {}
			BucketArr(int maxIndex) :maxGain(0) {
				buckets.resize(2 * maxIndex + 1);
			}
		};
		struct BucketStruct
		{
			List<BucketArr> bktArrList;

			using GainPtr = std::pair<int, std::list<int>::iterator>;
			List<List<GainPtr>> vptrList;

			BucketStruct() = default;
			BucketStruct(int partNum, int nodeNum, int maxIndex) {
				bktArrList.resize(partNum, BucketArr(maxIndex));
				vptrList.resize(nodeNum, List<GainPtr>(partNum));
				for (int v = 0; v < nodeNum; ++v) {
					for (int k = 0; k < partNum; ++k) {
						vptrList[v][k].first = -1;
					}
				}
			}

			// 将节点 nid 插入分区 pid 对应的桶数组结构中
			void insert(int gainIndex, int nid, int pid) {
				auto &gainPtr = vptrList[nid][pid];
				if (gainPtr.first >= 0) { // 如果已经存在
					if (gainPtr.first == gainIndex) { return; }
					else { remove(nid, pid); }
				}
				auto &buckets = bktArrList[pid].buckets;
				auto it = buckets[gainIndex].insert(buckets[gainIndex].begin(), nid);
				gainPtr = { gainIndex,it };
				if (gainIndex > bktArrList[pid].maxGain) {
					bktArrList[pid].maxGain = gainIndex;
				}
			}

			void remove(int nid, int pid) {
				auto &gainPtr = vptrList[nid][pid];
				if (gainPtr.first < 0) { return; } // 已经被删除

				auto &bkt = bktArrList[pid];
				bkt.buckets[gainPtr.first].erase(gainPtr.second);
				if (gainPtr.first == bkt.maxGain) {
					while (bkt.maxGain > 0 && bkt.buckets[bkt.maxGain].size() == 0) {
						--bkt.maxGain;
					}
				}
				gainPtr.first = -1; // 将指向链表节点的索引置为无效值
			}
		};

		struct GraphPartition {
			int nodeNum = 0;
			int partNum = 0;
			std::shared_ptr<GraphAdjList> p2G;
			List<int> nodesPart;
			GraphPartition(std::shared_ptr<GraphAdjList> p2G, int pn) :nodeNum(p2G->nodes.size()), partNum(pn), p2G(p2G) {
				nodesPart.resize(nodeNum);
			}
		};

		struct TabuStruct {
			int partNum, nodeNum;
			int maxIndex, maxPart; // 每个桶的最大索引, 拥有最大权重的分区
			int curObj;

			std::shared_ptr<GraphAdjList> G;
			List<int> vpmap;	// 各节点所在分区
			// 桶结构, 包含 k 个桶数组和一个指向桶中链表节点的迭代器的列表
			BucketStruct bktStruct;
			List<List<int>> tabuTable;
			List<int> mvFreq, partWgts;
			List<std::unordered_set<int>> curParts;       // 各分区包含的节点
			List<std::unordered_set<int>> borNodesOfPart; // 各分区相连的边界节点

			TabuStruct() = default;
			TabuStruct(const GraphPartition &gp, int mi, int obj) :partNum(gp.partNum), nodeNum(gp.nodeNum),
				maxIndex(mi), maxPart(0), curObj(obj), G(gp.p2G), vpmap(gp.nodesPart), bktStruct(gp.partNum, gp.nodeNum, mi) {
				tabuTable.resize(nodeNum, List<int>(partNum));
				mvFreq.resize(nodeNum);
				partWgts.resize(partNum);
				curParts.resize(partNum);
				borNodesOfPart.resize(partNum);
			}

			void initDataStrcut() {
				for (int i = 0; i < nodeNum; ++i) {
					curParts[vpmap[i]].insert(i); // i 所在分区添加 i 节点
					partWgts[vpmap[i]] += G->nodes[i].vWgt;
					if (partWgts[vpmap[i]] > maxPart) { maxPart = vpmap[i]; }
					for (auto pAdj = G->nodes[i].adj; pAdj; pAdj = pAdj->next) {
						if (vpmap[pAdj->adjId] != vpmap[i]) {
							// i 的邻居和 i 不在一个分区, 该邻居为 i 所在分区的边界节点
							borNodesOfPart[vpmap[i]].insert(pAdj->adjId);
						}
					}
				}
				// 计算各分区的边界点移动到分区中能够获得的 gain 值
				for (int k = 0; k < partNum; ++k) {
					for (int v : borNodesOfPart[k]) {
						int gainIndex = getGainIndex(v, k);
						bktStruct.insert(gainIndex, v, k);
					}
				}
			}

			int getGainIndex(int nid, int pid) {
				int oldCutWgt = 0, newCutWgt = 0; // 节点被移动前后贡献的切边权重
				for (auto p = G->nodes[nid].adj; p; p = p->next) {
					if (vpmap[p->adjId] != vpmap[nid]) { oldCutWgt += p->eWgt; }
					if (vpmap[p->adjId] != pid) { newCutWgt += p->eWgt; }
				}
				return oldCutWgt - newCutWgt + maxIndex;
			}
		};


		using Dvar = MpSolver::DecisionVar;
		using Expr = MpSolver::LinearExpr;

		// commmand line interface.
		struct Cli {
			static constexpr int MaxArgLen = 256;
			static constexpr int MaxArgNum = 32;

			static String InstancePathOption() { return "-p"; }
			static String SolutionPathOption() { return "-o"; }
			static String RandSeedOption() { return "-s"; }
			static String TimeoutOption() { return "-t"; }
			static String MaxIterOption() { return "-i"; }
			static String JobNumOption() { return "-j"; }
			static String RunIdOption() { return "-rid"; }
			static String EnvironmentPathOption() { return "-env"; }
			static String ConfigPathOption() { return "-cfg"; }
			static String LogPathOption() { return "-log"; }

			static String AuthorNameSwitch() { return "-name"; }
			static String HelpSwitch() { return "-h"; }

			static String AuthorName() { return "szx"; }
			static String HelpInfo() {
				return "Pattern (args can be in any order):\n"
					"  exe (-p path) (-o path) [-s int] [-t seconds] [-name]\n"
					"      [-iter int] [-j int] [-id string] [-h]\n"
					"      [-env path] [-cfg path] [-log path]\n"
					"Switches:\n"
					"  -name  return the identifier of the authors.\n"
					"  -h     print help information.\n"
					"Options:\n"
					"  -p     input instance file path.\n"
					"  -o     output solution file path.\n"
					"  -s     rand seed for the solver.\n"
					"  -t     max running time of the solver.\n"
					"  -i     max iteration of the solver.\n"
					"  -j     max number of working solvers at the same time.\n"
					"  -rid   distinguish different runs in log file and output.\n"
					"  -env   environment file path.\n"
					"  -cfg   configuration file path.\n"
					"  -log   activate logging and specify log file path.\n"
					"Note:\n"
					"  0. in pattern, () is non-optional group, [] is optional group\n"
					"     when -env option is not given.\n"
					"  1. an environment file contains information of all options.\n"
					"     explicit options get higher priority than this.\n"
					"  2. reaching either timeout or iter will stop the solver.\n"
					"     if you specify neither of them, the solver will be running\n"
					"     for a long time. so you should set at least one of them.\n"
					"  3. the solver will still try to generate an initial solution\n"
					"     even if the timeout or max iteration is 0. but the solution\n"
					"     is not guaranteed to be feasible.\n";
			}

			// a dummy main function.
			static int run(int argc, char *argv[]);
		};

		// controls the I/O data format, exported contents and general usage of the solver.
		struct Configuration {
			enum Algorithm { Greedy, TreeSearch, DynamicProgramming, LocalSearch, Genetic, MathematicallProgramming };


			Configuration() {}

			void load(const String &filePath);
			void save(const String &filePath) const;


			String toBriefStr() const {
				String threadNum(std::to_string(threadNumPerWorker));
				std::ostringstream oss;
				oss << "alg=" << alg
					<< ";job=" << threadNum;
				return oss.str();
			}


			Algorithm alg = Configuration::Algorithm::Greedy; // OPTIMIZE[szx][3]: make it a list to specify a series of algorithms to be used by each threads in sequence.
			int threadNumPerWorker = (std::min)(1, static_cast<int>(std::thread::hardware_concurrency()));
		};

		// describe the requirements to the input and output data interface.
		struct Environment {
			static constexpr int DefaultTimeout = (1 << 30);
			static constexpr int DefaultMaxIter = (1 << 30);
			static constexpr int DefaultJobNum = 0;
			// preserved time for IO in the total given time.
			static constexpr int SaveSolutionTimeInMillisecond = 1000;

			static constexpr Duration RapidModeTimeoutThreshold = 600 * static_cast<Duration>(Timer::MillisecondsPerSecond);

			static String DefaultInstanceDir() { return "Instance/"; }
			static String DefaultSolutionDir() { return "Solution/"; }
			static String DefaultVisualizationDir() { return "Visualization/"; }
			static String DefaultEnvPath() { return "env.csv"; }
			static String DefaultCfgPath() { return "cfg.csv"; }
			static String DefaultLogPath() { return "log.csv"; }

			Environment(const String &instancePath, const String &solutionPath,
				int randomSeed = Random::generateSeed(), double timeoutInSecond = DefaultTimeout,
				Iteration maxIteration = DefaultMaxIter, int jobNumber = DefaultJobNum, String runId = "",
				const String &cfgFilePath = DefaultCfgPath(), const String &logFilePath = DefaultLogPath())
				: instPath(instancePath), slnPath(solutionPath), randSeed(randomSeed),
				msTimeout(static_cast<Duration>(timeoutInSecond * Timer::MillisecondsPerSecond)), maxIter(maxIteration),
				jobNum(jobNumber), rid(runId), cfgPath(cfgFilePath), logPath(logFilePath), localTime(Timer::getTightLocalTime()) {}
			Environment() : Environment("", "") {}

			void load(const Map<String, char*> &optionMap);
			void load(const String &filePath);
			void loadWithoutCalibrate(const String &filePath);
			void save(const String &filePath) const;

			void calibrate(); // adjust job number and timeout to fit the platform.

			String solutionPathWithTime() const { return slnPath + "." + localTime; }

			String visualizPath() const { return DefaultVisualizationDir() + friendlyInstName() + "." + localTime + ".html"; }
			template<typename T>
			String visualizPath(const T &msg) const { return DefaultVisualizationDir() + friendlyInstName() + "." + localTime + "." + std::to_string(msg) + ".html"; }
			String friendlyInstName() const { // friendly to file system (without special char).
				auto pos = instPath.find_last_of('/');
				String filename = (pos == String::npos) ? instPath : instPath.substr(pos + 1);
				return filename.substr(0, filename.length() - 5); // drop ".json".
			}
			String friendlyLocalTime() const { // friendly to human.
				return localTime.substr(0, 4) + "-" + localTime.substr(4, 2) + "-" + localTime.substr(6, 2)
					+ "_" + localTime.substr(8, 2) + ":" + localTime.substr(10, 2) + ":" + localTime.substr(12, 2);
			}
			double timeoutInSecond() const { return msTimeout / Timer::MillisecondsPerSecond; }

			// essential information.
			String instPath;
			String slnPath;
			int randSeed;
			Duration msTimeout;

			// optional information. highly recommended to set in benchmark.
			Iteration maxIter;
			int jobNum; // number of solvers working at the same time.
			String rid; // the id of each run.
			String cfgPath;
			String logPath;

			// auto-generated data.
			String localTime;
		};

		struct Solution : public Problem::Output { // delivery plan.
			Solution(Solver *pSolver = nullptr) : solver(pSolver) {}
			Solver *solver;
		};
#pragma endregion Type

#pragma region Constant

#pragma endregion Constant

#pragma region Constructor
	public:
		Solver(const Problem::Input &inputData, const Environment &environment, const Configuration &config)
			: input(inputData), env(environment), cfg(config), rand(environment.randSeed),
			timer(std::chrono::milliseconds(environment.msTimeout)), iteration(1) {}
#pragma endregion Constructor

#pragma region Method
	public:
		void solve(); // return true if exit normally. solve by multiple workers together.

	protected:
		void coarsenGraph();
		void initialPartition(GraphPartition &gp);
		void uncoarsen(GraphPartition &gp);

		void its(GraphPartition &gp);
		void perturbation(TabuStruct &tss);
		void execMove(TabuStruct &tss, int node, int target, int gain = 0); // gain=0 表示不更新 obj
		List<int> selectSingleMove(int iter, TabuStruct &tss);
		List<int> selectDoubleMove(int iter, TabuStruct &tss);
		
		int getObj(GraphPartition &gp);
		int getObj(std::shared_ptr<GraphAdjList> &p2G, const List<int> &nodesPart);

#pragma endregion Method

#pragma region Field
	public:
		Problem::Input input;
		Problem::Output output;

		List<std::shared_ptr<GraphAdjList>> graphList;
		List<List<int>> nodeMap;

		Environment env;
		Configuration cfg;
		Random rand; // all random number in Solver must be generated by this.
		Timer timer; // the solve() should return before it is timeout.
		Iteration iteration;
#pragma endregion Field
	};

}


#endif // SMART_SZX_INVENTORY_ROUTING_SOLVER_H
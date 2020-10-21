////////////////////////////////
/// usage : 1.	the entry to the algorithm for the guillotine cut problem.
/// 
/// note  : 1.	
////////////////////////////////

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
			int adjId;   // 邻接点ID
			int eWgt; // 边权
			AdjNode *next;
			AdjNode(int id, int w, AdjNode *nx = nullptr) :adjId(id), eWgt(w), next(nx) {}
		};
		// 方案二: 用vector实现邻居集合
		struct Node {
			int vWgt;  // 节点权重
			int adjWgt;
			AdjNode *adj; // 节点邻居链表
		};
		struct GraphAdjList {
			List<Node> nodes;

			GraphAdjList() = default;
			GraphAdjList(int nodeNum) { nodes = List<Node>(nodeNum); }
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

		struct BucketNode {
			int nid;
			BucketNode *pre;
			BucketNode *next;

			BucketNode() = default;
			BucketNode(int id, BucketNode *pre = nullptr, BucketNode *nx = nullptr) :nid(id), pre(pre), next(nx) {}
		};
		struct Bucket
		{
			int maxGain;
			List<BucketNode*> bucket;

			Bucket() = default;
			Bucket(int maxWgt) {
				bucket.resize(2 * maxWgt + 1, nullptr);
			}
		};
		struct BucketArr
		{
			List<Bucket> buckets;
			List<List<BucketNode*>> nptrArr;

			BucketArr() = default;
			BucketArr(int partNum, int nodeNum, int maxWgt) {
				buckets.resize(partNum, Bucket(maxWgt));
				nptrArr.resize(nodeNum, List<BucketNode*>(partNum));
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
		bool solve(); // return true if exit normally. solve by multiple workers together.

	protected:
		void coarsenGraph();
		List<int> initialPartition();
		List<int> its(std::shared_ptr<GraphAdjList> pG, List<int> &curPart);
		List<int> selectSingleMove(int iter, std::shared_ptr<GraphAdjList> pG, const List<int> &curParts,
			const List<List<int>> &tabuList, const BucketArr &bucketArr, const List<int> &mvFreq);
		List<int> selectDoubleMove(int iter, std::shared_ptr<GraphAdjList> pG, const List<int> &curParts,
			const List<List<int>> &tabuList, const BucketArr &bucketArr, const List<int> &mvFreq);


		int optimizeCosGraph(List<List<int>> &parts, double timeoutInSec = 600);
		// 返回 gain >= bottomGain 的节点集合 highGainNodes
		int getHighGainNodes(List<int> &highGainNodes, int bottomGain = -2);
		// 当前被保留的节点集是否达到指定大小
		bool isUpperNonZero(const List<int> &reservedNodes, const List<List<int>> &parts, int upperNonZeros = 1e6);
		// 返回不被压缩的节点集合 reservedNodes
		int getReservedNodes(const List<int> &highGainNodes, List<int> &reservedNodes);

		

#pragma endregion Method

#pragma region Field
	public:
		Problem::Input input;
		Problem::Output output;

		List<std::shared_ptr<GraphAdjList>> graphList;
		List<List<int>> nodeMap;


		GraphAdjList graph;
		List<int> gidmap;
		List<List<int>> iniSolParts; // 初始解各分区及其节点

		struct {
			List<int> cosNodePart; // 压缩节点所在分区
			List<int> curNodes;
			Map<std::pair<int, int>, int> curEdges;

			List<int> highGainNodes, reservedNode;
			List<List<int>> parts;
		} aux;

		Environment env;
		Configuration cfg;
		Random rand; // all random number in Solver must be generated by this.
		Timer timer; // the solve() should return before it is timeout.
		Iteration iteration;
#pragma endregion Field
	};

}


#endif // SMART_SZX_INVENTORY_ROUTING_SOLVER_H
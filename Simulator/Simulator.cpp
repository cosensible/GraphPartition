#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <random>
#include <cstring>

#include "Simulator.h"
#include "ThreadPool.h"

using namespace std;

namespace szx {

	// EXTEND[szx][5]: read it from InstanceList.txt.
	static const vector<String> instList({
		"add20",    "data",        "3elt",      "uk",        "add32",     "bcsstk33",  "whitaker3",
		"crack",    "wing_nodal",  "fe_4elt2",  "vibrobox",  "bcsstk29",  "4elt",      "fe_sphere",
		"cti",      "memplus",     "cs4",       "bcsstk30",  "bcsstk31",  "fe_pwt",    "bcsstk32",
		"fe_body",  "t60k",        "wing",      "brack2",    "finan512",  "fe_tooth",  "fe_rotor",
		"598a",     "fe_ocean"
		});

	void Simulator::initDefaultEnvironment() {
		Solver::Environment env;
		env.save(Env::DefaultEnvPath());

		Solver::Configuration cfg;
		cfg.save(Env::DefaultCfgPath());
	}

	void Simulator::exe(const Task &task) {
		System::makeSureDirExist(SolutionDir());

		ostringstream oss;
		oss << ProgramName()
			<< " " << Cmd::InstancePathOption() << " " << InstanceDir() << task.instanceName()
			<< " " << Cmd::SolutionPathOption() << " " << SolutionDir() << task.solutionName();

		auto addOption = [&](const String &key, const String &value) {
			if (!value.empty()) { oss << " " << key << " " << value; }
		};

		addOption(Cmd::RandSeedOption(), task.randSeed);
		addOption(Cmd::TimeoutOption(), task.timeout);
		addOption(Cmd::MaxIterOption(), task.maxIter);
		addOption(Cmd::JobNumOption(), task.jobNum);
		addOption(Cmd::RunIdOption(), task.runId);
		addOption(Cmd::ConfigPathOption(), task.cfgPath);
		addOption(Cmd::LogPathOption(), task.logPath);

		System::exec(oss.str());
	}

	void Simulator::run(const Task &task) {
		System::makeSureDirExist(SolutionDir());

		char argBuf[Cmd::MaxArgNum][Cmd::MaxArgLen];
		char *argv[Cmd::MaxArgNum];
		for (int i = 0; i < Cmd::MaxArgNum; ++i) { argv[i] = argBuf[i]; }
		strcpy(argv[ArgIndex::ExeName], ProgramName().c_str());

		int argc = ArgIndex::ArgStart;

		strcpy(argv[argc++], Cmd::InstancePathOption().c_str());
		strcpy(argv[argc++], (InstanceDir() + task.instanceName()).c_str());
		strcpy(argv[argc++], Cmd::SolutionPathOption().c_str());
		strcpy(argv[argc++], (SolutionDir() + task.solutionName()).c_str());

		auto addOption = [&](const String &key, const String &value) {
			if (!value.empty()) {
				strcpy(argv[argc++], key.c_str());
				strcpy(argv[argc++], value.c_str());
			}
		};

		addOption(Cmd::RandSeedOption(), task.randSeed);
		addOption(Cmd::TimeoutOption(), task.timeout);
		addOption(Cmd::MaxIterOption(), task.maxIter);
		addOption(Cmd::JobNumOption(), task.jobNum);
		addOption(Cmd::RunIdOption(), task.runId);
		addOption(Cmd::ConfigPathOption(), task.cfgPath);
		addOption(Cmd::LogPathOption(), task.logPath);

		Cmd::run(argc, argv);
	}

	void Simulator::run(const String &envPath) {
		char argBuf[Cmd::MaxArgNum][Cmd::MaxArgLen];
		char *argv[Cmd::MaxArgNum];
		for (int i = 0; i < Cmd::MaxArgNum; ++i) { argv[i] = argBuf[i]; }
		strcpy(argv[ArgIndex::ExeName], ProgramName().c_str());

		int argc = ArgIndex::ArgStart;

		strcpy(argv[argc++], Cmd::EnvironmentPathOption().c_str());
		strcpy(argv[argc++], envPath.c_str());

		Cmd::run(argc, argv);
	}

	void Simulator::debug() {
		Task task;
		task.instSet = "";
		task.instId = "add20";
		task.timeout = "600";
		//task.randSeed = "1559429277";
		task.randSeed = to_string(Random::generateSeed());
		task.jobNum = "1";
		task.cfgPath = Env::DefaultCfgPath();
		task.logPath = Env::DefaultLogPath();
		task.runId = "0";

		run(task);
	}

	void Simulator::benchmark(int repeat) {
		Task task;
		task.instSet = "";
		task.timeout = "300";
		task.jobNum = "1";
		task.cfgPath = Env::DefaultCfgPath();
		task.logPath = Env::DefaultLogPath();

		for (int i = 0; i < repeat; ++i) {
			for (auto inst = instList.begin(); inst != instList.end(); ++inst) {
				task.instId = *inst;
				task.randSeed = to_string(Random::generateSeed());
				task.runId = to_string(i);
				run(task);
			}
		}
	}

	void Simulator::parallelrun(Task task, int repeat) {
		for (int i = 0; i < repeat; ++i) {
			task.runId = to_string(i);
			task.randSeed = to_string(Random::generateSeed());
			exe(task);
			this_thread::sleep_for(2s);
		}
	}

	void Simulator::parallelBenchmark(int repeat) {
		Task task;
		task.instSet = "";
		task.timeout = "3600";
		task.jobNum = "1";
		task.cfgPath = Env::DefaultCfgPath();
		task.logPath = Env::DefaultLogPath();

		ThreadPool<> tp(7);
		random_device rd;
		mt19937 rgen(rd());

		//for (auto inst = instList.begin(); inst != instList.end(); ++inst) {
		for (auto inst = instList.rbegin(); inst != instList.rend(); ++inst) {
			task.instId = *inst;
			tp.push([&, task, repeat]() { parallelrun(task, repeat); });
			this_thread::sleep_for(15s);
		}
	}

}

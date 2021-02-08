// Wrapper for UAI competition code  -- For qlou solvers (?)
// 

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <boost/program_options.hpp>

#include "graphmodel.h"
#include "wmbe.h"
#include "mcts.h"
#include "wmbsearch.h"
#include "searchsample.h"
#include "mmap.h"
#include "randheur.h"
#include "mmapIS.h"

using namespace std;
using mex::mxObject;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::vector;
using mex::graphModel;

using mex::timeSystem;

namespace po = boost::program_options;

#define c_log10 std::log(10)

double wmbMemLimit;
double dt;
int nOrders, nExtra;
size_t iBound;
double timeOrder;
double stopIter, stopTime;
int nSearch;

double dampTheta, stepWeights;
double addEpsilon;
int isAdaptive;
double memUseRandom = mex::infty(); // if memory *way* out of reach, just use random ordering...
int verbose; // output messages
//int max_nodes = 1;
// memory limit for search
double searchMemLimit;
// log ratio between upper and lower bounds
//double searchBoundRatio;
// output time unit
double timeUnit;
// output time ratio
// timeUnit * timeRatio^t is the actual timestamp
// if timeRatio > 1, its output is in exponential
double timeRatio;
// no. of samples
int nSample;
double DELTA;
// overall time budget
double timeBudget;
double searchTime;
// overall memory budget including wmb construction and search
double memoryBudget;
// overall memory control, let wmb construction use as memory as it can (within memory budget)
// the left memory is for other task
int isAllMemControl = 0;
// overall time control
int isAllTimeControl = 0;
// for generating a tree
int NSP;
// no. of nodes to expand in each iteration for searchsample
int NND;
// no. of samples to draw in each expansion of a frontier node
int KSP;
//
int lookahead_depth;
//
int sumVarCopies; // no. of copies of SUM variables, used in "MMAPIS"

// Caveat: MEX_ENUM seems buggy: it can not deal with space between the second last and the last options
MEX_ENUM( Task, MPE,PR,MAR,MMAP );
MEX_ENUM( Method, AOBFS,UBFS,MCTS,SESA,RDHR,MMAPIS );
//MEX_ENUM( Mode , DFS,SMA ); defined in wmbsearch.h
MEX_ENUM( Heur, UB,BOTH ); // upper bound heuristics, or both upper and lower heuristics, default BOTH

bool isAndOr = true; // whether using AND/OR structure
// whether we need lower heuristics
bool isLbHeur = true;
bool isSolved = false;

int main(int argc, char* argv[]) {
//  double timeStart = timeSystem();
	//if (argc < 4) { cout<<"Usage: "<<argv[0]<<" <file.uai> <seed> <PR|MPE|MAR>\n"; return 0; }
	const char* probName; // = argv[1];
	const char* taskName; //  = argv[3];
	const char* methodName;
	const char* pstName; // pseudo tree file
	const char* priorityName; // type of priority
	const char* searchModeName; // SMA or DFS
	const char* heurName;

	mex::vector<Factor> bel;

	Task task;
	Method method;
	Mode searchMode; // SMA or DFS
	Heur heur;

	po::options_description desc("Available options");
	desc.add_options()
			("help", "print help message")
			("file,f", po::value<std::string>(), "input problem filename")
			("evidence,e", po::value<std::string>(), "input evidence filename")
			("seed,S", po::value<int>(), "random number initial seed")
			("task,T", po::value<std::string>(), "inference task string: PR, MMAP")
			("heur,H", po::value<std::string>(), "heuristics, UB (upper bound) or BOTH (default)")
			("method", po::value<std::string>(), "inference methods: AOBFS, UBFS, MCTS, SESA, RDHR, MMAPIS")
			("orders,o", po::value<int>(&nOrders)->default_value(1), "number of variable orderings to try")
			("order-time,t", po::value<double>(&timeOrder)->default_value(1), "max time spend on variable orderings")
			("order-rand", po::value<int>(&nExtra)->default_value(0), "var order randomness; n=-1 (none), or among best+n")
			("order-file", po::value<std::string>(), "problem elimination ordering filename")
			("query-file", po::value<std::string>(), "query (MAX) variables, see uai-14 query file format")
			("ibound,i", po::value<size_t>(&iBound)->default_value(50), "initial i-bound")
			("memory,m", po::value<double>(&wmbMemLimit)->default_value(2 * 1024.0), "memory bound for WMB construction (MB)")
			("stopiter,n", po::value<double>(&stopIter)->default_value(1.0), "maximum number of iterations")
			("stoptime", po::value<double>(&stopTime)->default_value(20.0), "maximum iteration time (sec)")
			("damp", po::value<double>(&dampTheta)->default_value(1.0), "damping of factor reparameterization, real[0,1]: 0.0=no update, 1.0=no damping")
			("stepw", po::value<double>(&stepWeights)->default_value(0.1), "step size of weight update; 0.0=no update")
			("adaptive", po::value<int>(&isAdaptive)->default_value(1), "Adaptive step sizes; 0=no, 1=yes (default)")
			("eps", po::value<double>(&addEpsilon)->default_value(0.0), "constant to add to all factors (to avoid hard zeros). require for mmapIS sampling!")
			("nsearch,s", po::value<int>(&nSearch)->default_value(0), "max number of search nodes")
			("nsample", po::value<int>(&nSample)->default_value(1e8), "max number of samples")
			("delta", po::value<double>(&DELTA)->default_value(0.025), "delta  (1-confidence) for Empirical Bernstein bound")
			("pseudo-tree", po::value<bool>(&isAndOr)->default_value(1), "non-chain pseudo tree? default 1 (yes)")
			("timebudget", po::value<double>(&timeBudget)->default_value(60.0), "overall time budget for all processes")
			("priority", po::value<std::string>(), "priority type")
			("searchmemory", po::value<double>(&searchMemLimit)->default_value(2 * 1024.0), "memory limit for wmbsearch (MB)")
			("memorybudget", po::value<double>(&memoryBudget)->default_value(4 * 1024.0), "overall memory control")
			("timeunit", po::value<double>(&timeUnit)->default_value(60.0), "output period (second)")
			("searchtime", po::value<double>(&searchTime)->default_value(60.0), "time limit for wmbsearch (sec)")
			("pst_file", po::value<std::string>(), "output pseudo tree filename")
			("verbose,v", po::value<int>(&verbose)->default_value(0), "verbose, 0,1,2,3")
			("searchmode", po::value<std::string>(), "search mode in wmbsearch: SMA (default), DFS")
			("nsp", po::value<int>(&NSP)->default_value(-1), "no. of samples to draw in each iteration")
			("nnd", po::value<int>(&NND)->default_value(1), "no. of nodes to expand in each iteration")
			("ksp", po::value<int>(&KSP)->default_value(1), "no. of samples to draw in each expansion")
			("sumCopies,K", po::value<int>(&sumVarCopies)->default_value(1), "no. of copies of SUM variables")
			("lookahead", po::value<int>(&lookahead_depth)->default_value(std::numeric_limits<int>::infinity()),
					"depth of stochastic lookahead, 0 - adaptive")
			("timecontrol", po::value<int>(&isAllTimeControl)->default_value(1), "set 0 to disable overall time control")
			("memorycontrol", po::value<int>(&isAllMemControl)->default_value(0), "set > 0 to enable overall memory control")
			;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}
	if (vm.count("file")) {
		probName = vm["file"].as<std::string>().c_str();
	} else {
		std::cout << "Missing input problem file!\n";
		return 1;
	}
	if (vm.count("seed")) {
		mex::randSeed(vm["seed"].as<int>());
	}
	// task
	if (vm.count("task")) {
		taskName = vm["task"].as<std::string>().c_str();
		task = Task(taskName);
	} else {
		std::cout << "Missing task!\n";
		return 1;
	}
	if (task != Task::PR && task != Task::MMAP) {
		std::cout << "Task should be PR or MMAP" << std::endl;
		return 1;
	}
	// method
	if (vm.count("method")) {
		methodName = vm["method"].as<std::string>().c_str();
		method = Method(methodName);
	} else {
		std::cout << "Missing method!\n";
		return 1;
	}
	/*	if (method == Method::AND_OR) {
	 isAndOr = true;
	 }
	 if (method == Method::OR) {
	 isAndOr = false;
	 }*/
	// heuristics
	if (vm.count("heur")) {
		heurName = vm["heur"].as<std::string>().c_str();
		heur = Heur(heurName);
	} else {
		heur = Heur::BOTH; // default
	}
	if (heur == Heur::UB) {
		std::cout << "Only upper (bound) heuristics applied!\n";
		isLbHeur = false;
	} else {
//		assert(heur == Heur::BOTH);
		std::cout << "Both upper and Lower (bound) heuristics applied!\n";
		isLbHeur = true;
	}
	// file to output pseudotree
	if (vm.count("pst_file")) {
		pstName = vm["pst_file"].as<std::string>().c_str();
	} else {
		pstName = "";
	}
	// priority type
	if (vm.count("priority")) {
		priorityName = vm["priority"].as<std::string>().c_str();
	} else {
		priorityName = "";
	}

	if (vm.count("searchmode")) {
		searchModeName = vm["searchmode"].as<std::string>().c_str();
		searchMode = Mode(searchModeName);
	} else {
//		std::cout << "Default search mode (SMA) in wmbsearch applied!\n";
		searchMode = Mode::SMA;
	}

	if (isAllMemControl > 0) {
		wmbMemLimit = memoryBudget;
	}

	double processEnteringTime = timeSystem();

	/*** READ IN PROBLEM FILE **********************************************************/
	std::cout << "Reading model file: " << probName << "\n";
	ifstream is;
	is.open(probName);
	if (!is.is_open())
		throw std::runtime_error("Failed to open problem file");
	mex::vector<Factor> forig = Factor::readUai10(is);
	if (addEpsilon > 0.0)
		for (size_t f = 0; f < forig.size(); ++f)
			forig[f] += addEpsilon;
	//
	mex::wmbe wmbUB(forig);
//	mex::wmbe wmbLB(forig);
	mex::wmbe wmbLB;
	if (isLbHeur) {
		wmbLB = mex::wmbe(forig);
	}

	size_t nvar = wmbUB.nvar();
	bel.resize(nvar);
	mex::vector<size_t> evid(nvar, 0);
	std::cout << "Graphical model over " << wmbUB.nvar() << " variables with "
			<< forig.size() << " factors\n";

	/*** READ IN EVIDENCE FILE *********************************************************/
	VarSet evVar;
	ifstream is2;
	if (vm.count("evidence")) {
		is2.open(vm["evidence"].as<std::string>().c_str());
	}
	if (is2.is_open()) {
		std::cout << "Got evidence file\n";
		//int nEvid; is2 >> nEvid;
		int nEvid = 1;  // 2014 format: single evidence, no # of evidences entry
		std::cout << "Got " << nEvid << " evidences?\n";
		if (nEvid > 0) {
			int nEvidVar;
			is2 >> nEvidVar;
			for (size_t i = 0; i < nEvidVar; i++) {
				uint32_t vid;
				size_t vval;
				is2 >> vid >> vval;
				evid[vid] = vval;
				evVar |= wmbUB.var(vid);
				bel[vid] = Factor::delta(wmbUB.var(vid), vval);
			}
			std::cout << "Evidence on variables " << evVar << "\n";
			wmbUB.condition(evVar, evid);
			if (isLbHeur) {
				wmbLB.condition(evVar, evid);
			}
		}
	} else
		std::cout << "Evidence file not specified or not found\n";

	/*** Read in query (MAX) variables ************************************************************/
	VarSet maxVars;
	ifstream ismap;
//	auto queryFile = vm["query-file"].as<std::string>().c_str();
	if (vm.count("query-file")) {
//		ismap.open(queryFile);
		ismap.open(vm["query-file"].as<std::string>().c_str());
	}
	if (ismap.is_open()) {
		auto queryFile = vm["query-file"].as<std::string>().c_str();
		std::cout << "Reading query file: " << queryFile << std::endl;
//		std::cout << "Got query file" << std::endl;
		// the first entry is # of query (MAX) variables
		int nMax = 0;
		ismap >> nMax;
		std::cout << nMax << " MAX variables specified!" << std::endl;
		if (nMax > 0) {
			for (int i = 0; i < nMax; i++) {
				uint32_t vid;
				ismap >> vid;
				maxVars |= wmbUB.var(vid);
			}
			// Set MAX variables
			// this should not automatically set weights, thus safe for mmapIS
			wmbUB.setMaxVars(maxVars);
		}
		ismap.close();
	} else {
		std::cout << "Query (MAP) file not specified or not found" << std::endl;
	}

	/*** PREPARE REQUESTED TASK ********************************************************/
	std::cout << "Task is " << task << "\n";
	//if (task==Task::MPE) { std::cout<<"MPE task not supported\n"; return 1; }
	std::cout << "Method is " << method << "\n";

	/*** ELIMINATION ORDERS ************************************************************/

//    double startOrder = timeSystem();
	mex::VarOrder order;
	size_t InducedWidth, iOrder = 1;
	double score = wmbUB.order(mex::graphModel::OrderMethod::Random, order, 0,
			memUseRandom);
	const char *orderFile = NULL; // Check for pre-specified elimination order file
	if (vm.count("order-file")) {
		orderFile = vm["order-file"].as<std::string>().c_str();
	}
	ifstream orderIStream;
	if (orderFile != NULL)
		orderIStream.open(orderFile);
	if (orderIStream.is_open()) {
		// If we were given an input file with an elimination ordering, just use that
		std::cout << "Reading ordering file: " << orderFile << "\n";
		size_t ordersize;
		orderIStream >> ordersize;
		assert(ordersize == order.size());
		for (size_t i = 0; i < order.size(); ++i) {
			size_t tmp;
			orderIStream >> tmp;
			order[i] = tmp;
		};
		orderIStream.close();
		// double check, whether the order is valid  (SUM variables are eliminated first)
		if (maxVars.size() > 0) {
			bool maxStart = false;
			for (size_t i = 0; i < ordersize; ++i) {
				if (maxVars.contains(wmbUB.var(order[i]))) {
					maxStart = true;
				} else {
					if (maxStart) {
						std::cout
								<< "Invalid order: a SUM variable eliminated later than some MAX variable!"
								<< std::endl;
						return 1;
					}
				}
			}
		}
	} else {
		// current implementation of MMAP: input order file should be given
		if (task == Task::MMAP) {
			std::cout << "order file should be specified for MMAP task!"
					<< std::endl;
			return 1;
		}

		// Otherwise, calculate elimination order(s) ////////////////////////////////////
		double startOrder = timeSystem();
		iOrder = 0;
		// Try to build new orders until time or count limit reached ////////////////////
		while (iOrder < nOrders && (timeSystem() - startOrder < timeOrder)) {
			score = wmbUB.order(mex::graphModel::OrderMethod::WtMinFill, order,
					nExtra, score);
			++iOrder;
		}

		// If we were given an ordering file name but no file, write our order out
		ofstream orderOStream;
		if (orderFile != NULL)
			orderOStream.open(orderFile);
		if (orderOStream.is_open()) {
			std::cout << "Writing elimination order to " << orderFile << "\n";
			orderOStream << order.size();
			for (size_t i = 0; i < order.size(); ++i)
				orderOStream << " " << order[i];
			orderOStream << "\n";
			orderOStream.close();
		}
	}

	InducedWidth = wmbUB.inducedWidth(order);
	std::cout << "Best order of " << iOrder << " has induced width "
			<< InducedWidth << ", score " << score << "\n";
	wmbUB.setOrder(order);
	if (isLbHeur) {
		wmbLB.setOrder(order);
	}
	/*** SET IBOUND ********************************************************************/
	if (isAndOr) {
		std::cout << "non-chain structured pseudo tree used!" << std::endl;
	} else {
		std::cout << "chain structured pseudo tree used!" << std::endl;
	}
	double mem = 0.0;
	/*	if (isAllMemControl > 0) {
	 iBound = 50;
	 }*/
	if (iBound >= InducedWidth) {
		iBound = InducedWidth;
	}

	wmbUB.setIBound(iBound);
	if (isLbHeur) {
		wmbLB.setIBound(iBound);
	}
	for (;; --iBound) {
		wmbUB.init(isAndOr);
		if (isLbHeur) {
			wmbLB.init(isAndOr);
		}
		wmbUB.setIBound(iBound);
		if (isLbHeur) {
			wmbLB.setIBound(iBound);
		}
		wmbUB.build();
		if (isLbHeur) {
			wmbLB.build();
		}
		if (heur == Heur::BOTH) {
			mem = wmbUB.memory() + wmbLB.memory();
		} else {
			mem = wmbUB.memory();
		}
		if (mem <= wmbMemLimit || iBound == 1)
			break;
		if (verbose > 1) {
			std::cout << "Projected memory too high (" << iBound << "=> " << mem
					<< "MB)\n";
		}
	}
	std::cout << "The actual ibound used is: " << iBound << " => " << mem
			<< " MB" << std::endl;
	isSolved = iBound >= InducedWidth;
	if (isAllMemControl > 0) {
		searchMemLimit = memoryBudget - mem;
	}

	wmbUB.setTheta();
	if (method == Method::MMAPIS) {
		wmbUB.setMsgFwdRep(sumVarCopies);
	}

	if (isLbHeur) {
		wmbLB.setTheta();
		if (method == Method::MMAPIS)
			wmbLB.setMsgFwdRep(sumVarCopies);
	}
	switch (task) {
	case Task::PR:
		wmbUB.setElimType(mex::wmbe::ElimType::SumUpper);
		if (isLbHeur) {
			wmbLB.setElimType(mex::wmbe::ElimType::SumLower);
		}
		break;
	case Task::MMAP:
//		wmbUB.setElimType(mex::wmbe::ElimType::MaxSumUpper);
		// set elimination type for each variable
		for (size_t v = 0; v < wmbUB.nvar(); ++v) {
			auto var = wmbUB.var(v);
//			if (maxVars.contains(var))
			// note that for MMAPIS, weights for MAX variables still sums to 1.
			if ((method != Method::MMAPIS) && maxVars.contains(var)) {
				// MAX variable
				wmbUB.setElimType(var, mex::wmbe::ElimType::MaxUpper);
			} else {
				// SUM variable
				wmbUB.setElimType(var, mex::wmbe::ElimType::SumUpper);
			}
		}
		break;
	default:
		std::cout << "Task not supported\n";
		return 1;
	}

	//wmb.dump(std::cout);

	/*** RUN ITERATIONS ****************************************************************/

	double startIter = timeSystem();
	double UB = 0.0, LB = 0.0;
	size_t it = 0;
	assert(it < stopIter);

	//Caveat: we cannot do backward message passing in my current implementation
	if (method == Method::MMAPIS)
		stopIter = 1;

	while (true) {
		++it;
//	for (size_t it = 0; it < stopIter && (timeSystem() - startIter < stopTime); ++it)
		UB = 0.0;
		LB = 0.0;
		for (size_t b = 0; b < wmbUB.getOrder().size(); ++b) {
			if (method == Method::MMAPIS)
				UB += wmbUB.msgForwardAug(b, dampTheta, stepWeights,
						sumVarCopies);
			else
				UB += wmbUB.msgForward(b, dampTheta, stepWeights);
		}
		if (isLbHeur) {
			for (size_t b = 0; b < wmbLB.getOrder().size(); ++b) {
				if (method == Method::MMAPIS)
					LB += wmbLB.msgForwardAug(b, dampTheta, stepWeights,
							sumVarCopies);
				else
					LB += wmbLB.msgForward(b, dampTheta, stepWeights);
			}
		}

		//if (isAdaptive && obj > last) { dampTheta /= 2; stepWeights /= 2; }
		//last = obj;
		std::cout.precision(10);
		std::cout << "[" << std::setw(10) << timeSystem() - startIter << "] : "
				<< std::setw(10) << UB << " > ln Z > " << std::setw(10) << LB
				<< std::endl;
		// make sure we end up with forward message passing
		if (it == stopIter || timeSystem() - startIter > stopTime) {
//			std::cout << "Done " << it << " iterations of forward message passing!" << std::endl;
			break;
		}
		/*		if (timeSystem() - startIter > stopTime)
		 break;*/

		// I have not implemented "msgBackward" for "mmapIS"
		for (int b = wmbUB.getOrder().size() - 1; b >= 0; --b)
			wmbUB.msgBackward(b, dampTheta, stepWeights);
		if (isLbHeur) {
			for (int b = wmbLB.getOrder().size() - 1; b >= 0; --b)
				wmbLB.msgBackward(b, dampTheta, stepWeights);
		}
	}

	//wmbUB.dump(std::cout);
	//
	std::cout << "memory limit for WMB is: " << wmbMemLimit
			<< " (MB), actual used: " << mem << " (MB)!" << std::endl;
	startIter = timeSystem(); // just for search

	/*** Run Algorithms ****************************************************************/

	// task
	if (task == Task::MMAP) {
		searchTime = timeBudget - (timeSystem() - processEnteringTime);
		if (method == Method::MMAPIS) {
			std::cout << "[" << std::setw(10)
					<< timeSystem() - processEnteringTime
					<< "] : Run into MMAPIS..." << std::endl;
//			mmapIS(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt,
//					int nv, bool isSolved, int K, const mex::VarSet& maxVars, int nmax);
			mmapIS mmapisObj(wmbUB, UB, DELTA, verbose, searchTime,
					searchMemLimit, wmbUB.nvar(), isSolved, sumVarCopies,
					maxVars, maxVars.size());
			mmapisObj.start(NSP, NND);
		} else {
			// default: UBFS
			mmap mmapObj(wmbUB, verbose, wmbUB.nvar(), wmbUB.nMaxVar(), UB,
					searchMemLimit, searchTime);
			mmapObj.start();
		}
		return 0;
	}
	//

//	if ((method == Method::AND_OR) || (method == Method::OR))
	if (method == Method::AOBFS) {
		if (isAllTimeControl > 0) {
			searchTime = timeBudget - (timeSystem() - processEnteringTime);
			if (searchTime < 0) {
				std::cout << "no time for search!" << std::endl;
				return -1;
			}
		}
		if (searchMemLimit < 0) {
			std::cout << "no memory for search!" << std::endl;
			return -1;
		}
		if (isAndOr) {
			std::cout << "[" << std::setw(10)
					<< timeSystem() - processEnteringTime
					<< "] : Entering AND/OR search..." << std::endl;
		} else {
			std::cout << "[" << std::setw(10)
					<< timeSystem() - processEnteringTime
					<< "] : Entering OR search..." << std::endl;
		}
		wmbsearch andOrSearch(wmbLB, wmbUB, verbose, priorityName, startIter,
				timeUnit, timeRatio, LB, UB, wmbUB.nvar(), searchMemLimit,
				searchTime);
//		andOrSearch.start(searchTime, nSearch, pstName,
//				addEpsilon, probName, searchMemLimit);
		andOrSearch.startSearch(searchMode);
	} else if (method == Method::MCTS) {
		// mcts(const mex::wmbe& wmbeub, const mex::wmbe& wmbelb, int ns, double DELTA, double ub, double lb)
		std::cout << "run into mcts..." << endl;
		mcts mctsOR(wmbUB, wmbLB, nSample, DELTA, UB, LB);
//		mctsOR.start();
//		mctsOR.startTwoStageSampling();
		if (NSP < 0) {
//			mctsOR.startTwoStageSamplingWithFixedTree(searchTime, searchMemLimit, timeBudget);
			mctsOR.startTwoStageSamplingWithFixedTree(searchTime, memoryBudget,
					timeBudget);
		} else
			mctsOR.startTwoStageSamplingWithFixedTree(searchTime,
					searchMemLimit, timeBudget, NSP);
//			mctsOR.startGBFSWithNewGap(searchTime, searchMemLimit, timeBudget, NSP); // temporary use
	} else if (method == Method::SESA) {
		std::cout << "[" << std::setw(10) << timeSystem() - processEnteringTime
				<< "] : Run into searchsample..." << std::endl;
		searchsample SS(wmbUB, UB, nSample, DELTA, verbose,
				timeBudget - (timeSystem() - processEnteringTime),
				searchMemLimit, wmbUB.nvar(), isSolved);
		/*		if (nSearch > 0) {
		 SS.start(nSearch);
		 } else {
		 SS.start();
		 }*/
		SS.start(NSP, NND);
	} else if (method == Method::RDHR) {
		std::cout << "[" << std::setw(10) << timeSystem() - processEnteringTime
				<< "] : Run into RDHR..." << std::endl;
		//randheur(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt, int nv, bool isSolved, int ksp)
		if (lookahead_depth > nvar || lookahead_depth < 0) {
			// 0 for adaptive lookahead depth
			lookahead_depth = nvar;
		}
		randheur SS(wmbUB, UB, DELTA, verbose,
				timeBudget - (timeSystem() - processEnteringTime),
				searchMemLimit, wmbUB.nvar(), isSolved, KSP, lookahead_depth);
		SS.start(NSP, NND);
	} else {
		std::cout << "Method not supported" << endl;
		return 1;
	}
	//
	return 0;
}
// EOF


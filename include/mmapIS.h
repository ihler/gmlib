/*
 * mmapIS.h
 *
 *  Created on: Dec 1, 2017
 *      Author: qlou
*/

#ifndef MMAPIS_H_
#define MMAPIS_H_

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <set>

#include "mxUtil.h"
#include "graphmodel.h"
#include "wmbe.h"
#include "vector.h"
#include "../external/tree.hh"
#include "enum.h"
// namespace mex {

//MEX_ENUM(NODE_TYPE , MAX_NODE,SUM_NODE );
//MEX_ENUM(VAR_TYPE , MAX_VAR,SUM_VAR );
//typedef bool nodeType;
//typedef bool varType;

/*const bool OR_NODE = false;
const bool AND_NODE = true;
const bool SUM_VAR = false;
const bool MAX_VAR = true;*/

class mmapIS {
public:
	using nodeType = bool;
	using varType = bool;
	mmapIS(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt,
			int nv, bool isSolved, int K, const mex::VarSet& maxVars, int nmax);
	virtual ~mmapIS();
	// any default value in node is for AND node, MAX variable
	struct node {
		mmapIS::nodeType nodeType; // AND_NODE (default) or OR_NODE
		mmapIS::varType varType; // MAX_VAR (default) or SUM_VAR
		mex::Var X; // variable with id and state set
		int id = -1;  // unique id for the variable, "-1" to be safe
		int val = -1; // label, only meaningful for and node, "-1" to be safe
		/* Caveat: nsp is only useful for frontier MAX nodes and SUM nodes with a MAX parent
		 * for a frontier MAX node: times of its current best MAP sub-configuration below to be sampled
		 * for a SUM node with a MAX parent: no. of samples that pass this node
		 */
		unsigned nsp = 0;
		// initial deterministic bounds for AND and OR nodes.
		// for OR nodes, just summation of bounds of AND children.
		double ub = std::numeric_limits<double>::infinity();
//		double lb = -std::numeric_limits<double>::infinity(); // TODO: should be removed
//		double logEx = 0.0; // empirical mean
//		double logEx2 = 0.0;
		// estimate of the best MAP configuration value below this node
		// mainly used for tip MAP nodes
		double wSum = -std::numeric_limits<double>::infinity(); // \sum_i  Zhat_i/U_i
//		double sumInvUB = 0.0; // \sum_i 1/U_i = N/HM(U)
		// thus, the estimator is "wSum - sumInvUB" (in log)
		tree<node>::iterator toChild = NULL; // to the most promising child
		double down_ub = -std::numeric_limits<double>::infinity(); // used to track the local priority of the best frontier descendant
		double exact_val = -std::numeric_limits<double>::infinity(); // caveat: for "or ", should be initialized to -inf
		// for MAX nodes, if the optimal sub-configuration below is known; for SUM nodes, if
		bool solved = false;
		//  whether the value of the summation problem below is exact. it is equivalent to "solved" for SUM nodes, but not for MAX nodes.
		// Caveat: when "solved" = true for a MAX node, we will set "sumSolved" = true for the implementation convenience in sampling
		bool sumSolved = false;
		 // a configuration of its current best MAP descendants. <id, val> pairs
		// pair is ordered, i.e.,
		std::list<std::pair<int, int> > mapDesc;
		/*
		 * downEst:
		 * To track the estimate of the subproblem below this node
		 * Used for two type of nodes:
		 *  1) frontier AND-MAX node (i.e., AND-MAX node without descendants in the current search tree)
		 *  	it is the estimated value of its current best sub-MAP configuration.
		 *  2) OR-SUM node with a MAX parent.
		 *  	it is the current estimate of the SUM subproblem below (weighted mean: HM(U)/N\sum_i Z_i/U_i)
		 */
 		double downEst = -std::numeric_limits<double>::infinity();
		// track the size of the MAP space after pruning
		double mapSpaceLeft = -std::numeric_limits<double>::infinity(); // -inf for a SUM node, 0.0 for a leaf (or SOlVED) MAX node
	};
	// helper functions
	double logsumexp(std::vector<double>&);
	double logsumexp(std::list<double>&);
	double logsumexp(std::initializer_list<double>);
/*
 * functions
 */
	void setRootConfig(); // setup for root
	void setExactHeur(); // set exactHeur;
	void buildPseudoTree();
	void setMapSpaceSize();
//	void setDepth();
	void printNode(const node&);
	void checkSolved(node&);
	tree<node>::iterator pruneSolved(tree<node>::iterator);
	void updateBounds(tree<node>::iterator);
	// merge bound update and pruning into one function
	tree<node>::iterator updateBoundsPlusPruneSolved(tree<node>::iterator);
	void propMapSpace(tree<node>::iterator);
	void findChild(tree<node>::iterator); // find the best child
	int expandBest();
	tree<node>::iterator expandOneFrontierNode(tree<node>::iterator, mex::vector<uint32_t>& config);
	int backwardUpdate(tree<node>::iterator ptr);
//	void runSearch(const unsigned treeSizeLimit); // wmbsearch to create a search tree
	int runSearch(const unsigned nrnd, const long unsigned treeSizeLimit); // run several rounds of best-first search
	double logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done); // compute logF(x) over a (possibly partial) configuration of a set of variables
	tree<node>::iterator upperSampling(tree<node>::iterator); // sample from categorical distribution based on upper bounds
	std::list<mex::Var> getDescList(int id); // get descendants of a variable in the pseudo tree
	void getDesc(int id, std::list<mex::Var>& unDoneMax, std::list<mex::Var>& unDoneSum, bool maxDescOnly = false);
	std::pair<int,int> getDesc(int ID);
//	void addSample(double logFx, const std::vector<double>&); // add sample to the root
	void addSampleToRoot(double wt);
	// add one (or K) estistimator
	double addSampleToNode(tree<node>::iterator ptr, const mex::vector<uint32_t>& config, double est);
	void addMapDesc(tree<node>::iterator ptr, const mex::vector<uint32_t>& config);
	void updateGlobalEst();
//	std::pair<double, double> calcRootEBB(); // compute EBB for the root
//	std::pair<double, double> calcRootEBB(double ub); // ebb based on given ub
	// EBB based on the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
	std::pair<double, double> calcNormalizedEBB();
	bool calcMapEBB(double& ubZaug, double& lbZaug, double& lbMap); // an extension of "calcNormalizedEBB"
//	std::pair<double, double> calcHoeffding(); // Hoeffding's bounds at the root for convex combination of samples
	// update the weighted estimate: HM(U^2)*(\sum_i z_i/u_i^2)/n, U^2 = (u_1^2, ..., u_n^2), HM(U^2) is harmonic mean of U^2
	// hoeffding's bound is applicable to this estimate
//	void updateConvEstimate(double wt, double ub);
	// update the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
	// empirical bernstein bound is applicable to this estimate
	void updateNormalizedEstimate(double wt, double ub);
	void twoStepSampling();
	double twoStepSum(mex::vector<uint32_t>& config, tree<node>::iterator ptr); // two-step sampling for SUM variables by conditioning on their MAP ancestors
//	void doSampling();
	int doSampling(int, bool isInitUb = true); // whether we use the initial UB or updated UB to calculate EBB
//	void start();
//	void start(unsigned treesize); // for debug
	void start(const int nsp, const int nnd);
	//
	int display(); // display results on screen
	tree<node>::iterator DFS(tree<node>::iterator); // depth-first search
	int runDFS(const unsigned nnd); // expand a given number of frontier nodes via DFS

/*
 * data
 */
private:
	const int TIMEOUT = -2;
	const int MEMOUT = -1;
	const int SOLVED = 1;
	const bool OR_NODE = false;
	const bool AND_NODE = true;
	const bool SUM_VAR = false;
	const bool MAX_VAR = true; // true for max, do not switch!!!
	const double solvedThresh = 1e-10;
	const double epsilon = 1e-10;
	const bool DEBUG = false; // whether debug

	const double timeBudget;
	const double memoryBudget; // should be double
	const int _K; // K: no. of copies of SUM variables
	const double _logK; // log(K)

	mex::wmbe& wmbUB;
//	mex::wmbe& wmbLB;
	mex::VarOrder order; // v = order[i] is the ith variable eliminated
	mex::vector<size_t> priority; // i = priority[v] is the step in which variable v is eliminated
//	std::vector<int> height; // store height of each variable on the pseudo tree, useful for pruning
//	int iBound; // i-Bound
	const int nvar; // no. of variables
//	double nSample; // no. of samples
	double DELTA; // one-side 1-confidence, positive
	tree<node> exploredTree; // explored AND/OR subtree
	unsigned long GSIZE; // size of stored tree
	double UB, LB; // deterministic bounds, may not change
	int verbose;
	unsigned long NSP; // current no. of samples
	double startProcess;
	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
	std::vector<std::list<mex::graphModel::vindex> > ascList; // holds children for each variable
	std::vector<bool> exactHeur; // whether heuristics for variables are exact; initialized to "true"
	tree<node>::iterator root; // root = exploredTree.begin()
	mex::vector<uint32_t> tuple; // re-usable vector to accelerate, make sure reset it before use
	const double _initUB; // initial upper bound
	const bool _isExactHeur;
	// HM(U^2)/n*(\sum_i z_i/u_i^2), U^2 = (u_1^2, ..., u_n^2), HM(U^2) is harmonic mean of U^2
/*	double convEx = 0.0; //  convex combination of importance weights (in log)
	double sumSqrUB = -std::numeric_limits<double>::infinity(); // \sum_i u_i^2 across all samples (in log)
	double sumInvSqrUB = -std::numeric_limits<double>::infinity(); // \sum_i 1/u_i^2 across all samples (in log)
	*/
	double normalizedEx = 0.0; // (\sum_i y_i)/n, where y_i = z_i/u_i, z_i is a sample and u_i is the corresponding upper bound
	double normalizedEx2 = 0.0; // (\sum y_i^2)/n, where y_i = z_i/u_i
	double sumInvUB = 0.0; // \sum_i 1/U_i
	double hmUB = 0.0; // harmonic mean of (u_1, ..., u_n)
	//
	unsigned long _outFrequency = 1; // for output frequency
//
	const int _nMax; // no. of MAX variables
	mex::vector<uint32_t> mapConfig; // predicted MAP configuration, the size is still "nvar" for convenience
	// sum over all values in "mapValEst" get the current estimated value of the predicted MAP solution
	// entry for any SUM variable with a SUM parent is 0
//	std::vector<double> mapEstVec; // each entry corresponds to a state of a variable contributes to the estimate of the current MAP configuration.
	double mapSolEst = -std::numeric_limits<double>::infinity();
	bool isMapConfigSolved; // whether the value of the predicted MAP configuration is exact or not
	unsigned long _nExpansion = 0; // total number of expansions
	unsigned long _increment = 1;
	std::vector<varType> _varTypeVec;
	// for any variable X, all its descendants in the pseudo tree are consecutively located in the list
//	std::list<mex::Var> _descList; // should not changed after initialized
//	std::vector<int> _descPosVec;  // indicate position of each variable in _descList, i.e.,  loc = _descPosVec[X.label()] gives the location in _descList
//	const std::vector<int> _depthVec;
	std::vector<bool> isLeafVar; // indicates whether a variable is a leaf in the pseudo tree
//	std::stack<tree<node>::iterator> stackDFS; // stack used for DFS
	double _mapSpaceSize; // the size of the whole MAP space (in log)
	double curMapSpaceSize; // the current size of the MAP space (in log), update after pruning MAP nodes
	double sumWeightedMapSpaceSize = 0.0; // \sum_i  |A_i|/U_i, where |A_i| is the remaining size of the MAP space (in log)
	std::vector<double> _mapSpaceSizeVec; // size of MAP space below each variable
	tree<node>::iterator dfsFrontier = NULL; // current frontier node for DFS to expand
	mex::vector<uint32_t> dfsConfig; // configuration used for DFS to expand next
	bool isMapConfigNeverSampled = true; // a flag to make sure we assign a sampled MAP configuration to "mapConfig" for the first time.
/*
 *
 */
	std::vector<mex::Var> _pseudotree; // any variable between two descendants of a variable must be its descendant
	/* pair = _locPair[X], pair.first is the location of X in _pseudotree,  "pair.second-1" is the location of X's last descendant in _pseudotree
	 * i.e., indices in (pair.first, pair.second) are all X's descendants
	 */
	std::vector<std::pair<int,int> > _locPair;
};
// } /* namespace mex */
#endif  MMAPIS_H_

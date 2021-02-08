/*
 * randheur.h
 *
 *  Created on: Jun 21, 2017
 *      Author: qlou
 *
 *  function: test stochastic heuristics for AND/OR best first search
 */

#ifndef RANDHEUR_H_
#define RANDHEUR_H_

//namespace std {

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

MEX_ENUM( PRIOR , UB,UB_LB_DIFF,UB_EST_DIFF );

class randheur {
public:
	randheur(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt, int nv,
			bool isSolved, int ksp, int depth_lookahead = std::numeric_limits<int>::infinity());
	virtual ~randheur();
	// any default value in node is for AND node
	struct node {
		bool type; // AND_NODE or OR_NODE
		mex::Var X; // variable with id and state set
		int id = -1;  // unique id for the variable, "-1" to be safe
		int val = -1; // label, only meaningful for and node, "-1" to be safe
		int nsp = 0; // no. of samples at this node
		// initial deterministic bounds for AND and OR nodes.
		// for OR nodes, just summation of bounds of AND children.
		double ub = std::numeric_limits<double>::infinity();
//		double lb = -std::numeric_limits<double>::infinity(); // TODO: should be removed
//		// updated deterministic bounds for AND and OR nodes. improve as exploredTree grows
//		double UB = std::numeric_limits<double>::infinity();
//		double LB = -std::numeric_limits<double>::infinity();
//		double exact_val = 0.0; // caveat: for "or ", should be initialized to -inf
		// sum of weights (thetas) w(n,s) on the arcs n->s on the path from the root to current node
		// w(n,s) has actual weight if n is OR node, s is an AND child; otherwise w(n,s)==0;
//		double cost = 0.0;
		double logEx = -std::numeric_limits<double>::infinity(); // empirical mean
		double logEx2 = -std::numeric_limits<double>::infinity();
		// bias term in the empirical Bernstein bound, P( x < Ex+bias ) >= 1-delta
//		double logBias = 0.0; // only for AND node
		//
		tree<node>::iterator toChild = NULL; // to the most promising child
		double down_ub = std::numeric_limits<double>::infinity(); //
		double exact_val = 0.0; // caveat: for "or ", should be initialized to -inf
		bool solved = false; // whether all the paths from the current node have been instantiated.
		// initialized to -inf is important
		double down_est = -std::numeric_limits<double>::infinity(); // a randomized "down_lb"
		double est = -std::numeric_limits<double>::infinity();// estimate of the subproblem (with exact_val incorporated), a randomized "lb"

	};
	// helper functions
	double logsumexp(std::vector<double>&);
	double logsumexp(std::initializer_list<double>);
/*
 * functions
 */
	void setRootConfig(); // setup for root
	void setExactHeur(); // set exactHeur;
	void buildPseudoTree();
	void printNode(const node&);
	bool checkSolved(const node&); // check whether a node is solved or not
	tree<node>::iterator pruneSolved(tree<node>::iterator);
	void updateBounds(tree<node>::iterator);
//	void findChild(tree<node>::iterator); // find the best child
	int expandBest();
	int backwardUpdate(tree<node>::iterator ptr);
//	void runSearch(const unsigned treeSizeLimit); // wmbsearch to create a search tree
//	int runSearch(const unsigned nrnd, const unsigned treeSizeLimit); // run several rounds of best-first search
	// compute logF(x) over a (possibly partial) configuration of a set of variables
	double logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done);
	// generalize to the case where we only sample to some level
	double logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done, const std::list<mex::Var>& tipVars);
	tree<node>::iterator upperSampling(tree<node>::iterator); // sample from categorical distribution based on upper bounds
	std::list<mex::Var> getDescList(int id); // get descendants of a variable in the pseudo tree
	void addSample(double logFx, const std::vector<double>&); // add sample to the root
	void addSampleToNode(node& n, double wt); // add sample to a given node
	std::pair<double, double> calcRootEBB(); // compute EBB for the root
	std::pair<double, double> calcRootEBB(double ub); // ebb based on given ub
	// EBB based on the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
	std::pair<double, double> calcNormalizedEBB();
	std::pair<double, double> calcHoeffding(); // Hoeffding's bounds at the root for convex combination of samples
	// update the weighted estimate: HM(U^2)*(\sum_i z_i/u_i^2)/n, U^2 = (u_1^2, ..., u_n^2), HM(U^2) is harmonic mean of U^2
	// hoeffding's bound is applicable to this estimate
	void updateConvEstimate(double wt, double ub);
	// update the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
	// empirical bernstein bound is applicable to this estimate
	void updateNormalizedEstimate(double wt, double ub);
	void twoStepSampling(); // two-step sampling as described in the nips-17 paper
//	void doSampling();
	int doSampling(int nsp);
//	void start();
//	void start(unsigned treesize); // for debug
//	void start(const int nsp, const int nnd);

	/*
	 * helper functions
	 */
//	void setPseudoTree(); // set the pseudo tree for all variables
	void setDepth(); // set depth on the pseudo tree for each variable
	void setLookaheadDepth(); // set lookahead depth for each variable
	void getDescVars(int id, std::list<mex::Var>& desc, std::list<mex::Var>& tipVars);
	void setRandHeur(node&, std::list<mex::Var>& unDone, std::list<mex::Var>& tipVars, mex::vector<uint32_t>& config, std::vector<double>& logQxCond);
	void setRandHeur(node&, mex::vector<uint32_t>& config);
	void updateRandHeur(tree<node>::iterator); // update heuristics and priorities after a global sampling step
	double getPrior(node &);
	 // get cumulative estimate from children. TODO: may cache it in the node structure?
	std::pair<double,double> getSibBounds(tree<node>::iterator);
	double getChildEstimate(tree<node>::iterator);
	std::pair<double,double> getChildBounds(tree<node>::iterator); // est and ub
	void updateEsts(tree<node>::iterator); // update its estimate based on its children, always do after one draw of samples
	void setEst(tree<node>::iterator, mex::vector<uint32_t>&); // set estimate for newly generated node
	void findChild(tree<node>::iterator, double cumUb, double cumEst); // cumUb and cumEst are analogous to topUb, topLb
	void findChild(tree<node>::iterator);
	double robustLogDiff(double large, double small);
	int display(); // display results on screen
	int runSearch(const int nnd, const unsigned treeSizeLimit);
	void start(const int nsp, const int nnd);

/*
 * data
 */
private:
	const int TIMEOUT = -2;
	const int MEMOUT = -1;
	const int SOLVED = 1;
	const bool OR_NODE = false;
	const bool AND_NODE = true;
	const double solvedThresh = 1e-10;
	const double epsilon = 1e-10; // 1e-10
	const bool DEBUG = false; // whether debug
//	const bool DEBUG = true; // whether debug

	const double timeBudget;
	const double memoryBudget; // should be double
	const int _nsp = 1; // no. of samples created for each newly generated frontier node
	const PRIOR _priorType; // search priority type
	const int nvar; // no. of variables
	const int _depth_lookahead = std::numeric_limits<int>::infinity(); // default

	mex::wmbe& wmbUB;
//	mex::wmbe& wmbLB;
	mex::VarOrder order; // v = order[i] is the ith variable eliminated
	mex::vector<size_t> priority; // i = priority[v] is the step in which variable v is eliminated
//	std::vector<int> height; // store height of each variable on the pseudo tree, useful for pruning
//	int iBound; // i-Bound

//	double nSample; // no. of samples
	double DELTA; // one-side 1-confidence, positive
//	const int searchType; // ANDOR_SEARCH or OR_SEARCH
//	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
//	std::vector<std::list<mex::graphModel::vindex>> ascList; // holds children for each variable
	tree<node> exploredTree; // explored AND/OR subtree
	unsigned int GSIZE; // size of stored tree
	double UB, LB; // global determisitc bounds
	int verbose;
	long NSP; // current no. of samples
	double startProcess;
	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
	std::vector<std::list<mex::graphModel::vindex>> ascList; // holds children for each variable
	//
	std::vector<bool> exactHeur; // whether heuristics for variables are exact; initialized to "true"
	tree<node>::iterator root; // root = exploredTree.begin()
	mex::vector<uint32_t> tuple; // re-usable vector to accelerate, make sure reset it before use
	const double initUB; // initial upper bound
	const bool _isExactHeur;
	// HM(U^2)/n*(\sum_i z_i/u_i^2), U^2 = (u_1^2, ..., u_n^2), HM(U^2) is harmonic mean of U^2
	double convEx = 0.0; //  convex combination of importance weights (in log)
	double sumSqrUB = -std::numeric_limits<double>::infinity(); // \sum_i u_i^2 across all samples (in log)
	double sumInvSqrUB = -std::numeric_limits<double>::infinity(); // \sum_i 1/u_i^2 across all samples (in log)
	double normalizedEx = 0.0; // (\sum_i y_i)/n, where y_i = z_i/u_i, z_i is a sample and u_i is the corresponding upper bound
	double normalizedEx2 = 0.0; // (\sum y_i^2)/n, where y_i = z_i/u_i
	double sumInvUB = 0.0; // \sum_i 1/u_i
	double hmUB = 0.0; // harmonic mean of (u_1, ..., u_n)
	long _outFrequency = 1; // for output frequency
	const int _ksp = 1; // no. of samples drawn for each unsolved frontier node
//	std::list<double> topCumUb; // cumulative path ub from root to frontier
//	std::list<double> topCumEst; // cumulative path estimate from root to frontier
	//
	std::stack<std::pair<tree<node>::iterator, double>> cumUbStack; // cumulative path ub from root to frontier
//	std::stack<std::pair<tree<node>::iterator, double>> cumEstStack;
	unsigned long long _nExpansion = 0; // total number of expansions
	unsigned long long _increment = 1;
	std::vector<int> _depthVec; // depth of each variable in the pseudo tree, identified by its id number
	unsigned long _nsp_root = 0; // no. of samples drawn from the root via twoStep sampling
//	tree<mex::Var> _pseudotree; // pseudo tree of variables
	std::vector<int> _adaMaxDepthVec; // max depth of descendants for each variable, used for adaptive depth
};

//} /* namespace std */

#endif /* RANDHEUR_H_ */

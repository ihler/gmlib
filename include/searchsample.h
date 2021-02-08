/*
 * searchsample.h
 *
 *  Created on: Oct 10, 2016
 *      Author: qlou
 */

#ifndef SEARCHSAMPLE_H_
#define SEARCHSAMPLE_H_

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <set>
#include <initializer_list>

#include "mxUtil.h"
#include "graphmodel.h"
#include "wmbe.h"
#include "vector.h"
#include "../external/tree.hh"
// namespace mex {
//typedef bool nodeType;
//typedef bool varType;
//
class searchsample {
public:
	searchsample(mex::wmbe& wmbeub, double ub, int ns, double delta, int vb, double tbgt, double mbgt, int nv, bool isSolved);
	virtual ~searchsample();
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
		double logEx = 0.0; // empirical mean
		double logEx2 = 0.0;
		// bias term in the empirical Bernstein bound, P( x < Ex+bias ) >= 1-delta
//		double logBias = 0.0; // only for AND node
		//
		tree<node>::iterator toChild = NULL; // to the most promising child
		double down_ub = 0.0; // topPath_ub is necessary if we only use upper heuristic
		double exact_val = 0.0; // caveat: for "or ", should be initialized to -inf
		bool solved = false; // whether all the paths from the current node have been instantiated.
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
	void findChild(tree<node>::iterator); // find the best child
	int expandBest();
	int backwardUpdate(tree<node>::iterator ptr);
	void runSearch(const unsigned treeSizeLimit); // wmbsearch to create a search tree
	int runSearch(const unsigned nrnd, const unsigned treeSizeLimit); // run several rounds of best-first search
	double logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done); // compute logF(x) over a (possibly partial) configuration of a set of variables
	tree<node>::iterator upperSampling(tree<node>::iterator); // sample from categorical distribution based on upper bounds
	std::list<mex::Var> getDescList(int id); // get descendants of a variable in the pseudo tree
	void addSample(double logFx, const std::vector<double>&); // add sample to the root
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
	void twoStageSampling();
	void doSampling();
	int doSampling(int, bool isInitUb = true); // whether we use the initial UB or updated UB to calculate EBB
	void start();
	void start(unsigned treesize); // for debug
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
	const double epsilon = 1e-10;
	const bool DEBUG = false; // whether debug

	const double timeBudget;
	const double memoryBudget; // should be double

	mex::wmbe& wmbUB;
//	mex::wmbe& wmbLB;
	mex::VarOrder order; // v = order[i] is the ith variable eliminated
	mex::vector<size_t> priority; // i = priority[v] is the step in which variable v is eliminated
//	std::vector<int> height; // store height of each variable on the pseudo tree, useful for pruning
//	int iBound; // i-Bound
	int nvar; // no. of variables
	double nSample; // no. of samples
	double DELTA; // one-side 1-confidence, positive
//	const int searchType; // ANDOR_SEARCH or OR_SEARCH
//	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
//	std::vector<std::list<mex::graphModel::vindex>> ascList; // holds children for each variable
	tree<node> exploredTree; // explored AND/OR subtree
	unsigned int GSIZE; // size of stored tree
	double UB, LB; // deterministic bounds, may not change
	int verbose;
	long NSP; // current no. of samples
	double startProcess;
	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
	std::vector<std::list<mex::graphModel::vindex> > ascList; // holds children for each variable
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
};
// } /* namespace mex */

#endif /* SEARCHSAMPLE_H_ */

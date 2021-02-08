/*
 * mcts.h
 * Monte Carlo tree search (mcts) for bounding the partition function
 *
 *  Created on: May 26, 2016
 *      Author: qlou
 */

#ifndef MCTS_H_
#define MCTS_H_

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

//#ifndef OR_NODE
//#define OR_NODE false
//#endif
//#ifndef AND_NODE
//#define AND_NODE true
//#endif

class mcts {
public:
	mcts(mex::wmbe&, mex::wmbe&, int ns, double dlt, double ub, double lb);
	virtual ~mcts();
	struct node {
		bool type; // AND_NODE or OR_NODE
		int depth = -1; // 0 : nvar-1 for normal variables
		mex::Var X; // variable with id and state set
//		int id = -1;  // unique id for the variable, "-1" to be safe
		int val = -1; // label, only meaningful for and node, "-1" to be safe
		int nsp = 0; // no. of samples at this node
		double delta = 1.0; // confidence
		// initial deterministic bounds for AND and OR nodes.
		// for OR nodes, just summation of bounds of AND children.
		double ub = std::numeric_limits<double>::infinity();
		double lb = -std::numeric_limits<double>::infinity();
		// updated deterministic bounds for AND and OR nodes. improve as exploredTree grows
		double UB = std::numeric_limits<double>::infinity();
		double LB = -std::numeric_limits<double>::infinity();
//		double exact_val = 0.0; // caveat: for "or ", should be initialized to -inf
		// sum of weights (thetas) w(n,s) on the arcs n->s on the path from the root to current node
		// w(n,s) has actual weight if n is OR node, s is an AND child; otherwise w(n,s)==0;
		double cost = 0.0;
		double logEx = 0.0; // mean estimator
		double logEx2 = 0.0;
		// bias term in the empirical Bernstein bound, P( x < Ex+bias ) >= 1-delta
		double logBias = 0.0; // only for AND node

//		// used for GBFS
//		double priority = 0.0;
//		tree<node>::iterator toPar = NULL; // pointer to parent
//		bool operator<(const node& n) const {
//			return priority < n.priority;
//		};
	};
	void setRootConfig(); // setup for root
	void buildPseudoTree();
	void sample(int num); // sample
	void sampleOne(); // create one sample and do some necessary process
	// add the configuration to the explored tree
	void addSampleToTree(const mex::vector<uint32_t>& xhat, const std::vector<double>& condUB,
			const std::vector<double>& logQxVec, const std::vector<double>& thetaVec);
	tree<node>::iterator addChild(tree<node>::iterator, mex::Var X, const mex::vector<uint32_t>& config, double wt=0.0, double cost=0.0);
	void addSampleToNode(tree<node>::iterator, double wt=0.0);
	void calcNodeBias(tree<node>::iterator);
	double setConfidence(tree<node>::iterator); // set initial confidence
	void setBounds(node&, mex::vector<uint32_t>); // set initial deterministic bounds
	std::pair<double, double> calcEBB(tree<node>::iterator); // lower and upper bounds
	double logsumexp(std::vector<double>&); // helper function
	double logsumexp(std::initializer_list<double>);
	int aggregateBounds(int max_depth);
	int aggregateBounds(int min_sample, bool);
	void start();
/*
 * two-stage sampling with a varying tree (NOT quite successful)
 */
	void updatePathBounds(tree<node>::iterator); // update deterministic bounds via back propagation from current node
	// add one sample to the path defined by this sample
	void addSampleToPath(tree<node>::iterator, const mex::vector<uint32_t>& config, const std::vector<double>& logQxVec, double logFx);
	// expandTree if some conditions are satisfied
	tree<node>::iterator expandTree(tree<node>::iterator, const mex::vector<uint32_t>& config, const std::vector<double>& logQxVec, double logFx);
	void twoStageSampling(int);
	void startTwoStageSampling();

/*
 * two-stage sampling with a fixed tree. Use GBFS to generate a tree first
 */
	// GBFS
	double calcPriority(node& n) {
		return n.ub + std::log( 1 - std::exp( std::min(n.ub,n.lb)-n.ub ) );
	}
	struct myPair {
		double priority = -std::numeric_limits<double>::infinity();
		tree<node>::iterator ptr = NULL;
		bool operator<(const myPair& p) const {
					return priority < p.priority;
		};
	};

	void runGBFS(long treeSizeLimit, double timeLimit);
	// calculate EBB for the root, since we can use the refined upper bound for the root
	std::pair<double, double> calcRootEBB();
	tree<node>::iterator expandTree(tree<node>::iterator);
	tree<node>::iterator upperBoundBasedSampling(tree<node>::iterator); // sample from categorical distribution based on upper bounds
	void twoStageSamplingWithFixedTree();
	void startTwoStageSamplingWithFixedTree(double searchTimeLimit, double memBudget, double timeBudget);
	void printTree(); // print all the information of the explored tree. only used for debug

/*
 * two-stage sampling with a fixed tree. Use WMB-IS to generate a tree first
*/
	void genTreeBySampling(long treeSizeLimit, double timeLimit, int NSP); // generate a tree by importance sampling
	void startTwoStageSamplingWithFixedTree(double searchTimeLimit, double searchMemLimit, double timeBudget, int NSP);
/*
 * try a simple idea: use the difference between upper and estimated value as priority
 */
	tree<node>::iterator expandTree(tree<node>::iterator, mex::vector<uint32_t>& tuple);
	void runGBFSWithNewGap(long treeSizeLimit, double timeLimit, int NSP);
	void startGBFSWithNewGap(double searchTimeLimit, double searchMemLimit, double timeBudget, int NSP);
/*
 * data
 */
private:
	const bool OR_NODE = false;
	const bool AND_NODE = true;

	mex::wmbe& wmbUB;
	mex::wmbe& wmbLB;
	mex::VarOrder order; // v = order[i] is the ith variable eliminated
	mex::vector<size_t> priority; // i = priority[v] is the step in which variable v is eliminated
	int nvar; // no. of variables
	double nSample; // no. of samples
	double DELTA; // one-side confidence, positive
//	const int searchType; // ANDOR_SEARCH or OR_SEARCH
//	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
//	std::vector<std::list<mex::graphModel::vindex>> ascList; // holds children for each variable
	tree<node> exploredTree; // explored AND/OR subtree
	unsigned int GSIZE; // size of stored tree
	double UB, LB; // deterministic bounds, may not change
	int cur_nSample; // current no. of samples
	int verbose;
};

#endif /* MCTS_H_ */

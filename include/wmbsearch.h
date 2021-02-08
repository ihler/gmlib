/*
 * wmbsearch.h
 *
 *  Created on: Aug 31, 2015
 *      Author: qlou
 */

#ifndef WMBSEARCH_H_
#define WMBSEARCH_H_

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

//#ifndef OR_NODE
//#define OR_NODE false
//#endif
//#ifndef AND_NODE
//#define AND_NODE true
//#endif

//#define ANDOR_SEARCH 0
//#define OR_SEARCH 1
// it's possible that the graph is not connected even if the input is.
// thus, multiple roots can appear.
// be careful about iterator invalidation problem
// turn on the AndOr pseudo tree flag when initializing wmbe class
MEX_ENUM( Mode , DFS,SMA );

class wmbsearch {
public:
	wmbsearch(mex::wmbe& wmblb, mex::wmbe& wmbub, int v, std::string p, double st, double t, double r, double lb, double ub, size_t nv, double mb, double tb);
	virtual ~wmbsearch();
	// all values related to bounds are in log scale
	struct searchnode {
		bool type; // AND_NODE, OR_NODE
		mex::Var X; // variable with id and state set
		int id = -1;  // unique id for the variable, "-1" to be safe
		int val = -1; // label, only meaningful for "and" node, "-1" to be safe
//		int depth = 0; // depth in the pseudo tree, only consider "and" node
		bool solved = false; // whether all the paths from the current node have been instantiated.
//		whether children have been removed due to memory limit
		bool isReExpand = false; //
		// ub, lb: bounds for subtree rooted at the node,
		// pre_ub, pre_lb: help update ub, lb. not necessarily needed.
		// we should have path bounds from two directions
		double ub=0.0, lb=0.0; // best so far

		// it should be tied to current status
		double exact_val = 0.0; // caveat: for "or ", should be initialized to -inf
//		double tub=0.0, tlb=0.0; // these two should always be the same
//		double pre_ub=0.0, pre_lb=0.0;
//		double cur_ub = 0.0, cur_lb = 0.0; // current value of upper bounds, may not be the best
		// topPath_ub, topPath_lb are pathUB, pathLB from the root to current node
		// (for "or" nodes, including bounds from current node's siblings)
		// downPath_ub, downPath_lb are pathUB, pathLB from the leaf to current node
		// (for "or" nodes, not including bounds from current node's siblings)
		// (for "or" nodes, not including leaf_ub, leaf_lb)
		// topPath_ub, topPath_lb, downPath_ub, downPath_lb do not contain tub, tlb from the current "and" node
		double topPath_ub = 0.0, topPath_lb = 0.0, downPath_ub = 0.0, downPath_lb = 0.0;
//		double leaf_ub = 0.0, leaf_lb = 0.0; // bounds for the most promising descendant in OPEN
		// exact value for the subproblem. may only contain partial results if some child has not been solved
		// for pruning's numerical stability
		double downWorstPath_ub = 0.0, downWorstPath_lb = 0.0;
		// caveat: assuming append_child does not result in invalidation of iteration for the tree structure
		tree<searchnode>::iterator toChild = NULL; // to the most promising child
		tree<searchnode>::iterator toWorst = NULL; // to the worst child
		double priority = std::numeric_limits<double>::infinity(); // the best priority of your descendants you discovered before
	};
	/*struct searchNode {
		double priority;
		size_t depth;
		mex::vector<uint32_t> tuple;
		double ub, lb, costUB, costLB;
		bool operator< (const searchNode& n) const { return priority < n.priority; }
	};
	double calcPriority(searchNode n) {
		//return n.depth;
		//return (n.ub-n.lb)*std::log((double)n.depth);
		return n.ub + std::log( 1 - std::exp( std::min(n.ub,n.lb)-n.ub ) );
	}*/

//	void calcPriority(searchnode& n, const std::string priorType);
//	void calcPriority(searchnode& n, const std::string priorType, const mex::wmbe& wmbUB, const mex::wmbe& wmbLB);
	void setRootConfig(); // setup for root
	// update bounds for nodes along the current path, must run after changing bounds of a leaf node
//	void updateBounds(tree<searchnode>::post_order_iterator);
	void buildPseudoTree();
//	bool runSearch(tree<searchnode>::iterator toTop, const mex::wmbe& wmbUB, const mex::wmbe& wmbLB);
//	void lookAhead(searchnode& n, mex::vector<uint32_t>& tuple, const mex::wmbe& wmbUB, const mex::wmbe& wmbLB);
//	std::vector<double> getOracle(const mex::wmbe& wmbUB, const mex::wmbe& wmbLB, int depth=-1); // default, no computation
	void start(double, int, std::string,
			double addEpsilon, std::string probName, const double memLimit);
	void writeNodeToFile(unsigned, std::ostringstream&);
	void writePseudoTreeToFile(std::string);
	double calcPriority(const tree<searchnode>::iterator); // evaluate the (global or local) priority to determine which child is most promising
	tree<searchnode>::iterator forwardPass(); // forward pass to find the node in OPEN with the highest priority
//	int backwardPass(tree<searchnode>::iterator); // backward pass to update bounds and other information, (out of date)
	int debugMsg(const tree<searchnode>::iterator); // output information for a node, for debug purpose
	void printTree(); // print all the information of the explored tree. only used for debug
	void printNode(const searchnode& n);
	void printPath(tree<searchnode>::iterator); // print all the information from the root to the current node;
	double exploredRatio(); // the ratio: size of the explored tree / size of the overall search space
	bool isPruned(tree<searchnode>::iterator, bool isLeaf = false); // check whether pruning is needed
	tree<searchnode>::iterator DFS(const double); // when almost reach the memory budget, do DFS. The memory usage only proportional to the max depth
	void propagateBounds(tree<searchnode>::iterator); // propagate bounds upward for DFS
	tree<searchnode>::iterator doPruning(tree<searchnode>::iterator, tree<searchnode>::iterator);
	std::vector<size_t> sortIndexes(const std::vector<double> &); // return sorted indexes, ascending order
	bool checkSolved(searchnode &); // check whether a node is solved or not
	double logsumexp(std::vector<double>&); //
	int bottomUpPass(tree<searchnode>::iterator);
	void setDownPath(tree<searchnode>::iterator ptr, tree<searchnode>::iterator sib);
//	void SMA(); // simplified memory-bounded A*
	int expandBestNode();
	int backwardUpdate(tree<searchnode>::iterator);
	std::pair<double,double> getSibBounds(tree<searchnode>::iterator ptr);
	std::pair<double,double> getChildBounds(tree<searchnode>::iterator ptr);
	void updateBounds(tree<searchnode>::iterator);
	double calcPriority(const searchnode&);
	double calcWorstPriority(const searchnode&); // priority of the worst descendant in OPEN
	double calcFrontierPrior(const searchnode&); // calculate the node as if it is a frontier node
	tree<searchnode>::iterator prune(tree<searchnode>::iterator ptr);
	int removeWorstNode();
	int DFS_Beta(double startTime, double& outFrequency); // upgraded version of DFS
	void startSearch(Mode mode = Mode::SMA);

private:
	const bool OR_NODE = false;
	const bool AND_NODE = true;
	//
//	const int DFS_MODE = 0;
//	const int SMA_MODE = 1;
	//
	mex::wmbe& wmbUB;
	mex::wmbe& wmbLB;
//	mex::VarOrder order; // v = order[i] is the ith variable eliminated
	const size_t nvar;
	// a const thresh between its upper and lower bounds to determine whether we can call a variable "solved"
	const double solvedThresh = 1e-10;
	// log bound ratio, decrease by factor of 10
	double boundGap = 100;

	double LB, UB; // current upper and lower bounds for log partition function
//	int MAX_NODES; // max number of iterations;
//	int NUM_NODES; // no. of overall nodes that have been visited, i.e., the size of the explored tree

	const int verbose; // 0: no verbose; 1: regular output; 2: output in geometric growth; 3: detailed output. mainly for debug
	std::string priorityName; // type of priority

	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
	std::vector<std::list<mex::graphModel::vindex>> ascList; // holds children for each variable
	tree<searchnode> exploredTree; // explored And OR subtree
	unsigned GSIZE = 0; // size of the exploredTree, including all "and", "or" nodes

	const double startIter;
	double timeThresh = 0.0;
	// output frequency
	const double timeUnit;
	double timeRatio; // will be changed

	const double timeBudget;
	const double memoryBudget;
};

#endif  WMBSEARCH_H_

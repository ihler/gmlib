/*
 * mmap.h
 *
 *  Created on: Oct 28, 2016
 *      Author: qlou
 */

#ifndef MMAP_H_
#define MMAP_H_

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

typedef bool nodeType;
typedef bool varType;

const bool OR_NODE = false;
const bool AND_NODE = true;
const bool SUM_VAR = false;
const bool MAX_VAR = true;

class mmap {
public:
	mmap(mex::wmbe& wmbub, int vbs, int nv, int nmax, double ub, double mb, double tb);
	virtual ~mmap();
	struct node {
		nodeType nodetype = AND_NODE; // AND_NODE or OR_NODE
		varType vartype = MAX_VAR; // MAX_VAR or SUM_VAR, variable type
//		int depth = -1; // 0 : nvar-1 for normal variables
		mex::Var X; // variable with id and state set
		int id = -1;  // unique id for the variable, "-1" to be safe
		int val = -1; // label, only meaningful for and node, "-1" to be safe
		//		double lb = -std::numeric_limits<double>::infinity();
		double ub = std::numeric_limits<double>::infinity();
		double exact_val = 0.0; // caveat: for OR_NODE, should be set to -inf, for AND_NODE, it is heuristicTheta
		// initialization to inf is important
		double down_ub = std::numeric_limits<double>::infinity(); // for the best descendant ever, important to be initialized to be INF
//		double downMax_ub = std::numeric_limits<double>::infinity();
//		double downWorst_ub = std::numeric_limits<double>::infinity(); // TODO we may remove it and only use downWorstMax_ub to save memory
		double downWorstMax_ub = std::numeric_limits<double>::infinity(); // for the worst child on the fringe, max branch bound
		bool solved = false; // whether solved or not. a node is solved only when all its children are solved.
		bool expanded = false; // whether expanded before or not
		varType frontierType = MAX_VAR; // type of the best frontier descendant
		varType worstFrontierType = MAX_VAR;
		// Caveat: the following two quantity being "NULL" does not imply a node has no children
		tree<node>::iterator toChild = NULL; // to the most promising child
		tree<node>::iterator toWorst = NULL; // to the least promising child
	};
	void buildPseudoTree();
	void setRootConfig(); // setup for root
	void setExactHeur(); // set exactHeur;
	void printNode(const node&);
	void printTree(); // print the entire explored tree
	void printRemovableSet(); // print all removable nodes
	void printOPEN(); // print the open set
	void printMaxConfig();
	void getMaxConfig(); // get the max configuration based on current explored tree
	double logsumexp(std::list<double>& vec);
	double getSibUB(tree<node>::iterator);
	bool checkSolved(const node&); // check whether a node is solved or not
	tree<node>::iterator pruneSolved(tree<node>::iterator ptr, tree<node>::iterator ancestor); // stop at the ancestor
	void updateBounds(tree<node>::iterator);
	void findChild(tree<node>::iterator); // find best and worst child
	int expandBest();
	int removeWorst();
	int backwardUpdate(tree<node>::iterator);
	void greedyCompletion(); // greedy approach to complete the current (possibly partial) best MAX configuration
	void printGreedyCompletion(); // printing function of greedyCompletion
	void start();
	// debug functions

private:
	//
	const double epsilon = 1e-10;
	mex::wmbe& wmbUB;
	const int verbose;
	const int nvar; // no. of variables
	const int nMaxVar; // no. of max variables
	double UB; // upper bound on the current best (sub-)solution tree
	const double memoryBudget;
	const double timeBudget;
	std::list<mex::graphModel::vindex> rts; // nodes without parent, possibly multiple
	std::vector<std::list<mex::graphModel::vindex>> ascList; // children for each variable
	mex::vector<uint32_t> maxConfig; // current best sub-configuration of MAX variables, size: nvar, initialized to be 0 instead of -1
	std::vector<varType> varTypeVec; // varTypeVec[n.id] = MAX_VAR or SUM_VAR
	mex::vector<uint32_t> tuple; // re-usable vector to accelerate, make sure reset it before use
	// indicate whether current heuristics (below that variable) are exact for a variable,
	// only applicable for SUM variables currently since we set non-zero weights for max variables
	std::vector<bool> exactHeur; // initialized to "true"
//	std::vector<unsigned> varHeightVec; // height of each variable in the pseudo tree
//	mex::VarOrder order; // v = order[i] is the ith variable eliminated
//	mex::vector<size_t> priority; // i = priority[v] is the step in which variable v is eliminated
	tree<node> exploredTree; // explored AND/OR subtree
	tree<node>::iterator root; // root = exploredTree.begin()
	unsigned int GSIZE; // size of stored tree
	double startTime;
};

#endif /* MMAP_H_ */

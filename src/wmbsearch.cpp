/*
 * wmbsearch.cpp
 *
 *  Created on: Aug 31, 2015
 *      Author: qlou
 */

#include "wmbsearch.h"

wmbsearch::wmbsearch(mex::wmbe& wmblb, mex::wmbe& wmbub, int v, std::string p, double st, double t, double r, double lb, double ub, size_t nv, double mb, double tb) :
wmbLB(wmblb), wmbUB(wmbub),
LB(lb), UB(ub), nvar(nv), verbose(v),
priorityName(p), ascList(wmbUB.nvar()), startIter(st),
timeUnit(t), timeRatio(r), memoryBudget(mb), timeBudget(tb){
	// TODO Auto-generated constructor stub
	if ( priorityName.empty() ) {
		priorityName = "absolute_gap"; // default
	}

	if (timeRatio < 1.0) {
		if (verbose > 2) {
			std::cout << "timeRatio = " << timeRatio << " < 1.0, reset to 1.0!" << std::endl;
		}
		timeRatio = 1.0;
	}
	std::cout<<" ==== Running WMB search process ====\n";
}
wmbsearch::~wmbsearch() {
	// TODO Auto-generated destructor stub
	std::cout<<" ==== Done WMB search process ====\n";
}
void wmbsearch::buildPseudoTree() {
	mex::vector<mex::graphModel::vindex> parents = wmbUB.getPseudotree();

	assert(ascList.size()==nvar); // debug
	assert((ascList[0]).empty()); // debug
	assert((ascList[nvar-1]).empty()); // debug
	for (size_t i=0; i<parents.size(); ++i) {
		mex::graphModel::vindex par = parents[i];
		if (par < 0 || par >= nvar) {
			// nodes without parents
			rts.push_back(i);
			continue;
		}
		(ascList[par]).push_back(i);
	}
	std::cout<<" ==== Pseudo tree built ====\n";

	if (verbose > 2) {
		std::cout << "no. of nodes without parents: " << rts.size() << "/"
				<< nvar << "\n";
		for (size_t i = 0; i < ascList.size(); ++i) {
			std::cout << "variable " << i << " has no. of children: " << (ascList[i]).size() << "\n";
		}
	}
}

void wmbsearch::setRootConfig() {
	// initialize root
	// this is a virtual root node, nodes without parents are children of this virtual root node
	searchnode vtRoot; // virtual root
	vtRoot.type = AND_NODE; // virtual AND_NODE node
	vtRoot.id = -1; // flag
	// for safety
//	vtRoot.ub = std::numeric_limits<double>::max();
//	vtRoot.lb = std::numeric_limits<double>::lowest();
//	vtRoot.cur_lb = vtRoot.lb = LB;
//	vtRoot.cur_ub = vtRoot.ub = UB;
	vtRoot.lb = LB;
	vtRoot.ub = UB;
	//
	exploredTree.insert(exploredTree.begin(), vtRoot);
	GSIZE = 1;

	std::cout<<" ==== Root configuration set ====\n";
}

void wmbsearch::writeNodeToFile(unsigned id, std::ostringstream& oss)  {
	oss << "(" << id;
	for (std::list<mex::graphModel::vindex>::iterator it = (ascList[id]).begin(); it != (ascList[id]).end(); ++it) {
		writeNodeToFile(*it, oss);
	}
	oss << ")";
}

void wmbsearch::writePseudoTreeToFile(std::string of_name)  {
	std::ostringstream oss;
	int root = -1; // virtual root

	oss << "(" << root;
	for (auto iter = rts.begin(); iter != rts.end(); ++iter) {
		writeNodeToFile(*iter, oss);
	}
	oss << ")";

	std::ofstream of;
	of.open(of_name.c_str(), std::ios_base::out);
	of << oss.str() << std::endl;
	of.close();
	std::cout << "Pseudo tree has been written into " << of_name <<"!\n";
}

void wmbsearch::printTree() {
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~printing the whole explored tree~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	for (auto iter = exploredTree.begin(); iter != exploredTree.end(); ++iter) {
		auto n = *iter;

		std::cout <<"variable "<< n.id << " with value " << n.val << ", solved = " << n.solved << "\n";
//		std::cout << n.ub << ", " << n.leaf_ub << ", " << n.topPath_ub << ", " << n.downPath_ub << ", " << n.exact_tub << ", " <<"\n";
//		std::cout << n.lb << ", " << n.leaf_lb << ", " << n.topPath_lb << ", " << n.downPath_lb << ", " << n.exact_tlb << ", " <<"\n";
//
		std::cout << n.ub << ", " << n.topPath_ub << ", " << n.downPath_ub << ", " << n.exact_val << ", " <<"\n";
		std::cout << n.lb << ", " << n.topPath_lb << ", " << n.downPath_lb << ", " << n.exact_val << ", " <<"\n";
	}

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}
void wmbsearch::printNode(const searchnode& n) {
	/*
	 * print information of a node
	 */
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	auto type = (n.type==AND_NODE)? "AND" : "OR";
	auto isReExpand  = n.isReExpand? "re-expanded" : "non-re-expanded";
	std::cout << isReExpand << " " << type << " node, "  << "stored priority = " << n.priority
			<< ", priority = " << calcPriority(n) << ", worstPrior = " << calcWorstPriority(n) <<  std::endl;
	std::cout <<"variable "<< n.id << " with value " << n.val << ", solved = " << n.solved << "\n";
	std::cout << n.lb << ", " << n.ub << ", " << n.exact_val << ",\n";
//	std::cout << n.cur_lb << ", " << n.cur_ub << ", " << n.exact_val << ",\n";
	std::cout << n.topPath_lb << ", " << n.topPath_ub << ",\n";
	std::cout << n.downPath_lb << ", "  << n.downPath_ub << ",\n";
	std::cout << n.downWorstPath_lb << ", " << n.downWorstPath_ub << ".\n";
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void wmbsearch::printPath(tree<searchnode>::iterator ptr) {
	while ((*ptr).toChild != NULL) {
		ptr = (*ptr).toChild;
	}

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~printing the path from the leaf to the root~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	while ( (*ptr).id > -1 ) {
		auto ptrPar = exploredTree.parent(ptr);
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptrPar); sib != exploredTree.end(ptrPar); ++sib) {
			std::cout <<"variable "<< (*sib).id << " with value " << (*sib).val << ", solved = " << (*sib).solved << "\n";
			std::cout << (*sib).ub << ",  " << (*sib).topPath_ub << ", " << (*sib).downPath_ub << "\n";
			std::cout << (*sib).lb << ",  " << (*sib).topPath_lb << ", " << (*sib).downPath_lb << "\n";
		}
		ptr = ptrPar;
	}
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

double wmbsearch::exploredRatio() {
	double logRatio = 0.0; // log ratio
	double logSpaceSize = 0.0; // in log
	int exploredSize = 0; // not in log,
	for (int id = 0; id < int(wmbUB.nvar()); ++id) {
		auto X = wmbUB.var(id);
		logSpaceSize += std::log(double(X.states()));
	}
	for (auto it = exploredTree.begin(); it != exploredTree.end(); ++it) {
		if ((*it).type == AND_NODE) {
			++ exploredSize;
		}
	}
	-- exploredSize; // exclude the root

	logRatio = std::log(exploredSize) - logSpaceSize;

	if (verbose) {
		std::cout << std::scientific;
		std::cout << "The explored graph has size: " << exploredSize
				<< ", compared to search space size: " << std::exp(logSpaceSize)
				<< ", ratio = " << std::exp(logRatio) << "!\n";
	}
	return logRatio;
}

int wmbsearch::debugMsg(const tree<searchnode>::iterator ptr) {
	double eps = 1e-10;
	auto n = *ptr;
	// greater by some margin
	bool flag = n.lb > (n.ub+eps)  || n.topPath_lb > (n.topPath_ub+eps) || n.downPath_lb > (n.downPath_ub+eps);
	flag = flag || std::isnan(n.lb) || std::isnan(n.ub);
	flag = flag || std::isnan(n.topPath_lb) || std::isnan(n.topPath_ub) || std::isnan(n.downPath_lb) || std::isnan(n.downPath_ub);
	flag = flag || std::isinf(n.lb) || std::isinf(n.ub);
	flag = flag || std::isinf(n.topPath_lb) || std::isinf(n.topPath_ub) || std::isinf(n.downPath_lb) || std::isinf(n.downPath_ub);
//	flag = flag || n.solved;
	if (flag) {
		std::cout <<"variable "<< n.id << " with value " << n.val << ", solved = " << n.solved << "\n";
		std::cout << n.ub <<", "<< n.exact_val << ", " << n.topPath_ub << ", " << n.downPath_ub
				<< ", " << (*ptr).downPath_ub + (*ptr).topPath_ub + (*ptr).exact_val << "\n";
		std::cout << n.lb <<", "<< n.exact_val << ", " << n.topPath_lb << ", " << n.downPath_lb
				<< ", " <<  (*ptr).downPath_lb + (*ptr).topPath_lb + (*ptr).exact_val << "\n";
		std::cout << "prior = " << calcPriority(ptr) << "\n";
		//
		if ( n.lb > n.ub ) std::cout << "warning: n.lb > n.ub !\n";
//		if ( n.pre_lb > n.pre_ub ) std::cout << "warning: n.pre_lb > n.pre_ub !\n";
		if (n.topPath_lb > n.topPath_ub) std::cout << "warning: n.topPath_lb > n.topPath_ub !\n";
		if (n.downPath_lb > n.downPath_ub) std::cout << "warning: n.downPath_lb > n.downPath_ub !\n";

//		std::cout << "bounds for this node may not be correct!\n";
		printPath(ptr);
		throw std::runtime_error("debug: bounds are problematic!");
		return 1;
	}
	return 0;
}

bool wmbsearch::isPruned(tree<searchnode>::iterator ptr, bool isLeaf) {
	// check whether a subtree rooted at the given node has to be deleted from the exploredTree
	double pruneThresh = 1e-10; // threshold for pruning
	bool flag = (*ptr).solved  || ( (*ptr).ub - (*ptr).lb < pruneThresh );
	flag = flag || ( ~isLeaf  && (*ptr).toChild == NULL );  // TODO qlou ~isLeaf evaluates to true always?
	return flag;
}

double wmbsearch::calcPriority(const tree<searchnode>::iterator ptr) {
	// evaluate the local priority to determine which child is most promising
	// should work for both AND_NODE and OR_NODE nodes

	// not defined for the virtual root node
	assert( (*ptr).id > -1 );


	if ( (*ptr).solved ) {
		// if solved, return minimum
		return std::numeric_limits<double>::lowest();
//		return std::log(0.0); // -inf
	}
	//
	double prior = 0.0, ub = 0.0, lb = 0.0;
	if (priorityName == "absolute_gap") {
		// absolute gap between upper and lower bounds, i.e., "upper_minus_lower"
		if ((*ptr).toChild != NULL) {
			// for AND_NODE nodes,  topPath_ub, topPath_lb, downPath_ub, downPath_lb do not contain tub, tlb.
			// for OR_NODE nodes, tub = tlb = 0.
//			ub = (*ptr).leaf_ub + (*ptr).downPath_ub + (*ptr).topPath_ub + (*ptr).tub;
//			lb = (*ptr).leaf_lb + (*ptr).downPath_lb + (*ptr).topPath_lb + (*ptr).tlb;
			// tub (= tlb) is part of exact_val
			if ( (*ptr).type == AND_NODE ) {
				lb = (*ptr).downPath_lb + (*ptr).topPath_lb + (*ptr).exact_val;
				ub = (*ptr).downPath_ub + (*ptr).topPath_ub + (*ptr).exact_val;
			} else {
				// for OR_NODE node, exact_val = -inf for initialization
				lb = (*ptr).downPath_lb + (*ptr).topPath_lb;
				ub = (*ptr).downPath_ub + (*ptr).topPath_ub;
			}
		} else {
			// if it is a leaf AND_NODE node in OPEN
			// for leaf AND_NODE nodes in OPEN, leaf_ub, leaf_lb have tub, tlb to be part of them; and downPath_ub = downPath_lb = 0.
//			ub = (*ptr).leaf_ub + (*ptr).topPath_ub;
//			lb = (*ptr).leaf_lb + (*ptr).topPath_lb;
			// (*ptr).ub = (*ptr).leaf_ub for leaf AND_NODE node,
			// for frontier OR_NODE node, existing in the follow-up DFS, the following formula are more general

			// note that exact_val is part of n.ub, n.lb
			lb = (*ptr).lb + (*ptr).topPath_lb;
			ub = (*ptr).ub + (*ptr).topPath_ub;
		}
		prior = lb - ub;
		prior = (prior < 0.0) ? prior : 0.0;
		prior = ub + std::log(1 - std::exp(prior));
	}
	else if (priorityName == "absolute_upper") {
		// absolute upper bounds
		if ((*ptr).toChild != NULL) {
//			ub = (*ptr).leaf_ub + (*ptr).downPath_ub + (*ptr).topPath_ub + (*ptr).tub;
			if ( (*ptr).type == AND_NODE ) {
				ub = (*ptr).downPath_ub + (*ptr).topPath_ub + (*ptr).exact_val;
			} else {
				ub = (*ptr).downPath_ub + (*ptr).topPath_ub;
			}

		}
		else {
//			ub = (*ptr).leaf_ub + (*ptr).topPath_ub;
			ub = (*ptr).ub + (*ptr).topPath_ub;
		}
		prior = ub;
	}
	// do not support those priority types any  more
//	else if (priorityName == "relative_gap") {
//		ub = (*ptr).leaf_ub;
//		lb = (*ptr).leaf_lb;
//
//		if ( (*ptr).toChild == NULL ) {
//			// frontier nodes, possibly OR_NODE nodes in the follow-up DFS
//			ub = (*ptr).ub;
//			lb = (*ptr).lb;
//		}
//
//		prior = lb - ub;
//		prior = (prior < 0.0) ? prior : 0.0;
//		prior = ub + std::log(1 - std::exp(prior));
//	}
//	else if (priorityName == "relative_upper") {
//		// absolute gap between upper bounds
//
//		ub = (*ptr).leaf_ub;
//		if ( (*ptr).toChild == NULL ) {
//			// frontier nodes, possibly OR_NODE nodes in the follow-up DFS
//			ub = (*ptr).ub;
//		}
//
//		prior = ub;
//	}
	else {
		throw std::runtime_error("no such priorityName defined!");
	}
	if (std::isnan(prior)) {
		std::cout << "prior = " << prior << "\n";
		std::cout <<"variable "<< (*ptr).id << " with value " << (*ptr).val << ", solved = " << (*ptr).solved << "\n";
		std::cout << (*ptr).ub << ", " << (*ptr).topPath_ub << ", " << (*ptr).downPath_ub << "\n";
		std::cout << (*ptr).lb << ", " << (*ptr).topPath_lb << ", " << (*ptr).downPath_lb << "\n";
		throw std::runtime_error("priority is not proper!");
	}
	return prior;
}

tree<wmbsearch::searchnode>::iterator wmbsearch::forwardPass() {
	// forward pass to find the node in OPEN with the highest priority
	mex::vector<uint32_t> tuple(nvar);
	auto ptr = exploredTree.begin();
	// thetas on the path
//	double costUB = 0.0, costLB = 0.0;
//	// bounds for those OR_NODE siblings of AND_NODE nodes along the path, i.e., branches
//	double branchUB = 0.0, branchLB = 0.0;
	// bounds for current OR_NODE nodes on the path
	// forward pass, also update topPath_ub, topPath_lb
	if (verbose>2) {
//		debugMsg(ptr);
	}
	//

	while ( (*ptr).toChild != NULL ) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		auto ptrPar = ptr;
		ptr = (*ptr).toChild;

		if (  (*ptr).type == AND_NODE ) {
			// ptr is AND_NODE
			// ptrPar is OR_NODE
			// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
			(*ptr).topPath_ub = (*ptrPar).topPath_ub;
			(*ptr).topPath_lb = (*ptrPar).topPath_lb;

			tuple[(*ptr).X] = (*ptr).val;
//			// current tub, tlb are not part of topPath_ub or toPath_lb
		} else {
			// ptr is OR_NODE
			// ptrPar is AND_NODE
			// tub, tlb already included in ub, lb
			(*ptr).topPath_ub = (*ptrPar).topPath_ub + (*ptrPar).ub - (*ptr).ub;
			(*ptr).topPath_lb = (*ptrPar).topPath_lb + (*ptrPar).lb - (*ptr).lb;


		}
		// debug
//		if (verbose>2) debugMsg(ptr);
	}
	// debug
	assert( ((*ptr).type == AND_NODE)   &&  ((*ptr).toChild == NULL));
	//
	std::list<mex::graphModel::vindex> childrenList;
	if ( (*ptr).id  > -1 ) {
		childrenList = ascList[(*ptr).id];
	}
	else {
		// if it is the virtual root
		childrenList = rts;
	}
	// if the most promising node is a node without children,
	// the whole search progress should be terminated, although it's possible that the bound gap remains non-zero.
	// it depends on what priority you use
	if (childrenList.empty()) {
		std::cout <<"variable "<< (*ptr).id << " with value " << (*ptr).val << " has no children!\n";
		std::cout << (*ptr).ub << ", "  << (*ptr).topPath_ub << ", " << (*ptr).downPath_ub << "\n";
		std::cout << (*ptr).lb << ", "  << (*ptr).topPath_lb << ", " << (*ptr).downPath_lb << "\n";
		std::cout << "priority = " << calcPriority(ptr) << "\n";

		while ( (*ptr).id > -1 ) {
			std::cout <<"variable "<< (*ptr).id << " with value " << (*ptr).val << ", solved = " << (*ptr).solved << "\n";
			ptr = exploredTree.parent(ptr);
		}
		std::cout <<"variable "<< (*ptr).id << " with value " << (*ptr).val << ", solved = " << (*ptr).solved << "\n";
		return NULL;
	}
	// update bounds later
	(*ptr).ub = (*ptr).lb = (*ptr).exact_val;

	//
	bool isProp = false;
	//
	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
		searchnode child;
		child.type = OR_NODE;
		child.id = *it;
		child.X = wmbUB.var(child.id);
		// for OR_NODE node
		child.exact_val = -std::numeric_limits<double>::infinity();
		// append to the explored tree
		auto ptrChild = exploredTree.append_child(ptr, child);
		++ GSIZE;
		//
		mex::Factor dUB(child.X,0.0);
		mex::Factor dLB(child.X,0.0);
		//
		mex::vector<uint32_t> childTuple(tuple);
		//
		int nMini = wmbUB.duplicateInBucket(child.X, childTuple);
		isProp = isProp || (nMini > 1);

		if (verbose > 2) {

			if ( nMini >1 ) {
				std::cout<< "Bounds should be improved when expanding this node, no. of mini-buckets: " << nMini << "\n";
			}
		}


		for (size_t v = 0; v < (child.X).states(); ++v) {
			searchnode kid;
			kid.type = AND_NODE;
			kid.X = child.X;
			kid.id = child.id;
			kid.val = v;

			childTuple[child.X] = v;
			// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
			// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
			kid.exact_val = wmbUB.heuristicTheta(kid.X, childTuple);
			// actually, we can speed up this evaluation by computing each factor only once
			// for a leaf node, leaf_ub, leaf_lb are its ub, lb
//			kid.leaf_ub = kid.ub = dUB[v] = wmbUB.heuristicIn(kid.X, childTuple) + kid.tub;
//			kid.leaf_lb = kid.lb = dLB[v] = wmbLB.heuristicIn(kid.X, childTuple) + kid.tlb;
			kid.downPath_lb = wmbLB.heuristicIn(kid.X, childTuple);
			kid.downPath_ub = wmbUB.heuristicIn(kid.X, childTuple);

			kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
			kid.ub = dUB[v] = kid.downPath_ub + kid.exact_val;

			// whether solved
			kid.solved = checkSolved(kid);
			// if solved, no need to append to the exploredTree
			if (kid.solved) {
				double mx = std::max( (*ptrChild).exact_val, kid.ub );
				(*ptrChild).exact_val = mx + log( exp( (*ptrChild).exact_val - mx ) + exp( kid.ub -mx ) );
				continue;
			}
			// append
			exploredTree.append_child(ptrChild, kid);
			++ GSIZE;
		}
		//
		(*ptrChild).ub = dUB.logsumexp();
		(*ptrChild).lb = dLB.logsumexp();
		//
		(*ptr).ub += (*ptrChild).ub;
		(*ptr).lb += (*ptrChild).lb;
		// if no child left, safely delete this OR_NODE node
		if ( exploredTree.begin(ptrChild) == exploredTree.end(ptrChild) ) {
			(*ptr).exact_val += (*ptrChild).ub;
			exploredTree.erase(ptrChild);
			-- GSIZE;
		}
	}
	// propagate bounds if necessary
	if (isProp) propagateBounds(ptr);

	// this node is solved only when all its children has been solved
	// note: none of the solved children has been appended
	(*ptr).solved = ( exploredTree.begin(ptr) == exploredTree.end(ptr) );
	if ( (*ptr).solved ) return ptr;

	// set best child here, local forward, backward pass
	double bestPrior =  -std::numeric_limits<double>::infinity();
	tree<searchnode>::iterator toChild = NULL;
	for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		tree<searchnode>::iterator toKid = NULL;
		double prior = -std::numeric_limits<double>::infinity();
		// tub is part of ub
		(*sib).topPath_ub = (*ptr).topPath_ub + (*ptr).ub - (*sib).ub;
		(*sib).topPath_lb = (*ptr).topPath_lb + (*ptr).lb - (*sib).lb;

		for (tree<searchnode>::sibling_iterator nb = exploredTree.begin(sib); nb != exploredTree.end(sib); ++nb) {
			// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
			(*nb).topPath_ub = (*sib).topPath_ub;
			(*nb).topPath_lb = (*sib).topPath_lb;
			double p = calcPriority(nb);
			if ( (toKid == NULL) || ( (toKid != NULL) && (p > prior) ) ) {
				toKid = nb;
				prior = p;
			}
		}
		//
//		(*sib).solved = nb_solved;
		//
		(*sib).toChild = toKid;
//		(*sib).leaf_ub = (*toKid).leaf_ub;
//		(*sib).leaf_lb = (*toKid).leaf_lb;

		// first level OR node
		sib->downPath_lb = toKid->lb;
		sib->downPath_ub = toKid->ub;
		if ( (toChild == NULL) || ( (toChild != NULL) && (prior > bestPrior) ) ) {
			toChild = sib;
			bestPrior = prior;
		}
	}
//	(*ptr).solved = sib_solved;

	// if not solved
	(*ptr).toChild = toChild;
//	(*ptr).leaf_ub = (*toChild).leaf_ub;
//	(*ptr).leaf_lb = (*toChild).leaf_lb;

	// not including exact_val
//	(*ptr).downPath_ub = (*ptr).ub - (*ptr).exact_val - (*toChild).ub;
//	(*ptr).downPath_lb = (*ptr).lb - (*ptr).exact_val - (*toChild).lb;
	ptr->downPath_lb = toChild->downPath_lb + ptr->lb - ptr->exact_val - toChild->lb;
	ptr->downPath_ub = toChild->downPath_ub + ptr->ub - ptr->exact_val - toChild->ub;
	// return
	return ptr;
}
int wmbsearch::bottomUpPass(tree<searchnode>::iterator ptr) {
	// this function is an updated version of backwardPass, run after DFS, will substitute backwardPass
	if (ptr == NULL) {
		std::cout << "NULL occurs in the bottomUpDFS, problematic or in purpose!\n";
		return -1;
	}
	// double check

//	std::cout << "(*ptr).type = " << (*ptr).type << "\n";
//	assert((*ptr).type == AND_NODE);
	// the returned ptr is not solved
	// so if run this function after DFS, we will get input ptr's ancestor
	// doPruning will possibly make the returned node points to an arbitrary child.
	// no worries, we will re-calculate the best child here

	// caveat: make sure bounds have been fully updated before pruning
//	propagateBounds(ptr);
	ptr = doPruning(ptr, exploredTree.begin());
	// check whether the whole tree has been solved or not
	auto rt = exploredTree.begin();
	if ( rt->solved ) {
		std::cout << "the root has been solved" << std::endl;
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		return 1;
	}

	// after pruning, it's possible to be either AND_NODE or OR_NODE
	//
	// no need to do pruning in the while loop
	// also, no need to update "solved" in the while loop because if any ancestor is solved,
	// it should be already deleted in the pruning step.
	while (true) {
		// actually, we do not necessarily update bounds here
		// propagateBounds should have done the job
		// debug
		//if (verbose>2) debugMsg(ptr);
		//
		// initialize with the smallest value
		double bestPrior = -std::numeric_limits<double>::infinity();
		tree<searchnode>::sibling_iterator toChild = NULL;
		// possibly to be assigned later
		(*ptr).toChild = NULL;
		// traverse all the children to determine which child to point to.
		// priority can be recovered at any ancestor
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//			solved = solved && (*sib).solved;
			if ((*ptr).type == AND_NODE) {
				// AND_NODE node, exact_val is part of ub, lb
				(*sib).topPath_ub = (*ptr).topPath_ub + (*ptr).ub - (*sib).ub;
				(*sib).topPath_lb = (*ptr).topPath_lb + (*ptr).lb - (*sib).lb;
				double prior = calcPriority(sib);

				if ( (toChild == NULL) || ( (toChild != NULL) && (prior > bestPrior) )  ) {
					toChild = sib;
					bestPrior = prior;
					(*ptr).toChild = sib;
//					(*ptr).leaf_ub = (*sib).leaf_ub;
//					(*ptr).leaf_lb = (*sib).leaf_lb;
					// including bounds from other OR_NODE children
					// downward not including "tub", "tlb" for AND_NODE nodes
//					(*ptr).downPath_ub = (*sib).downPath_ub + (*ptr).ub - (*ptr).tub - (*sib).ub;
//					(*ptr).downPath_lb = (*sib).downPath_lb + (*ptr).lb - (*ptr).tlb - (*sib).lb;

					// tub (= tlb) is part of exact_val
					(*ptr).downPath_ub = (*sib).downPath_ub + (*ptr).ub - (*ptr).exact_val - (*sib).ub;
					(*ptr).downPath_lb = (*sib).downPath_lb + (*ptr).lb - (*ptr).exact_val - (*sib).lb;
				}
			}
			else {
				// ptr is OR node
				// use most update-to-date information, topPath for AND_NODE siblings are the same
				// equal to those of the OR_NODE parent
				(*sib).topPath_ub = (*ptr).topPath_ub;
				(*sib).topPath_lb = (*ptr).topPath_lb;
				double prior = calcPriority(sib);

				if ( (toChild == NULL) || ( (toChild != NULL) && (prior > bestPrior) ) ) {
					toChild = sib;
					bestPrior = prior;
					(*ptr).toChild = sib;
//					(*ptr).leaf_ub = (*sib).leaf_ub;
//					(*ptr).leaf_lb = (*sib).leaf_lb;
					// including tub, tlb from the AND_NODE children
					// caveat: the AND_NODE child should NOT be a node in OPEN
					if ( (*sib).toChild != NULL ) {
//						(*ptr).downPath_ub = (*sib).tub + (*sib).downPath_ub;
//						(*ptr).downPath_lb = (*sib).tlb + (*sib).downPath_lb;
						// tub (=tlb) is part of exact_val
						(*ptr).downPath_ub = (*sib).exact_val + (*sib).downPath_ub;
						(*ptr).downPath_lb = (*sib).exact_val + (*sib).downPath_lb;
					}
					else {
						// if sib is a frontier AND node
//						(*ptr).downPath_ub = (*sib).downPath_ub; // = 0.0
//						(*ptr).downPath_lb = (*sib).downPath_lb; // = 0.0

						ptr->downPath_lb = sib->lb;
						ptr->downPath_ub = sib->ub;
					}
				}
			}
		}

		// debug
//		if ( verbose > 2 ) {
//			if ( (*ptr).id > -1 )
//				std::cout <<"node type: "<<  (*ptr).type << ", id (val): "
//				<< (*ptr).id << " ( "<< (*ptr).val << " )"
//				<< ", priority: " <<  calcPriority(ptr) << "\n";
////			if ((*ptrPar).toChild != ptr) {
////				std::cout << "path changed!\n";
////			}
//		}
		// go upward
		if ( ptr == exploredTree.begin() ) {
			break;
		}
		ptr = exploredTree.parent(ptr);
	}

//	(*ptr).solved = true;
//	// be careful, make sure the root has been expanded before
//	for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//		if ( ! (*sib).solved ) {
//			(*ptr).solved = false;
//			break;
//		}
//	}

	return 0;
}

tree<wmbsearch::searchnode>::iterator wmbsearch::DFS(const double timeLimit) {
	/* do DFS when the memory budget is (almost) used up
	 *	DFS implemented by stack, O(branch_factor * max_depth) memory required
	 *	DFS may be implemented s.t. the memory usage is only linear to max_depth
	 *	the order of the children is sorted (ascending order) by priority when pushed into the stack
	 */
	mex::vector<uint32_t> tuple(nvar); // a little tricky to use in DFS
	auto ptr = exploredTree.begin();
	// thetas on the path

	while ( (*ptr).toChild != NULL ) {
			// note that the virtual root is also AND_NODE with id = -1
			// first node in the loop would be OR_NODE node
			auto ptrPar = ptr;
			ptr = (*ptr).toChild;

			if (  (*ptr).type == AND_NODE ) {
				// ptr is AND_NODE
				// ptrPar is OR_NODE
				// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
				(*ptr).topPath_ub = (*ptrPar).topPath_ub;
				(*ptr).topPath_lb = (*ptrPar).topPath_lb;

				tuple[(*ptr).X] = (*ptr).val;
	//			// current tub, tlb are not part of topPath_ub or toPath_lb
			} else {
				// ptr is OR_NODE
				// ptrPar is AND_NODE
				(*ptr).topPath_ub = (*ptrPar).topPath_ub + (*ptrPar).ub - (*ptr).ub;
				(*ptr).topPath_lb = (*ptrPar).topPath_lb + (*ptrPar).lb - (*ptr).lb;
			}
			// debug
//			if (verbose>2) debugMsg(ptr);
	}
	// safety check
	assert( ((*ptr).type == AND_NODE)   &&  ((*ptr).toChild == NULL) );
	(*ptr).solved = checkSolved( *ptr );
	if ( (*ptr).solved ) {
		std::cout << "The node is already solved before doing DFS" << std::endl;
		return NULL;
	}

	// push current best leaf node to the stack
	// this node won't be deleted during pruning inside while loop
	std::stack<tree<searchnode>::iterator> stack;
	assert(stack.empty());
	stack.push(ptr);
	std::list<mex::graphModel::vindex> childrenList;
	auto ancestor = ptr;

	// we push both AND_NODE, OR_NODE nodes to the stack
	while (!stack.empty()) {
		ptr = stack.top();
		stack.pop();

		if (verbose > 2) {
//			std::cout<< "node type: " << (*ptr).type << "\n";
			std::cout.precision(10);
			std::cout << "DFS before bound propagation: " << std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree Size: " << GSIZE << std::endl;
		}

		if ((*ptr).type == AND_NODE) {
			// AND_NODE nodes, expand its OR_NODE children, update bounds
			if ((*ptr).id > -1) {
				childrenList = ascList[(*ptr).id];
			} else {
				// if it is the virtual root
				childrenList = rts;
			}

			if ( (*ptr).solved ) {
				doPruning(ptr, ancestor);
				continue;
			}
			// add to tuple, be careful
			if ((*ptr).id > -1) tuple[ (*ptr).X ] = (*ptr).val;

			// update bounds later

			(*ptr).ub = (*ptr).exact_val;
			(*ptr).lb = (*ptr).exact_val;

			// whether propagate
			bool isProp = false;

			for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
				searchnode child;
				child.type = OR_NODE;
				child.id = *it;
				child.X = wmbUB.var(child.id);
				child.exact_val = -std::numeric_limits<double>::infinity();
				//
				mex::Factor dUB(child.X, 0.0);
				mex::Factor dLB(child.X, 0.0);
				//
				mex::vector<uint32_t> childTuple(tuple);
				isProp = isProp || (wmbUB.duplicateInBucket(child.X, childTuple) > 1);

				for (size_t v = 0; v < (child.X).states(); ++v) {
					searchnode kid;
					kid.type = AND_NODE;
					kid.X = child.X;
					kid.id = child.id;
					kid.val = v;

					childTuple[child.X] = v;
					// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
					// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
//					kid.tub = wmbUB.heuristicTheta(kid.X, childTuple);
//					kid.tlb = wmbLB.heuristicTheta(kid.X, childTuple);
//					//
//					kid.exact_val = kid.tub;

					kid.exact_val = wmbUB.heuristicTheta(kid.X, childTuple);
					// actually, we can speed up this evaluation by computing each factor only once
					// for a leaf node, leaf_ub, leaf_lb are its ub, lb
//					kid.leaf_ub = kid.ub = dUB[v] =
//							wmbUB.heuristicIn(kid.X, childTuple) + kid.tub;
//					kid.leaf_lb = kid.lb = dLB[v] =
//							wmbLB.heuristicIn(kid.X, childTuple) + kid.tlb;
					kid.downPath_lb = wmbLB.heuristicIn(kid.X, childTuple);
					kid.downPath_ub = wmbUB.heuristicIn(kid.X, childTuple);

					kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
					kid.ub = dUB[v] = kid.downPath_ub + kid.exact_val;
				}

				// OR_NODE node as a leaf node
				child.downPath_ub = child.ub = dUB.logsumexp();
				child.downPath_lb = child.lb = dLB.logsumexp();

				(*ptr).ub += child.ub;
				(*ptr).lb += child.lb;

				child.solved = checkSolved(child);
				if (child.solved) {
					(*ptr).exact_val += child.ub; // child.ub = child.lb
					continue;
				}
				// append to the explored tree if not solved yet
				// no solved nodes will be pushed into the stack
				exploredTree.append_child(ptr, child);
				++ GSIZE;
			}
			// propagate bounds to the root before pruning,
			// only need for AND_NODE
			// propagate bounds only when necessary
			if (isProp) {
				propagateBounds(ptr);
				if (verbose > 2) {
					std::cout.precision(10);
					std::cout << "DFS after bound propagation: " << std::setw(10) << LB << " < ln Z < " << UB << "\n";
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
				}
			}
			// if no children appended, node is solved, do pruning
			if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
				(*ptr).solved = true;
				doPruning(ptr, ancestor);
				continue;
			}

			// calculate the priority
			// sort the OR_NODE children based on priority
			// all those nodes are not yet solved
			std::vector<double> childrenPriorities;
			std::vector<tree<searchnode>::iterator> childrenPtrs;
			assert( childrenPriorities.empty() && childrenPtrs.empty() );
			for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				//
				(*sib).topPath_ub = (*ptr).topPath_ub + (*ptr).ub - (*sib).ub;
				(*sib).topPath_lb = (*ptr).topPath_lb + (*ptr).lb - (*sib).lb;

				childrenPriorities.push_back( calcPriority(sib) );
				childrenPtrs.push_back(sib);
			}
			// ascending order
			auto sortedIds = sortIndexes(childrenPriorities);
			// not quite necessary in DFS
			(*ptr).toChild = childrenPtrs[ sortedIds.back() ];

			// thus high-priority nodes will be popped out first
			for (auto id : sortedIds) {
				stack.push(childrenPtrs[id]);
			}
		}
		else {
			// OR_NODE node, re-expand its AND_NODE children
			// no need to update bounds because we've already done when generating the OR_NODE node
			std::vector<double> childrenPriorities;
			std::vector<tree<searchnode>::iterator> childrenPtrs;
			assert( childrenPriorities.empty() && childrenPtrs.empty() );
			for (size_t v = 0; v < ((*ptr).X).states(); ++v) {
				searchnode kid;
				kid.type = AND_NODE;
				kid.X = (*ptr).X;
				kid.id = (*ptr).id;
				kid.val = v;

				tuple[(*ptr).X] = v;
				// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
				// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
//				kid.tub = wmbUB.heuristicTheta(kid.X, tuple);
//				kid.tlb = wmbLB.heuristicTheta(kid.X, tuple);
//				// actually, kid.tub = kid.tlb
//				kid.exact_val = kid.tub;
				kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
				// actually, we can speed up this evaluation by computing each factor only once
				// for a leaf node, leaf_ub, leaf_lb are its ub, lb
//				kid.leaf_ub = kid.ub = wmbUB.heuristicIn(kid.X,
//						tuple) + kid.tub;
//				kid.leaf_lb = kid.lb = wmbLB.heuristicIn(kid.X,
//						tuple) + kid.tlb;
				kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
				kid.downPath_ub = wmbUB.heuristicIn(kid.X, tuple);

				kid.lb = kid.downPath_lb + kid.exact_val;
				kid.ub = kid.downPath_ub + kid.exact_val;

				kid.topPath_ub = (*ptr).topPath_ub;
				kid.topPath_lb = (*ptr).topPath_lb;
				//
				kid.solved = checkSolved(kid);
				if (kid.solved) {
					// no need to update bounds here
					double mx = std::max( (*ptr).exact_val, kid.ub );
					(*ptr).exact_val = mx + log( exp( (*ptr).exact_val - mx ) + exp( kid.ub -mx ) );
					continue;
				}

				auto toKid = exploredTree.append_child(ptr, kid);
				++ GSIZE;

				childrenPriorities.push_back( calcPriority(toKid) );
				childrenPtrs.push_back( toKid );
			}
			// no active children
			if (childrenPriorities.empty()) {
				(*ptr).solved = true;
				doPruning(ptr, ancestor);
				continue;
			}

			auto sortedIds = sortIndexes(childrenPriorities);
			for (auto id : sortedIds) {
				stack.push(childrenPtrs[id]);
			}
		}

		if ( mex::timeSystem() - startIter > timeThresh) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//			timeThresh += timeUnit;
			timeThresh += timeUnit*timeRatio;
			timeRatio *= timeRatio;
		}


		if (mex::timeSystem() - startIter > timeLimit) {
			return NULL;
		}
	}
	// none of anscetor's descendants should exist
	assert(exploredTree.begin(ancestor) == exploredTree.end(ancestor));
	(*ancestor).solved = true;

	return ancestor;
}

std::vector<size_t> wmbsearch::sortIndexes(const std::vector<double> &v) {
	// return index according to values's ascending order
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

bool wmbsearch::checkSolved(searchnode& n) {
	// works for both AND_NODE and OR_NODE nodes
	// caveat: may only apply to frontier nodes!
	bool solved = false;
	if (n.id > -1) {
		// if not root
		solved = (ascList[n.id]).empty() || ( n.ub - n.lb <= solvedThresh );
	}
	else {
		solved = (n.ub - n.lb <= solvedThresh) ;
	}
	return solved;
}

tree<wmbsearch::searchnode>::iterator wmbsearch::doPruning(tree<searchnode>::iterator ptr, tree<searchnode>::iterator ancestor) {
/* do bottom-up pruning if necessary,
 * make sure do it from the frontiers to avoid memory issue
 * always be very careful about deleting nodes, make sure those node pointers in the stack will not be wide pointers
 * ancestor is the upper limit node. deletion stops once we reach this node even if it is solved, default: the root
 * ancestor must be an ancestor on the path
 * return the highest not being solved
 * currently, it should only apply to DFS related procedure
*/
	if ( !(*ptr).solved || (ptr == ancestor) ) return ptr;
	tree<searchnode>::iterator highestDelete = ptr;
	//
	ptr = exploredTree.parent(ptr);

	while ( ptr != ancestor ) {
		bool solved = true;
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// a node is solved only when all its children are solved
			if( !(*sib).solved ) {
				solved = false;
				break;
			}
		}
		(*ptr).solved = solved;
		// a node is solved only when all its children are solved
		if (!(*ptr).solved) {
			break;
		}
		//
		highestDelete = ptr;
		ptr = exploredTree.parent(ptr);
	}
	// erase the whole subtree rooted at this node (including itself)
	// store exact value in its parent's exact_val
	// if not the root
	ptr = exploredTree.parent(highestDelete);
	if ((*ptr).type == AND_NODE) {
		(*ptr).exact_val += (*highestDelete).ub; // at this moment, ub=lb
	} else {
		double maxval = std::max((*ptr).exact_val, (*highestDelete).ub);
		(*ptr).exact_val = maxval
				+ log(
						exp((*ptr).exact_val - maxval)
								+ exp((*highestDelete).ub - maxval));
	}
	// the parent points to some child for safety
	if ((*ptr).toChild == highestDelete) {
		(*ptr).toChild = NULL;
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			if (sib != highestDelete) {
				// no need to pick the best child at this moment, we'll do it in backward pass
				(*ptr).toChild = sib;
				break;
			}
		}
	}

	// finally delete
	int GSIZE_pre = GSIZE;
	// size gives you the size of the subtree rooted at given node ( including that node )
	GSIZE -= exploredTree.size(highestDelete);
	exploredTree.erase(highestDelete);

	if (verbose > 2) {
		std::cout << "pruning done: " << GSIZE_pre - GSIZE << " nodes deleted!\n" ;
	}

	// return parent of the highestDelete
	return ptr;
}

void wmbsearch::propagateBounds(tree<searchnode>::iterator ptr) {
	// re-calculate bounds from existing children and exact_val
	assert( ptr != NULL );
//	if (verbose >2 ) std::cout << "propagating bounds....\n";

	while( ptr != exploredTree.begin() ) {
		ptr = exploredTree.parent(ptr);
		// update bounds for the parent
//		(*ptr).pre_ub = (*ptr).ub;
//		(*ptr).pre_lb = (*ptr).lb;
		// update bounds by re-calculation.
		// more time-consuming than simple update, but more numerically stable
		if ((*ptr).type == AND_NODE) {
			(*ptr).ub = (*ptr).exact_val;
			(*ptr).lb = (*ptr).exact_val;
			// possibly some children have been pruned
			for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				(*ptr).ub += (*sib).ub;
				(*ptr).lb += (*sib).lb;
			}
		}
		else {
			// OR_NODE node
			std::vector<double> ub_vec {(*ptr).exact_val}, lb_vec {(*ptr).exact_val};
			for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				ub_vec.push_back( (*sib).ub );
				lb_vec.push_back( (*sib).lb );
			}
			(*ptr).ub = logsumexp(ub_vec);
			(*ptr).lb = logsumexp(lb_vec);
		}

		// debug
//		std::cout << "type: " << (*ptr).type << std::endl;
//		std::cout.precision(10);
//		std::cout << "lb: " << std::setw(10) << (*ptr).pre_lb << " < " << (*ptr).lb << "\n";
//		std::cout << "ub: " << std::setw(10) << (*ptr).pre_ub << " > " << (*ptr).ub << "\n";
	}

	auto rt = exploredTree.begin();
	LB = (*rt).lb;
	UB = (*rt).ub;

	if (UB-LB < boundGap) {
		std::cout.precision(10);
		std::cout << "[" << mex::timeSystem() - startIter << "]: "
				<< std::setw(10) << LB << " < ln Z < " << UB << ", UB-LB = " << UB-LB << " < " << boundGap << "\n";
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		boundGap = boundGap/10.0;
	}
}

double wmbsearch::logsumexp(std::vector<double>& vec) {
	assert(!vec.empty()); // cannot be empty
	double mx = *(std::max_element(vec.begin(), vec.end()));
	double r = 0;
	if (mx <= std::numeric_limits<double>::lowest()) {
		return mx;
	}
	for (double& val : vec) {
		r += std::exp(val - mx);
	}
	return (std::log(r) + mx);
}

/*void wmbsearch::setDownPath(tree<searchnode>::iterator ptr, tree<searchnode>::iterator sib) {
	// update information when ptr change toChild to sib
	assert( ptr == exploredTree.parent(sib) );
	(*ptr).toChild  = sib;

	if ((*ptr).type == AND_NODE) {
		// AND_NODE node, exact_val is part of ub, lb
		(*sib).topPath_ub = (*ptr).topPath_ub + (*ptr).ub - (*sib).ub;
		(*sib).topPath_lb = (*ptr).topPath_lb + (*ptr).lb - (*sib).lb;

		(*ptr).leaf_ub = (*sib).leaf_ub;
		(*ptr).leaf_lb = (*sib).leaf_lb;
		// including bounds from other OR_NODE children
		// downward not including "tub", "tlb" for AND_NODE nodes
//					(*ptr).downPath_ub = (*sib).downPath_ub + (*ptr).ub - (*ptr).tub - (*sib).ub;
//					(*ptr).downPath_lb = (*sib).downPath_lb + (*ptr).lb - (*ptr).tlb - (*sib).lb;

		// tub (= tlb) is part of exact_val
		(*ptr).downPath_ub = (*sib).downPath_ub + (*ptr).ub - (*ptr).exact_val
				- (*sib).ub;
		(*ptr).downPath_lb = (*sib).downPath_lb + (*ptr).lb - (*ptr).exact_val
				- (*sib).lb;
	} else {
		// OR_NODE node
		// use most update-to-date information, topPath for AND_NODE siblings are the same
		// equal to those of the OR_NODE parent
		(*sib).topPath_ub = (*ptr).topPath_ub;
		(*sib).topPath_lb = (*ptr).topPath_lb;

		(*ptr).leaf_ub = (*sib).leaf_ub;
		(*ptr).leaf_lb = (*sib).leaf_lb;
		// including tub, tlb from the AND_NODE children
		// caveat: the AND_NODE child should NOT be a node in OPEN
		if ((*sib).toChild != NULL) {
//						(*ptr).downPath_ub = (*sib).tub + (*sib).downPath_ub;
//						(*ptr).downPath_lb = (*sib).tlb + (*sib).downPath_lb;
			// tub (=tlb) is part of exact_val

			(*ptr).downPath_ub = (*sib).exact_val + (*sib).downPath_ub;
			(*ptr).downPath_lb = (*sib).exact_val + (*sib).downPath_lb;
		} else {
			(*ptr).downPath_ub = (*sib).downPath_ub; // = 0.0
			(*ptr).downPath_lb = (*sib).downPath_lb; // = 0.0
		}

	}
}*/

void wmbsearch::start(double timeLimit, int nSearch, std::string pstName, double addEpsilon, std::string probName, const double memLimit) {
	// build pseudo tree from the induced graph
	buildPseudoTree();
	// visualize the pseudo tree if required
	if (!pstName.empty()) {
		writePseudoTreeToFile(pstName);
	}
	// set configuration for the root node
	setRootConfig();
//	std::cout << " ==== Beginning WMB search process ====\n";
	std::cout << "Priority type is " << priorityName << "\n";
	//
//	int nodes = 0;
//	int nLines = 1e3; // no. of lines to output
//	int dividant = std::max(nSearch/nLines, 1);
//	if (verbose>2) dividant = 1;

	int indicator = 0; // whether we have to terminate the program
	tree<searchnode>::iterator ptr = NULL;
	auto ptrRoot = exploredTree.begin();
	// in byte
	int memUnit = sizeof(*ptrRoot);
	// MemLimit is in megabyte
	auto treeSizeLimit = long(1024*1024*memLimit/memUnit);

	std::cout << "Time limit (sec): " << timeLimit << "\n";
	std::cout << "Memory limit (MB): " << memLimit << "\n";
	std::cout << "Tree size limit: " << treeSizeLimit << "\n";
	// for geometric index
//	int geoIndex = 0;

	// once we start to do DFS, we will never go back
	// actually we can do in a mixed way. if after pruning the tree size is smaller than the threshold,
	// we may go back to do best-first search again.
	bool isDFS = false;

	// debug
//	nSearch = 1e3;
//	while (nodes < nSearch) {
	while(true) {
		if ( (*ptrRoot).solved ) {
			std::cout << "All nodes been instantiated!\n";
			break;
		}



		if ( GSIZE <= treeSizeLimit && !isDFS) {
			// debug

			if (verbose > 2) {
				std::cout << "~~~~~~~~~~~~~~~~~now forward pass~~~~~~~~~~~~~~~~~\n";
				std::cout << "Tree size: " << GSIZE << std::endl;
			}
			// forward pass
			ptr = forwardPass();
//		printTree();
			if (verbose > 2) {
				std::cout << "~~~~~~~~~~~~~~~~~now backward pass~~~~~~~~~~~~~~~~~\n";
				std::cout << "Tree size: " << GSIZE << std::endl;
			}
			// backward pass
//			indicator = backwardPass(wmbUB, wmbLB, ptr);
			// bottomUpPass is a more up-to-date alternative
			indicator = bottomUpPass(ptr);
//		printTree();
			if (indicator < 0) {
				std::cout<< "The most promising node in OPEN has no children, loop terminated!\n";
				break;
			}

		} else {
			// reach tree size limit (memory limit)
			if (!isDFS) {
				// only prompt at the first time switch to DFS
				std::cout << "reach the memory limit (MB): "<< memLimit << ", switch to DFS!\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE <<"\n";
				std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
			}
			isDFS = true;

			if (verbose > 2) {
					std::cout << "~~~~~~~~~~~~~~~~~now DFS~~~~~~~~~~~~~~~~~\n";
			}
			ptr = DFS(timeLimit);

			if (verbose > 2) {
					std::cout << "~~~~~~~~~~~~~~~~~now bottomUpPass~~~~~~~~~~~~~~~~~\n";
			}
			indicator = bottomUpPass(ptr);

			if ( indicator == 1 ) {
				std::cout << "problem solved in DFS!\n";
				break;
			}
//
//			if ( indicator == -1  ) {
//				break;
//			}

		}

		UB = (*ptrRoot).ub;
		LB = (*ptrRoot).lb;


//
//		if ( ((nodes % dividant == 0) && (verbose > 0) && (verbose !=2 )) || ( (nodes == std::pow(2,geoIndex) ) && (verbose == 2) ) ) {
//
//			++ geoIndex;
//			std::cout<< "nodes explored: " << nodes << "\n";
//			std::cout.precision(10);
//			std::cout << "[" << mex::timeSystem() - startIter << "]: "
//			<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
//			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//		}
//		++ nodes;
//		if ( mex::timeSystem() - startIter > )

		if ( mex::timeSystem() - startIter > timeThresh) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//			timeThresh += timeUnit;
			timeThresh += timeUnit*timeRatio;
			timeRatio *= timeRatio;
		}


		if (UB - LB < boundGap) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << ", UB-LB = "
					<< UB - LB << " < " << boundGap << "\n";
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
			boundGap = boundGap / 10.0;
		}


		if ( UB - LB < solvedThresh ) {
			std::cout << "the gap between UB and LB is less than "<< solvedThresh << ", quit!\n";
			break;
		}

		if ( mex::timeSystem() - startIter > timeLimit ) {
			std::cout << "reach the time limit (sec): "<< timeLimit  << ", quit!\n";
			break;
		}
	}
	std::cout << "Done iterations!\n";
//	std::cout<<"nodes explored: " << nodes << "\n";
	std::cout << "Tree size: " << GSIZE <<"\n";
	std::cout.precision(10);
	std::cout<<"["<<mex::timeSystem()-startIter<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<"\n";
//	exploredRatio(wmbUB);
}

int wmbsearch::expandBestNode() {
/*
 * forward pass to expand the best node in OPEN
 */
	if (verbose>2) {
		std::cout << "expandBestNode: start running..." << std::endl;
	}
	mex::vector<uint32_t> tuple(nvar);
	auto ptr = exploredTree.begin();
	auto ptrPar = ptr;
	// forward pass, also update topPath_ub, topPath_lb
	while ( (*ptr).toChild != NULL ) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		ptrPar = ptr;
		ptr = (*ptr).toChild;

		if (  (*ptr).type == AND_NODE ) {
			// ptr is AND_NODE
			// ptrPar is OR_NODE
			// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
			(*ptr).topPath_ub = (*ptrPar).topPath_ub;
			(*ptr).topPath_lb = (*ptrPar).topPath_lb;

			tuple[(*ptr).X] = (*ptr).val;
//			// current tub, tlb are not part of topPath_ub or toPath_lb
		} else {
			// ptr is OR_NODE
			// ptrPar is AND_NODE
			// use the best bounds of siblings
			auto bd = getSibBounds(ptr);
			ptr->topPath_lb = ptrPar->topPath_lb + bd.first;
			ptr->topPath_ub = ptrPar->topPath_ub + bd.second;
		}
	}
	if (verbose>2) {
		std::cout << "expandBestNode: best child found!" << std::endl;
		printNode(*ptr);
	}
	//
//	ptr->solved = checkSolved(*ptr);
	if ( ptr->solved  ) {
		std::cout << "The best frontier node has been solved, quit!" << std::endl;
		return 3;
	}
	// now an OR node can be frontier node when all its children have been deleted before
	if (ptr->type == AND_NODE) {
		double cur_lb=0.0, cur_ub = 0.0;
		std::list<mex::graphModel::vindex> childrenList;
		if ( (*ptr).id  > -1 ) {
			childrenList = ascList[(*ptr).id];
		}
		else {
			// if it is the virtual root
			childrenList = rts;
		}
		// if the most promising node is a node without children,
		// the whole search progress should be terminated, although it's possible that the bound gap remains non-zero.
		// it depends on what priority you use
		if (childrenList.empty()) {
			std::cout <<"The best node: variable "<< (*ptr).id << " with value " << (*ptr).val << " has no children, quit!\n";
			printNode(*ptr);
			if (verbose > 1) {
				while ( (*ptr).id > -1 ) {
					printNode(*ptr);
					ptr = exploredTree.parent(ptr);
				}
				printNode(*ptr);
			}
			return 2;
		}
		// update bounds later
		bool isProp = false;
//	(*ptr).ub = (*ptr).lb = (*ptr).exact_val;
//		if (ptr->isExpanded) {
		if (ptr->isReExpand) {
			// if expanded before, reset
			// also no need to propagate bounds
			ptr->exact_val = wmbUB.heuristicTheta(ptr->X, tuple);
//			isProp = true;
		}
//		ptr->cur_lb = ptr->cur_ub = ptr->exact_val;
		cur_lb = cur_ub = ptr->exact_val;
		//
		for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
			searchnode child;
			child.type = OR_NODE;
			child.id = *it;
			child.X = wmbUB.var(child.id);
			// for OR_NODE node
			child.exact_val = -std::numeric_limits<double>::infinity();
			// append to the explored tree
			auto ptrChild = exploredTree.append_child(ptr, child);
			++GSIZE;
			//
			mex::Factor dUB(child.X, 0.0);
			mex::Factor dLB(child.X, 0.0);
			//
			int nMini = wmbUB.duplicateInBucket(child.X, tuple);
			isProp = isProp || (nMini > 1);
			if (verbose > 2) {
				if (nMini > 1) {
					std::cout
							<< "Bounds should be improved when expanding this node, no. of mini-buckets: "
							<< nMini << std::endl;
				}
			}
			for (size_t v = 0; v < (child.X).states(); ++v) {
				searchnode kid;
				kid.type = AND_NODE;
				kid.X = child.X;
				kid.id = child.id;
				kid.val = v;

				tuple[child.X] = v;
				// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
				// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
				kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
				// actually, we can speed up this evaluation by computing each factor only once
				// for a leaf node, leaf_ub, leaf_lb are its ub, lb
				kid.downWorstPath_lb = kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
				kid.downWorstPath_ub = kid.downPath_ub = wmbUB.heuristicIn(kid.X, tuple);

//				kid.cur_lb = kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
//				kid.cur_ub = kid.ub = dUB[v] = kid.downPath_ub + kid.exact_val;
				kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
				kid.ub = dUB[v] = kid.downPath_ub + kid.exact_val;
				// whether solved
				kid.solved = checkSolved(kid);
				// if solved, no need to append to the exploredTree
				if (kid.solved) {
					double mx = std::max((*ptrChild).exact_val, kid.ub);
					(*ptrChild).exact_val = mx
							+ log(
									exp((*ptrChild).exact_val - mx)
											+ exp(kid.ub - mx));
					continue;
				}
				// append
				exploredTree.append_child(ptrChild, kid);
				++GSIZE;
			}
			//
//			ptrChild->cur_lb = ptrChild->lb = dLB.logsumexp();
//			ptrChild->cur_ub = ptrChild->ub = dUB.logsumexp();

			ptrChild->lb = dLB.logsumexp();
			ptrChild->ub = dUB.logsumexp();

//		(*ptr).lb += (*ptrChild).lb;
//		(*ptr).ub += (*ptrChild).ub;
//			ptr->cur_lb += ptrChild->cur_lb;
//			ptr->cur_ub += ptrChild->cur_ub;
			cur_lb += ptrChild->lb;
			cur_ub += ptrChild->ub;
			// if no child left, safely delete this OR_NODE node
			if (exploredTree.begin(ptrChild) == exploredTree.end(ptrChild)) {
				(*ptr).exact_val += (*ptrChild).ub;
				exploredTree.erase(ptrChild);
				--GSIZE;
			}
		}
		//
//		if (!ptr->isReExpand) {
//			// if has not been expanded before
////			ptr->lb = ptr->cur_lb = cur_lb;
////			ptr->ub = ptr->cur_ub = cur_ub;
//			ptr->lb = std::max(ptr->lb, cur_lb);
//			ptr->ub = std::min(ptr->ub, cur_ub);
//		}
		// propagate bounds if necessary
		if (isProp && !ptr->isReExpand) {
//			ptr->lb = std::max(ptr->lb, cur_lb);
//			ptr->ub = std::min(ptr->ub, cur_ub);
			ptr->lb = cur_lb;
			ptr->ub = cur_ub;
			updateBounds(ptr);
		}
		// this node is solved only when all its children has been solved
		// note: none of the solved children has been appended
		ptr->solved = (exploredTree.begin(ptr)==exploredTree.end(ptr));
//		ptr->solved = (exploredTree.begin(ptr)==exploredTree.end(ptr)) || checkSolved(*ptr);
		if ((*ptr).solved) {
			return backwardUpdate(ptr);
		}
		// set best child here, local forward, backward pass
		double bestPrior = -std::numeric_limits<double>::infinity();
		double worstPrior = std::numeric_limits<double>::infinity();
		tree<searchnode>::iterator toChild = NULL;
		ptr->toWorst=NULL; // initialization   // TODO qlou was ==NULL?
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			tree<searchnode>::iterator toKid = NULL;
			double bestKidPrior = -std::numeric_limits<double>::infinity();
			double worstKidPrior = std::numeric_limits<double>::infinity();
			// tub is part of ub
//			sib->topPath_lb = ptr->topPath_lb + ptr->cur_lb - sib->cur_lb;
//			sib->topPath_ub = ptr->topPath_ub + ptr->cur_ub - sib->cur_ub;
			sib->topPath_lb = ptr->topPath_lb + cur_lb - sib->lb;
			sib->topPath_ub = ptr->topPath_ub + cur_ub - sib->ub;

			for (tree<searchnode>::sibling_iterator nb = exploredTree.begin(
					sib); nb != exploredTree.end(sib); ++nb) {
				// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
				(*nb).topPath_ub = (*sib).topPath_ub;
				(*nb).topPath_lb = (*sib).topPath_lb;
				// works for both the worst and the best
				double prior = calcFrontierPrior(*nb);
				//
				nb->priority = prior;
				//
				if ((toKid == NULL) || (prior > bestKidPrior) ) {
					toKid = nb;
					bestKidPrior = prior;
				}
//				if ( (sib->toWorst==NULL) || (prior < worstKidPrior) ) {
//					sib->toWorst = nb;
//					worstKidPrior = prior;
//				}
			}
			// since those children are all frontier nodes, sib's best = sib'worst since we delete all sib's children immediately
			worstKidPrior = bestKidPrior;
			sib->toWorst = toKid;
			//
			sib->toChild = toKid;
			// update with best child
			sib->downPath_lb = toKid->lb;
			sib->downPath_ub = toKid->ub;
			// update with best kid's priority
			sib->priority = bestKidPrior;
			//
			if ( (toChild == NULL) || (bestKidPrior > bestPrior) ) {
				toChild = sib;
				bestPrior = bestKidPrior;
			}
			// update with worst child
			sib->downWorstPath_lb = sib->toWorst->lb;
			sib->downWorstPath_ub = sib->toWorst->ub;
			if ( (ptr->toWorst == NULL) || (worstKidPrior < worstPrior) ) {
				ptr->toWorst = sib;
				worstPrior = worstKidPrior;
			}
		}
		// if not solved
		ptr->toChild = toChild;
		// take the min
		ptr->priority = std::min(bestPrior, ptr->priority); // initialized to be inf
		// not including exact_val
//	(*ptr).downPath_ub = (*ptr).ub - (*ptr).exact_val - (*toChild).ub;
//	(*ptr).downPath_lb = (*ptr).lb - (*ptr).exact_val - (*toChild).lb;
//	ptr->downPath_lb = toChild->downPath_lb + ptr->lb - ptr->exact_val - toChild->lb;
//	ptr->downPath_ub = toChild->downPath_ub + ptr->ub - ptr->exact_val - toChild->ub;

		ptr->downPath_lb = toChild->downPath_lb + cur_lb - ptr->exact_val - toChild->lb;
		ptr->downPath_ub = toChild->downPath_ub + cur_ub - ptr->exact_val - toChild->ub;
		//
		ptr->downWorstPath_lb = ptr->toWorst->downWorstPath_lb + cur_lb - ptr->exact_val - ptr->toWorst->lb;
		ptr->downWorstPath_ub = ptr->toWorst->downWorstPath_ub + cur_ub - ptr->exact_val - ptr->toWorst->ub;
		// now being expanded
//		ptr->isExpanded = true;
	} else {
		// OR node
		assert(ptr->isReExpand); // children must be deleted before for OR node

		auto X = ptr->X;
		ptr->exact_val = -std::numeric_limits<double>::infinity();
		ptr->toChild = ptr->toWorst = NULL;
//		mex::Factor dUB(X, 0.0);
//		mex::Factor dLB(X, 0.0);

		double bestPrior = -std::numeric_limits<double>::infinity();
//		double worstPrior = std::numeric_limits<double>::infinity();
		double prior = 0.0;

		for (size_t v = 0; v < X.states(); ++v) {
			// do it first
			tuple[X] = v;

			searchnode kid;
			kid.type = AND_NODE;
			kid.X = X;
			kid.id = ptr->id;
			kid.val = v;
			// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
			// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
			kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
			// actually, we can speed up this evaluation by computing each factor only once
			kid.downWorstPath_lb = kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
			kid.downWorstPath_ub = kid.downPath_ub = wmbUB.heuristicIn(kid.X, tuple);

//			kid.cur_lb = kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
//			kid.cur_ub = kid.ub = dUB[v] = kid.downPath_ub + kid.exact_val;
			kid.lb = kid.downPath_lb + kid.exact_val;
			kid.ub = kid.downPath_ub + kid.exact_val;
			// whether solved
			kid.solved = checkSolved(kid);
			// if solved, no need to append to the exploredTree
			if (kid.solved) {
				double mx = std::max(ptr->exact_val, kid.ub);
				ptr->exact_val = mx + log( exp(ptr->exact_val - mx) + exp(kid.ub - mx));
				continue;
			}
			// assign path values
			kid.topPath_lb = ptr->topPath_lb;
			kid.topPath_ub = ptr->topPath_ub;
			// append
			auto ptrKid = exploredTree.append_child(ptr, kid);
			++GSIZE;
			// work both for the worst and the best
			prior = calcFrontierPrior(kid);
			//
			kid.priority = prior;

			if ( (ptr->toChild==NULL) || (prior>bestPrior) ) {
				ptr->toChild = ptrKid;
				bestPrior = prior;
			}
//			if ( (ptr->toWorst==NULL) || (prior<worstPrior) ) {
//				ptr->toWorst = ptrKid;
//				worstPrior = prior;
//			}
		}
		//
		ptr->priority = std::min(bestPrior,  ptr->priority);

		// since we delete OR node's children at the same time set toWorst = toBest
		ptr->toWorst = ptr->toChild;
		// check solved, actually it is not quite possible
		ptr->solved = (exploredTree.begin(ptr) == exploredTree.end(ptr));
		if (!ptr->solved) {
			ptr->downPath_lb = ptr->toChild->lb;
			ptr->downPath_ub = ptr->toChild->ub;

			ptr->downWorstPath_lb = ptr->toWorst->lb;
			ptr->downWorstPath_ub = ptr->toWorst->ub;
		}
//		cur_lb = dLB.logsumexp();
//		cur_ub = dUB.logsumexp();
//		// update current bounds if necessary
//		// ptr->cur_lb, ptr->cur_ub cannot be worse than cur_lb, cur_ub
//		if ( (ptr->cur_lb > cur_lb + solvedThresh) || (ptr->cur_ub < cur_ub - solvedThresh) ) {
//			ptr->cur_lb = cur_lb;
//			ptr->cur_ub = cur_ub;
//			updateBounds(ptr);
//		}
	}
	if (verbose>2) {
		std::cout << "the best frontier node after expansion" << std::endl;
		printNode(*ptr);
	}
	return backwardUpdate(ptr);
}
int wmbsearch::backwardUpdate(tree<searchnode>::iterator ptr) {
/*
 * backward update after expanding the best or removing the worst frontier node
 */
	if (verbose > 2) {
		std::cout << "running backwardUpdate..." << std::endl;
	}

	if (ptr == NULL) {
		std::cout << "backwardUpdate: NULL occurs, problematic or in purpose!" << std::endl;
		return -1;
	}
	// double check,
	// it's possible to be OR node
	// the returned ptr is not solved
	// so if run this function after DFS, we will get input ptr's ancestor
	// doPruning will possibly make the returned node points to an arbitrary child.
	// no worries, we will re-calculate the best child here
	// caveat: make sure bounds have been fully updated before pruning
	ptr = prune(ptr);
	// check whether the whole tree has been solved or not
	auto rt = exploredTree.begin();
	if ( rt->solved ) {
		std::cout << "The root has been solved!\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		return 1;
	}
	// after pruning, it's possible to be either AND_NODE or OR_NODE
	//
	// no need to do pruning in the while loop
	// also, no need to update "solved" in the while loop because if any ancestor is solved,
	// it should be already deleted in the pruning step.
	std::pair<double,double> bd;
	while (true) {
		// initialize with the smallest value
		double bestPrior = -std::numeric_limits<double>::infinity();
		double worstPrior = std::numeric_limits<double>::infinity();
		double bePrior = 0.0;
		double woPrior = 0.0;
//		assert(ptr->toChild != NULL);
//		if (ptr->isReExpand) {
//			ptr->priority = std::min( ptr->priority, calcPriority(*ptr) );
//		}
		// set to be NULL
		ptr->toChild = NULL;
		ptr->toWorst = NULL;
		if ((*ptr).type == AND_NODE) {
			bd = getChildBounds(ptr);
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				// AND_NODE node, exact_val is part of ub, lb
				// use the best bounds to set priority
				sib->topPath_lb = ptr->topPath_lb + bd.first - sib->lb;
				sib->topPath_ub = ptr->topPath_ub + bd.second - sib->ub;
				bePrior = calcPriority(*sib);
				//
				if (bePrior > sib->priority){
					bePrior = sib->priority;
				} else {
					sib->priority = bePrior;
				}
				//
				woPrior = calcWorstPriority(*sib);
				//
				woPrior = std::min(woPrior, sib->priority);
				//
				if ( (ptr->toChild==NULL) || (bePrior > bestPrior) ) {
					bestPrior = bePrior;
					ptr->toChild = sib;
				}
				if ( (ptr->toWorst==NULL) || (woPrior < worstPrior) ) {
					// if sib is OR node and has no children, should not be considered
					if (exploredTree.number_of_children(sib) < 1) continue;
					worstPrior = woPrior;
					ptr->toWorst = sib;
				}
			}
			// update with the best child
			auto sib = ptr->toChild;
			if (sib != NULL) {
				ptr->downPath_lb = sib->downPath_lb + bd.first - ptr->exact_val - sib->lb;
				ptr->downPath_ub = sib->downPath_ub + bd.second - ptr->exact_val - sib->ub;
				//
				ptr->priority = std::min(ptr->priority,  sib->priority);
			}
			sib = ptr->toWorst;
			if (sib != NULL) {
				ptr->downWorstPath_lb = sib->downWorstPath_lb + bd.first - ptr->exact_val - sib->lb;
				ptr->downWorstPath_ub = sib->downWorstPath_ub + bd.second - ptr->exact_val - sib->ub;
			} else {
				// it is possible that none of ptr's children have any children
				ptr->toWorst = ptr->toChild;
				ptr->downWorstPath_lb = ptr->downPath_lb;
				ptr->downWorstPath_ub = ptr->downPath_ub;
			}
		}
		else {
			// ptr is OR node
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				// ptr is OR node
				// use most update-to-date information, topPath for AND_NODE siblings are the same
				// equal to those of the OR_NODE parent
				(*sib).topPath_ub = (*ptr).topPath_ub;
				(*sib).topPath_lb = (*ptr).topPath_lb;
				bePrior = calcPriority(*sib);
				//
				if (bePrior > sib->priority){
					bePrior = sib->priority;
				} else {
					sib->priority = bePrior;
				}
				//
				woPrior = calcWorstPriority(*sib);
				// pick the smaller one
				woPrior = std::min(woPrior, sib->priority);
				//
				if ( (ptr->toChild==NULL) || (bePrior > bestPrior) ) {
					bestPrior = bePrior;
					(*ptr).toChild = sib;
				}
				if ( (ptr->toWorst==NULL) || (woPrior < worstPrior) ) {
					// if sib is AND node and has no children, should not be considered
					if (exploredTree.number_of_children(sib) < 1) continue;
					worstPrior = woPrior;
					ptr->toWorst = sib;
				}
			}
			// update with the best child
			auto sib = ptr->toChild;
			if ( sib != NULL ) {
				if ((*sib).toChild != NULL) {
					(*ptr).downPath_lb = (*sib).exact_val + (*sib).downPath_lb;
					(*ptr).downPath_ub = (*sib).exact_val + (*sib).downPath_ub;
				} else {
					// if sib is a frontier AND node
					ptr->downPath_lb = sib->lb;
					ptr->downPath_ub = sib->ub;
				}
				//
				ptr->priority = std::min(ptr->priority,  sib->priority);
			}
			// update with the worst child
			sib = ptr->toWorst;
			if ( sib != NULL ) {
				if ((*sib).toChild != NULL) {
					(*ptr).downWorstPath_lb = (*sib).exact_val + (*sib).downWorstPath_lb;
					(*ptr).downWorstPath_ub = (*sib).exact_val + (*sib).downWorstPath_ub;
				} else {
					std::runtime_error("sib without children should not be considered!");
					// if sib is a frontier AND node
					ptr->downWorstPath_lb = sib->lb;
					ptr->downWorstPath_ub = sib->ub;
				}
			}
			else {
				// it is possible that none of ptr's children have any children
				ptr->toWorst = ptr->toChild;
				ptr->downWorstPath_lb = ptr->downPath_lb;
				ptr->downWorstPath_ub = ptr->downPath_ub;
			}
		}
		// go upward
		if ( ptr == rt ) {
			break;
		}
		ptr = exploredTree.parent(ptr);
	}
	return 0;
}
std::pair<double,double> wmbsearch::getSibBounds(tree<searchnode>::iterator ptr) {
/*
 * accumulate bounds from sublings (deleted or non-deleted), edge weight also included
 * if getAll is true, return all bounds including ptr itself, default, false
 * use the best bounds of its sublings
 * currently only for OR node
 */
	assert( ptr->type == OR_NODE);
	auto ptrPar = exploredTree.parent(ptr);
	double lb = ptrPar->exact_val;
	double ub = lb;
	for (auto sib=exploredTree.begin(ptrPar); sib != exploredTree.end(ptrPar); ++sib) {
		if (sib != ptr ) {
			lb += sib->lb;
			ub += sib->ub;
		}
	}
	return std::make_pair(lb,ub);
}
std::pair<double,double> wmbsearch::getChildBounds(tree<searchnode>::iterator ptr) {
/*
 * accumulate best bounds of it's children (deleted or non-deleted), edge weight also included
 * children's best bound may not correspond to the parent's best bounds
 * currently only for AND node
 */
//	assert( ptr!=exploredTree.begin());
	assert( ptr->type == AND_NODE);
	double lb = ptr->exact_val;
	double ub = lb;
	for (auto sib=exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		lb += sib->lb;
		ub += sib->ub;
	}
	return std::make_pair(lb,ub);
}
void wmbsearch::updateBounds(tree<searchnode>::iterator ptr) {
	// re-calculate bounds from existing children and exact_val
	if (verbose > 2) {
		std::cout << "now update bounds: " << std::endl;
	}
	assert( ptr != NULL );
	auto rt = exploredTree.begin();
	double cur_lb, cur_ub;
	while( ptr != rt ) {
		ptr = exploredTree.parent(ptr);
		// update bounds by re-calculation.
		// more time-consuming than simple update, but more numerically stable
		if ((*ptr).type == AND_NODE) {
//			(*ptr).cur_lb = (*ptr).exact_val;
//			(*ptr).cur_ub = (*ptr).exact_val;
			cur_lb = (*ptr).exact_val;
			cur_ub = (*ptr).exact_val;
			// possibly some children have been pruned
			for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//				(*ptr).cur_lb += (*sib).cur_lb;
//				(*ptr).cur_ub += (*sib).cur_ub;
				cur_lb += (*sib).lb;
				cur_ub += (*sib).ub;
			}
		}
		else {
			// OR_NODE node
			std::vector<double> ub_vec {(*ptr).exact_val}, lb_vec {(*ptr).exact_val};
			for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//				lb_vec.push_back( (*sib).cur_lb );
//				ub_vec.push_back( (*sib).cur_ub );
				lb_vec.push_back( (*sib).lb );
				ub_vec.push_back( (*sib).ub );
			}
//			(*ptr).cur_lb = logsumexp(lb_vec);
//			(*ptr).cur_ub = logsumexp(ub_vec);
			cur_lb = logsumexp(lb_vec);
			cur_ub = logsumexp(ub_vec);
		}
//		if ( ptr->cur_lb > ptr->lb || ptr->cur_ub < ptr->ub ) {
//			ptr->lb = ptr->cur_lb;
//			ptr->ub = ptr->cur_ub;
//		}
		if ( cur_lb > ptr->lb || cur_ub < ptr->ub ) {
			ptr->lb = cur_lb;
			ptr->ub = cur_ub;
		} else {
			// no need to move upward
			break;
		}
	}
	// update bounds for the root
	LB = rt->lb;
	UB = rt->ub;
//
//	if (UB-LB < boundGap) {
//		std::cout.precision(10);
//		std::cout << "[" << mex::timeSystem() - startIter << "]: "
//				<< std::setw(10) << LB << " < ln Z < " << UB << ", UB-LB = " << UB-LB << " < " << boundGap << "\n";
//		std::cout << "Tree size: " << GSIZE <<"\n";
//		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//		boundGap = boundGap/10.0;
//	}
}
tree<wmbsearch::searchnode>::iterator wmbsearch::prune(tree<searchnode>::iterator ptr) {
/* do bottom-up pruning if necessary,
 * make sure do it from the frontiers to avoid memory issue
 * always be very careful about deleting nodes, make sure those node pointers in the stack will not be wide pointers
 * ancestor is the upper limit node. deletion stops once we reach this node even if it is solved, default: the root
 * ancestor must be an ancestor on the path
 * return the highest not being solved
 * currently, it should only apply to DFS related procedure
*/
	auto ancestor = exploredTree.begin();
	if (ptr==NULL) return ptr;
	if ( !(*ptr).solved || (ptr == ancestor) ) return ptr;
	tree<searchnode>::iterator highestDelete = ptr;
	//
	ptr = exploredTree.parent(ptr);

	while ( ptr != ancestor ) {
		bool solved = true;
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// a node is solved only when all its children are solved
			if( !(*sib).solved ) {
				solved = false;
				break;
			}
		}
		// caveat: unsolved children cannot be removed alone even in SMA*
//		solved = solved && checkSolved(*ptr);

		(*ptr).solved = solved;
		// a node is solved only when all its children are solved
		if (!(*ptr).solved) {
			break;
		}
		//
		highestDelete = ptr;
		ptr = exploredTree.parent(ptr);
	}
	// erase the whole subtree rooted at this node (including itself)
	// store exact value in its parent's exact_val
	// if not the root
	ptr = exploredTree.parent(highestDelete);
	if ((*ptr).type == AND_NODE) {
		(*ptr).exact_val += (*highestDelete).ub; // at this moment, ub=lb
	} else {
		double maxval = std::max((*ptr).exact_val, (*highestDelete).ub);
		(*ptr).exact_val = maxval
				+ log(
						exp((*ptr).exact_val - maxval)
								+ exp((*highestDelete).ub - maxval));
	}
	// the parent points to some child for safety
	if ((*ptr).toChild == highestDelete) {
		(*ptr).toChild = NULL;
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			if (sib != highestDelete) {
				// no need to pick the best child at this moment, we'll do it in backward pass
				(*ptr).toChild = sib;
				break;
			}
		}
	}
	if ((*ptr).toWorst == highestDelete) {
		(*ptr).toWorst = NULL;
		for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			if (sib != highestDelete) {
				// no need to pick the worst child at this moment, we'll do it in backward pass
				(*ptr).toWorst = sib;
				break;
			}
		}
	}

	// finally delete
	int GSIZE_pre = GSIZE;
	// size gives you the size of the subtree rooted at given node ( including that node )
	GSIZE -= exploredTree.size(highestDelete);
	exploredTree.erase(highestDelete);

	if (verbose > 2) {
		std::cout << "pruning done: " << GSIZE_pre - GSIZE << " nodes deleted!\n" ;
	}
	// return parent of the highestDelete
	return ptr;
}

int wmbsearch::removeWorstNode() {
/*
 * remove the worst node
 * tricky, be careful
 */
	if (verbose > 2) {
		std::cout << "removeWorstNode: starting running..." << std::endl;
	}
	mex::vector<uint32_t> tuple(nvar);
	auto ptr = exploredTree.begin();
	// forward pass, also update topPath_ub, topPath_lb
	auto ptrPar = ptr;
	while ( (*ptr).toWorst != NULL ) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		ptrPar = ptr;
		ptr = ptrPar->toWorst;

		if (  (*ptr).type == AND_NODE ) {
			// ptr is AND_NODE
			// ptrPar is OR_NODE
			// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
			//			// current tub, tlb are not part of topPath_ub or toPath_lb
			(*ptr).topPath_ub = (*ptrPar).topPath_ub;
			(*ptr).topPath_lb = (*ptrPar).topPath_lb;

			tuple[(*ptr).X] = (*ptr).val;
		} else {
			// ptr is OR_NODE
			// ptrPar is AND_NODE
			// use the best bounds of siblings
			auto bd = getSibBounds(ptr);
			ptr->topPath_lb = ptrPar->topPath_lb + bd.first;
			ptr->topPath_ub = ptrPar->topPath_ub + bd.second;
		}
	}
	if (verbose>2) {
		std::cout << "the worst frontier node found" << std::endl;
		printNode(*ptr);
	}
//	ptr->solved = checkSolved(*ptr);
	if ( ptr->solved && (verbose > 2) ) {
		std::cout << "Hmm, the worst node is solved!" << std::endl;
//		return 0;
	}
	// can be AND or OR node
	if (ptr == exploredTree.begin()) {
		std::cout << "The worst node is the root, quit!" << std::endl;
		return 3; // rare case when tiny memory
	}
	ptrPar = exploredTree.parent(ptr);
	if (verbose > 2) {
		std::cout << "No. of children of the worst node's parent is " << exploredTree.number_of_children(ptrPar) << std::endl;
	}
	if (ptr->type == AND_NODE) {
//		--GSIZE;
//		if (exploredTree.number_of_children(ptrPar) > 1) {
//			// more than 1 child left, easy
////			if (ptrPar->toChild == ptr) {
////				ptrPar->toChild = NULL;
////			}
//			exploredTree.erase(ptr);
//			// no need to assign toChild and toWorst here, will be updated in backwardUpdate
//			// for safety, pick random one
//			ptrPar->toChild = ptrPar->toWorst = exploredTree.begin(ptrPar);
//		} else {
//			// if only one child left
//			ptrPar->downWorstPath_lb = ptrPar->downPath_lb = ptr->lb;
//			ptrPar->downWorstPath_ub = ptrPar->downPath_ub = ptr->ub;
//			ptrPar->toChild = ptrPar->toWorst = NULL;
//			exploredTree.erase(ptr);
//		}
	// remove all its siblings at the same time
		GSIZE -= (exploredTree.size(ptrPar)-1);
		exploredTree.erase_children(ptrPar);
		ptrPar->toChild = ptrPar->toWorst = NULL;

	} else {
		// ptr is OR node, delete all its sublings and their descendants
		// size gives you the size of the subtree rooted at given node ( including that node )
		GSIZE -= (exploredTree.size(ptrPar)-1);
		exploredTree.erase_children(ptrPar);
		ptrPar->toChild = ptrPar->toWorst = NULL;
	}
	// set to be re-expanded
	ptrPar->isReExpand = true;
	if (verbose>2) {
		std::cout << "removeWorstNode: The worst frontier node removed, print its parent before backtrack:"<< std::endl;
		printNode(*ptrPar);
	}
	return backwardUpdate(ptrPar);
}
double wmbsearch::calcPriority(const searchnode& n) {
	// the priority function is non-decreasing from child to parent when treated as a frontier node.
	// for a frontier OR node, should calculate its best child
//	assert( n.id > -1 );
	if ( n.solved ) {
		return -std::numeric_limits<double>::infinity();
	}
	//
	double prior = 0.0, ub = 0.0, lb = 0.0;
	if (priorityName == "absolute_gap") {
		// absolute gap between upper and lower bounds, i.e., "upper_minus_lower"
//		if (n.toChild != NULL) {
//			if ( n.type == AND_NODE ) {
//				lb = n.downPath_lb + n.topPath_lb + n.exact_val;
//				ub = n.downPath_ub + n.topPath_ub + n.exact_val;
//			} else {
//				// for OR_NODE node, exact_val = -inf for initialization
//				lb = n.downPath_lb + n.topPath_lb;
//				ub = n.downPath_ub + n.topPath_ub;
//			}
//		} else {
//			// if it is a leaf AND_NODE node in OPEN
//			// for frontier OR_NODE node, existing in the follow-up DFS, the following formula are more general
//			// note that exact_val is part of n.ub, n.lb
//
//			ub = n.ub + n.topPath_ub;
//			lb = n.lb + n.topPath_lb;
//
//		}
		if (n.type == AND_NODE) {
			if (n.toChild != NULL) {
				lb = n.downPath_lb + n.topPath_lb + n.exact_val;
				ub = n.downPath_ub + n.topPath_ub + n.exact_val;
			} else {
				lb = n.lb + n.topPath_lb;
				ub = n.ub + n.topPath_ub;
			}
		} else {
			// for OR_NODE node, exact_val = -inf for initialization
			// also work for frontier OR_NODE node,
			lb = n.downPath_lb + n.topPath_lb;
			ub = n.downPath_ub + n.topPath_ub;
		}

		prior = lb - ub;
		prior = (prior < 0.0) ? prior : 0.0;
		prior = ub + std::log(1 - std::exp(prior));
	}
	else if (priorityName == "absolute_upper") {
		// absolute upper bounds
//		if (n.toChild != NULL) {
////			ub = n.leaf_ub + n.downPath_ub + n.topPath_ub + n.tub;
//			if ( n.type == AND_NODE ) {
//				ub = n.downPath_ub + n.topPath_ub + n.exact_val;
//			} else {
//				ub = n.downPath_ub + n.topPath_ub;
//			}
//
//		}
//		else {
////			ub = n.leaf_ub + n.topPath_ub;
//			ub = n.ub + n.topPath_ub;
//		}
		if ( n.type == AND_NODE ) {
			if (n.toChild != NULL) {
				ub = n.downPath_ub + n.topPath_ub + n.exact_val;
			} else {
				ub = n.ub + n.topPath_ub;
			}
		}
		else {
//			ub = n.leaf_ub + n.topPath_ub;
			// even work for frontier OR node
			ub = n.downPath_ub + n.topPath_ub;
		}
		prior = ub;
	}
	// do not support those priority types any  more
	else {
		throw std::runtime_error("no such priorityName defined!");
	}
	if (std::isnan(prior)) {
		std::cout << "prior = " << prior << "\n";
		std::cout <<"variable "<< n.id << " with value " << n.val << ", solved = " << n.solved << "\n";
		std::cout << n.ub << ", " << n.topPath_ub << ", " << n.downPath_ub << "\n";
		std::cout << n.lb << ", " << n.topPath_lb << ", " << n.downPath_lb << "\n";
		throw std::runtime_error("priority is not proper!");
	}
	// to make sure that an AND node's priority is always larger than its descendants
	// for OR nodes, we actually calculates its best child's priority
	if ((n.type == AND_NODE) && (n.toChild!=NULL) ) {
		// pick the min one
		prior = std::min(prior, calcFrontierPrior(n));
	}
	return prior;
}
double wmbsearch::calcFrontierPrior(const searchnode& n) {
/*
 * treat n as if it is a leaf node, only AND_NODE can be frontier node
 */
	assert(n.type == AND_NODE);
	if ( n.solved ) {
		return -std::numeric_limits<double>::infinity();
	}
	double prior = 0.0, ub = 0.0, lb = 0.0;
	if (priorityName == "absolute_gap") {
		// absolute gap between upper and lower bounds, i.e., "upper_minus_lower"
		lb = n.lb + n.topPath_lb;
		ub = n.ub + n.topPath_ub;
		prior = lb - ub;
		prior = (prior < 0.0) ? prior : 0.0;
		prior = ub + std::log(1 - std::exp(prior));
	}
	else if (priorityName == "absolute_upper") {
		// absolute upper bounds
		ub = n.ub + n.topPath_ub;
		prior = ub;
	}
	else {
		throw std::runtime_error("no such priorityName defined!");
	}
	if (std::isnan(prior)) {
		std::cout << "prior = " << prior << "\n";
		std::cout <<"variable "<< n.id << " with value " << n.val << ", solved = " << n.solved << "\n";
		std::cout << n.ub << ", " << n.topPath_ub << ", " << n.downPath_ub << "\n";
		std::cout << n.lb << ", " << n.topPath_lb << ", " << n.downPath_lb << "\n";
		throw std::runtime_error("priority is not proper!");
	}
	return prior;
}
double wmbsearch::calcWorstPriority(const searchnode& n) {
	// the priority function is non-decreasing from child to parent when treated as a frontier node.
	// for a frontier OR node, should calculate its best child
//	assert( n.id > -1 );
	if ( n.solved ) {
		return -std::numeric_limits<double>::infinity();
	}
	//
	double prior = 0.0, ub = 0.0, lb = 0.0;
	if (priorityName == "absolute_gap") {
		// absolute gap between upper and lower bounds, i.e., "upper_minus_lower"
		if (n.type == AND_NODE) {
			if (n.toChild != NULL) {
				lb = n.downWorstPath_lb + n.topPath_lb + n.exact_val;
				ub = n.downWorstPath_ub + n.topPath_ub + n.exact_val;
			} else {
				lb = n.lb + n.topPath_lb;
				ub = n.ub + n.topPath_ub;
			}
		} else {
			// for OR_NODE node, exact_val = -inf for initialization
			// also work for frontier OR_NODE node,
			lb = n.downWorstPath_lb + n.topPath_lb;
			ub = n.downWorstPath_ub + n.topPath_ub;
		}

		prior = lb - ub;
		prior = (prior < 0.0) ? prior : 0.0;
		prior = ub + std::log(1 - std::exp(prior));
	}
	else if (priorityName == "absolute_upper") {
		// absolute upper bounds
		if ( n.type == AND_NODE ) {
			if (n.toChild != NULL) {
				ub = n.downWorstPath_ub + n.topPath_ub + n.exact_val;
			} else {
				ub = n.ub + n.topPath_ub;
			}
		}
		else {
//			ub = n.leaf_ub + n.topPath_ub;
			// even work for frontier OR node
			ub = n.downWorstPath_ub + n.topPath_ub;
		}
		prior = ub;
	}
	// do not support those priority types any  more
	else {
		throw std::runtime_error("no such priorityName defined!");
	}
	if (std::isnan(prior)) {
		std::cout << "prior = " << prior << "\n";
		std::cout <<"variable "<< n.id << " with value " << n.val << ", solved = " << n.solved << "\n";
		std::cout << n.ub << ", " << n.topPath_ub << ", " << n.downPath_ub << "\n";
		std::cout << n.lb << ", " << n.topPath_lb << ", " << n.downPath_lb << "\n";
		throw std::runtime_error("priority is not proper!");
	}
	// to make sure that an AND node's priority is always larger than its descendants
	// for OR nodes, we actually calculates its best child's priority
	if ( (n.type==AND_NODE) && (n.toChild!=NULL) ) {
		// pick the min one
		prior = std::min(prior, calcFrontierPrior(n));
	}
	return prior;
}
int wmbsearch::DFS_Beta(double startTime, double& outFrequency) {
	/* upgrade DFS to serve as a baseline
	 * do DFS when the memory budget is (almost) used up
	 *	DFS implemented by stack, O(branch_factor * max_depth) memory required
	 *	DFS may be implemented s.t. the memory usage is only linear to max_depth
	 *	the order of the children is sorted (ascending order) by priority when pushed into the stack
	 */
	if (verbose>2) {
		std::cout << "DFS: start running..." << std::endl;
	}
	mex::vector<uint32_t> tuple(nvar);
	auto ptr = exploredTree.begin();
	auto ptrPar = ptr;
	// forward pass, also update topPath_ub, topPath_lb
	while ( (*ptr).toChild != NULL ) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		ptrPar = ptr;
		ptr = (*ptr).toChild;

		if (  (*ptr).type == AND_NODE ) {
			// ptr is AND_NODE
			// ptrPar is OR_NODE
			// for the AND_NODE node, topPath bounds are the same as those of its OR_NODE parent
			(*ptr).topPath_ub = (*ptrPar).topPath_ub;
			(*ptr).topPath_lb = (*ptrPar).topPath_lb;

			tuple[(*ptr).X] = (*ptr).val;
//			// current tub, tlb are not part of topPath_ub or toPath_lb
		} else {
			// ptr is OR_NODE
			// ptrPar is AND_NODE
			// use the best bounds of siblings
			auto bd = getSibBounds(ptr);
			ptr->topPath_lb = ptrPar->topPath_lb + bd.first;
			ptr->topPath_ub = ptrPar->topPath_ub + bd.second;
		}
	}
	if (verbose>2) {
		std::cout << "DFS: best child found!" << std::endl;
		printNode(*ptr);
	}
	//
//	ptr->solved = checkSolved(*ptr);
	if ( ptr->solved  ) {
		std::cout << "DFS: the best node has been solved, quit!" << std::endl;
		return 4;
	}
	// push current best leaf node to the stack
	// this node won't be deleted during pruning inside while loop
	std::stack<tree<searchnode>::iterator> stack;
	assert(stack.empty());
	stack.push(ptr);
	std::list<mex::graphModel::vindex> childrenList;
	auto ancestor = ptr;
	// we push both AND_NODE, OR_NODE nodes to the stack
	while (!stack.empty()) {
		ptr = stack.top();
		stack.pop();

		if (verbose > 2) {
			std::cout.precision(10);
			std::cout << "DFS before bound propagation: " << std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree Size: " << GSIZE << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}

		if ((*ptr).type == AND_NODE) {
			// AND_NODE nodes, expand its OR_NODE children, update bounds
			if ((*ptr).id > -1) {
				childrenList = ascList[(*ptr).id];
			} else {
				// if it is the virtual root
				childrenList = rts;
			}

			if ( (*ptr).solved ) {
				doPruning(ptr, ancestor);
				continue;
			}
			// add to tuple, be careful
			if ((*ptr).id > -1) tuple[ (*ptr).X ] = (*ptr).val;
			// update bounds later
			(*ptr).lb = (*ptr).exact_val;
			(*ptr).ub = (*ptr).exact_val;
			// whether propagate
			bool isProp = false;
			for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
				searchnode child;
				child.type = OR_NODE;
				child.id = *it;
				child.X = wmbUB.var(child.id);
				child.exact_val = -std::numeric_limits<double>::infinity();
				//
				mex::Factor dUB(child.X, 0.0);
				mex::Factor dLB(child.X, 0.0);
				//
//				mex::vector<uint32_t> childTuple(tuple);
//				isProp = isProp || (wmbUB.duplicateInBucket(child.X, childTuple) > 1);
				isProp = isProp || (wmbUB.duplicateInBucket(child.X, tuple) > 1);

				for (size_t v = 0; v < (child.X).states(); ++v) {
					searchnode kid;
					kid.type = AND_NODE;
					kid.X = child.X;
					kid.id = child.id;
					kid.val = v;

//					childTuple[child.X] = v;
					tuple[child.X] = v;
//					kid.exact_val = wmbUB.heuristicTheta(kid.X, childTuple);
					kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
					// actually, we can speed up this evaluation by computing each factor only once
//					kid.downPath_lb = wmbLB.heuristicIn(kid.X, childTuple);
//					kid.downPath_ub = wmbUB.heuristicIn(kid.X, childTuple);
					kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
					kid.downPath_ub = wmbUB.heuristicIn(kid.X, tuple);

					kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
					kid.ub = dUB[v] = kid.downPath_ub + kid.exact_val;
				}

				// OR_NODE node as a leaf node
				child.downPath_lb = child.lb = dLB.logsumexp();
				child.downPath_ub = child.ub = dUB.logsumexp();

				(*ptr).lb += child.lb;
				(*ptr).ub += child.ub;

				child.solved = checkSolved(child);
				if (child.solved) {
					(*ptr).exact_val += child.ub; // child.ub = child.lb
					continue;
				}
				// append to the explored tree if not solved yet
				// no solved nodes will be pushed into the stack
				exploredTree.append_child(ptr, child);
				++ GSIZE;
			}
			// propagate bounds to the root before pruning,
			// only need for AND_NODE
			// propagate bounds only when necessary
			if (isProp) {
				propagateBounds(ptr);
				if (verbose > 2) {
					std::cout.precision(10);
					std::cout << "DFS after bound propagation: " << std::setw(10) << LB << " < ln Z < " << UB << "\n";
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				}
			}
			// if no children appended, node is solved, do pruning
			if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
				(*ptr).solved = true;
				doPruning(ptr, ancestor);
				continue;
			}

			// calculate the priority
			// sort the OR_NODE children based on priority
			// all those nodes are not yet solved
			std::vector<double> childrenPriorities;
			std::vector<tree<searchnode>::iterator> childrenPtrs;
//			assert( childrenPriorities.empty() && childrenPtrs.empty() );
			for (tree<searchnode>::sibling_iterator sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				//
				(*sib).topPath_lb = (*ptr).topPath_lb + (*ptr).lb - (*sib).lb;
				(*sib).topPath_ub = (*ptr).topPath_ub + (*ptr).ub - (*sib).ub;

				childrenPriorities.push_back( calcPriority(sib) );
				childrenPtrs.push_back(sib);
			}
			// ascending order
			auto sortedIds = sortIndexes(childrenPriorities);
			// not quite necessary in DFS
			(*ptr).toChild = childrenPtrs[ sortedIds.back() ];
			// thus high-priority nodes will be popped out first
			for (auto id : sortedIds) {
				stack.push(childrenPtrs[id]);
			}
		}
		else {
			// OR_NODE node, re-expand its AND_NODE children
			// no need to update bounds because we've already done when generating the OR_NODE node
			std::vector<double> childrenPriorities;
			std::vector<tree<searchnode>::iterator> childrenPtrs;
//			assert( childrenPriorities.empty() && childrenPtrs.empty() );
			for (size_t v = 0; v < ((*ptr).X).states(); ++v) {
				searchnode kid;
				kid.type = AND_NODE;
				kid.X = (*ptr).X;
				kid.id = (*ptr).id;
				kid.val = v;

				tuple[(*ptr).X] = v;
				kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
				// actually, we can speed up this evaluation by computing each factor only once
				// for a leaf node, leaf_ub, leaf_lb are its ub, lb
				kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
				kid.downPath_ub = wmbUB.heuristicIn(kid.X, tuple);

				kid.lb = kid.downPath_lb + kid.exact_val;
				kid.ub = kid.downPath_ub + kid.exact_val;

				kid.topPath_ub = (*ptr).topPath_ub;
				kid.topPath_lb = (*ptr).topPath_lb;
				//
				kid.solved = checkSolved(kid);
				if (kid.solved) {
					// no need to update bounds here
					double mx = std::max( (*ptr).exact_val, kid.ub );
					(*ptr).exact_val = mx + log( exp( (*ptr).exact_val - mx ) + exp( kid.ub -mx ) );
					continue;
				}

				auto toKid = exploredTree.append_child(ptr, kid);
				++GSIZE;

				childrenPriorities.push_back( calcPriority(toKid) );
				childrenPtrs.push_back( toKid );
			}
			// no active children
			if (childrenPriorities.empty()) {
				(*ptr).solved = true;
				doPruning(ptr, ancestor);
				continue;
			}

			auto sortedIds = sortIndexes(childrenPriorities);
			for (auto id : sortedIds) {
				stack.push(childrenPtrs[id]);
			}
		}

		if ( (verbose>2) && (mex::timeSystem() - startIter > timeThresh) ) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			timeThresh += timeUnit*timeRatio;
			timeRatio *= timeRatio;
		}
		if (mex::timeSystem() - startTime > timeBudget) {
//			std::cout << "DFS runs out of time, quit!" << std::endl;
			std::cout << "Reach search time limit (sec): "<< timeBudget  << ", quit!\n";
			return 3;
		}
		if (mex::timeSystem() - startTime > outFrequency ) {
			outFrequency += outFrequency; // exponential
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startTime << "]: "
			<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree size: " << GSIZE <<std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}
	}
	// none of anscetor's descendants should exist
	assert(exploredTree.begin(ancestor) == exploredTree.end(ancestor));
	ancestor->solved = true;
	if (verbose>2) {
		std::cout << "DFS: best child solved!" << std::endl;
		printNode(*ancestor);
	}
	return bottomUpPass(ancestor);
}
void wmbsearch::startSearch(Mode mode){
//	std::cout << "Memory used for wmbLB + wmbUB  (MB): " << wmbLB.memory() + wmbUB.memory() << std::endl;
	double timeLimit = timeBudget;
	double memLimit = memoryBudget;
	double startTime = mex::timeSystem();
	std::string searchMode = (mode==Mode::SMA)? "SMA" : "DFS";
	buildPseudoTree();
	setRootConfig();
	std::cout << "Priority type: " << priorityName << std::endl;
	std::cout << "Search mode: " << searchMode << std::endl;
	int indicator = 0; // whether we have to terminate the program
	auto rt = exploredTree.begin();
	// in byte
	// caveat very coarse control, should be modified later
	int memUnit = sizeof(*rt);
	// MemLimit is in megabyte
	auto treeSizeLimit = unsigned(1024*1024*memLimit/memUnit);
	// debug
//	treeSizeLimit = 2e3;

	std::cout << "Search time limit (sec): " << timeLimit << "\n";
	std::cout << "Search memory limit (MB): " << memLimit << "\n";
	std::cout << "Tree size limit: " << treeSizeLimit << std::endl;

	double outFrequency = 0.1; // for output

	if (verbose > 2) {
		std::cout.precision(10);
		std::cout << "[" << mex::timeSystem() - startTime << "]: "
				<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
		std::cout << "Tree size: " << GSIZE << "\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
				<< std::endl;
	}

	bool isFirstReachMemLimit = false;
	while(true) {
		if ( rt->solved ) {
			std::cout << "The root is solved, quit!"<<std::endl;
			break;
		}
//		UB = rt->ub;
//		LB = rt->lb;
		if (mex::timeSystem() - startTime > outFrequency ) {
			outFrequency += outFrequency; // exponential
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startTime << "]: "
			<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
			std::cout << "Tree size: " << GSIZE <<std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}

		if ( GSIZE < treeSizeLimit) {
			// debug
			if (verbose > 2) {
				std::cout << "~~~~~~~~~~~~~~~~~before expanding the best frontier~~~~~~~~~~~~~~~~~\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
				<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE <<"\n";
				std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				std::cout << "root information:" << std::endl;
				printNode(*rt);
				std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			}
			indicator = expandBestNode();
			if (indicator != 0) {
//				std::cout<< "The best frontier has no children, loop terminated!\n";
				break;
			}
			if (verbose > 2) {
				std::cout << "~~~~~~~~~~~~~~~~~after expanding the best frontier~~~~~~~~~~~~~~~~~\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
				<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE <<"\n";

				// debug
//				if (GSIZE != exploredTree.size()) std::runtime_error("size not match");

				std::cout << "root information:" << std::endl;
				printNode(*rt);
				std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			}
		} else {
			// reach tree size limit (memory limit)
			if (!isFirstReachMemLimit) {
				// only prompt at the first reach memory limit
				std::cout << "First time reach tree size limit: " << treeSizeLimit << " ( memory limit: "<< memLimit << " MB )\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE <<"\n";
				std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				isFirstReachMemLimit = true;
			}

			if (mode == Mode::SMA) {
				if (verbose > 2) {
					std::cout << "~~~~~~~~~~~~~~~~~before removing the worst frontier~~~~~~~~~~~~~~~~~\n";
					std::cout.precision(10);
					std::cout << "[" << mex::timeSystem() - startTime << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
					std::cout << "Tree size: " << GSIZE <<"\n";
					std::cout << "root information:" << std::endl;
					printNode(*rt);
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				}
				indicator = removeWorstNode();
				if (verbose > 2) {
					std::cout << "~~~~~~~~~~~~~~~~~after removing the worst frontier~~~~~~~~~~~~~~~~~\n";
					std::cout.precision(10);
					std::cout << "[" << mex::timeSystem() - startTime << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
					std::cout << "Tree size: " << GSIZE <<"\n";
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
					std::cout << "root information:" << std::endl;
					printNode(*rt);
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				}
			}
			else if (mode == Mode::DFS) {
				if (verbose > 2) {
					std::cout << "~~~~~~~~~~~~~~~~~before DFS~~~~~~~~~~~~~~~~~\n";
					std::cout.precision(10);
					std::cout << "[" << mex::timeSystem() - startTime << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
					std::cout << "Tree size: " << GSIZE <<"\n";
					std::cout << "root information:" << std::endl;
					printNode(*rt);
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				}
				indicator = DFS_Beta(startTime, outFrequency);
				if (verbose > 2) {
					std::cout << "~~~~~~~~~~~~~~~~~after DFS~~~~~~~~~~~~~~~~~\n";
					std::cout.precision(10);
					std::cout << "[" << mex::timeSystem() - startTime << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
					std::cout << "Tree size: " << GSIZE <<"\n";
					std::cout << "root information:" << std::endl;
					printNode(*rt);
					std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
				}
			}
			else {
				std::runtime_error("Unsupported search mode!");
			}
			// break
			if ( indicator != 0 ) {
//				std::cout <<"all nodes have been removed!" <<std::endl;
				break;
			}

		}
		if (UB - LB < boundGap) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startTime << "]: "
					<< std::setw(10) << LB << " < ln Z < " << UB << ", UB-LB = "
					<< UB - LB << " < " << boundGap << "\n";
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			boundGap /= 10.0;
		}
		if ( UB - LB < solvedThresh ) {
			std::cout << "UB - LB < "<< solvedThresh << ", quit!\n";
			break;
		}

		if ( mex::timeSystem() - startTime > timeLimit ) {
			std::cout << "Reach search time limit (sec): "<< timeLimit  << ", quit!\n";
			break;
		}
	}
	LB = rt->lb;
	UB = rt->ub;
	std::cout.precision(10);
	std::cout<<"["<<mex::timeSystem()-startTime<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
	std::cout << "Tree size: " << GSIZE <<"\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
}

//EOF

/*
 * searchsample.cpp
 *
 *  Created on: Oct 10, 2016
 *      Author: qlou
 */

#include "searchsample.h"

// namespace mex {

searchsample::searchsample(mex::wmbe& wmbeub, double ub, int ns, double delta, int vb, double tbgt, double mbgt, int nv, bool isSolved):
wmbUB(wmbeub), ascList(wmbeub.nvar()), UB(ub), nSample(ns), DELTA(delta), verbose(vb), timeBudget(tbgt), memoryBudget(mbgt),
exactHeur(nv, true), tuple(nv,-1), initUB(ub), _isExactHeur(isSolved) {
	startProcess = mex::timeSystem();
	order = wmbeub.getOrder();
	priority = wmbeub.getPriority();
	nvar = wmbeub.nvar();
	GSIZE = 0;
	LB = -std::numeric_limits<double>::infinity();
//	iBound = wmbeub.getIBound();
	NSP = 0;
	buildPseudoTree();
	setRootConfig();
	setExactHeur();
}
searchsample::~searchsample() {
	// TODO Auto-generated destructor stub
	std::cout<<"==== Done Search+Sampling ====\n";
}
void searchsample::buildPseudoTree() {
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
	std::cout<<"==== Pseudo tree built ====\n";

	if (verbose > 2) {
		std::cout << "no. of nodes without parents: " << rts.size() << "/"
				<< nvar << "\n";
		for (size_t i = 0; i < ascList.size(); ++i) {
			std::cout << "variable " << i << " has no. of children: " << (ascList[i]).size() << "\n";
		}
	}
}

void searchsample::setRootConfig() {
	// initialize root
	// this is a virtual root node, nodes without parents are children of this virtual root node
	node vtRoot; // virtual root
	vtRoot.type = AND_NODE; // virtual AND_NODE node
	vtRoot.id = -1; // flag
//	vtRoot.lb = LB;
	vtRoot.ub = UB;
	//
	exploredTree.insert(exploredTree.begin(), vtRoot);
	root = exploredTree.begin();
	GSIZE = 1;
	std::cout<<"==== Root configuration set ====\n";
}
void searchsample::setExactHeur() {
/*
 * set exactHeur via breadth-first search from leaf nodes of the pseudo tree
 * exactHeur tells whether the heuristics below that node is exact or not
 * note that exactHeur is initialized to be true
 */
	std::vector<size_t> nVis(nvar,0); // no. of children visited so far
	auto parents = wmbUB.getPseudotree();
	std::queue<mex::graphModel::vindex> queue; // variable id
	for (int i=0; i<nvar; ++i) {
		// add leaf nodes
		if (ascList[i].empty()) queue.push(i);
	}
	while(!queue.empty()) {
		auto top = queue.front(); queue.pop();
		auto par =  parents[top];
		if ( par >= nvar  ||  par < 0 )  {
			// no parent in pseudo tree
			continue;
		}
		++nVis[par];
		if ( wmbUB.getNumberOfMiniBuckets( wmbUB.var(top) ) > 1  || !exactHeur[top] ) {
			// multiple mini-buckets or already non-exact at the lower level
			exactHeur[par] = false;
		}
		// push a node into queue after visiting its last child
		// this prevents adding the same node multiple times;
		if (nVis[par] == ascList[par].size()) {
			queue.push(par);
		}
	}
	if (verbose > 2) {
		std::cout<<"~~~~~~~~~~~~~~~Printing exactHeur~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		for (int i=0; i<nvar; ++i) {
			std::cout << "[" << i << "," << exactHeur[i] << "], ";
			assert(nVis[i] == ascList[i].size());
		}
		std::cout<<"\n~~~~~~~~~~~~~~~Done printing~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	}
}
// helper functions
double searchsample::logsumexp(std::vector<double>& vec) {
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
// c++ 11 feature
double searchsample::logsumexp(std::initializer_list<double> vec) {
	double mx = *(std::max_element(vec.begin(), vec.end()));
	double r = 0;
	if (mx <= std::numeric_limits<double>::lowest()) {
		return mx;
	}
	for (auto val : vec) {
		r += std::exp(val - mx);
	}
	return (std::log(r) + mx);
}
void searchsample::printNode(const node& n) {
	return;
}
bool searchsample::checkSolved(const node& n) {
	// works for both AND_NODE and OR_NODE nodes
	// caveat: may only apply to frontier nodes!
	bool solved = false;
	if (n.id > -1) {
		// if not root
//		solved = (ascList[n.id]).empty() || ( n.ub - n.lb <= solvedThresh );
		solved = exactHeur[n.id]; // for leaf node in pseudo tree, this is always true
	}
//	else {
//		solved = (n.ub - n.lb <= solvedThresh) ;
//	}
	return solved;
}
tree<searchsample::node>::iterator searchsample::pruneSolved(tree<node>::iterator ptr) {
/*
 * Must run from the newly expanded node
 * Caveat: a solved node will be deleted only when all its siblings are solved
 * this is crucial (to Rao-blackwellisation) here
 * return the highest UNsolved node
 */
//	auto rt = exploredTree.begin();
	if ( (ptr==NULL) || (!ptr->solved) ) return ptr;
	auto highestSolved = ptr;
	ptr = exploredTree.parent(ptr); // since ptr is solved, start from its parent

	auto rtPar = exploredTree.parent(root);
	while ( ptr != rtPar ) {
		bool solved = true;
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// a node is solved only when all its children are solved
			if( !sib->solved ) {
				solved = false;
				break;
			}
		}
		// caveat: unsolved children cannot be removed alone even in SMA*
		ptr->solved = solved;
		// a node is solved only when all its children are solved
		if (!ptr->solved) {
			break;
		}
		//
		highestSolved = ptr;
		ptr = exploredTree.parent(ptr);
	}
	assert(highestSolved->solved);
	// erase the whole subtree rooted at this node (NOT including itself)
	// if not the root
/*	ptr = exploredTree.parent(highestSolved);
	if ((*ptr).type == AND_NODE) {
		(*ptr).exact_val += (*highestSolved).ub; // at this moment, ub=lb
	} else {
		double maxval = std::max((*ptr).exact_val, (*highestSolved).ub);
		(*ptr).exact_val = maxval
				+ log(
						exp((*ptr).exact_val - maxval)
								+ exp((*highestSolved).ub - maxval));
	}
	// the parent points to some child for safety
	if ((*ptr).toChild == highestSolved) {
		(*ptr).toChild = NULL;
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			if (sib != highestSolved) {
				// no need to pick the best child at this moment, we'll do it in backward pass
				(*ptr).toChild = sib;
				break;
			}
		}
	}*/
	// delete child one by one
	// no need to back up those values in exact_val since we delete all solved children at once.
/*	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		if ( highestSolved->type == AND_NODE ) {
			highestSolved->exact_val += sib->ub;
		} else {
			double maxval = std::max(highestSolved->exact_val, sib->ub);
			highestSolved->exact_val = maxval + log(exp(highestSolved->exact_val - maxval) + exp(sib->ub - maxval));
		}
	}*/
	// finally delete
	int GSIZE_pre = GSIZE;
	// size gives you the size of the subtree rooted at given node ( including that node )
//	GSIZE -= exploredTree.size(highestSolved);
//	exploredTree.erase(highestSolved);
	GSIZE -= exploredTree.size(highestSolved); // size of the subtree rooted at the given node (including that node)
	++GSIZE; // do not count highestSolved
	exploredTree.erase_children(highestSolved);

	if (verbose > 2) {
		std::cout << "pruning done: " << GSIZE_pre - GSIZE << " nodes deleted!\n" ;
	}
//	return ptr;
	// return parent if not root
	if ( highestSolved == root ) return root;
	return exploredTree.parent(highestSolved);
}
void searchsample::updateBounds(tree<node>::iterator ptr) {
	// re-calculate bounds from existing children and exact_val
	if (verbose > 2) {
		std::cout << "now update bounds: " << std::endl;
	}
	assert( ptr != NULL );
//	auto rt = exploredTree.begin();
//	double cur_lb,
	ptr = exploredTree.parent(ptr);
	auto top = exploredTree.parent(root);
	double cur_ub = 0.0;

	while( ptr != top ) {
		// update bounds by re-calculation.
		// more time-consuming than simple update, but more numerically stable
		if (ptr->type == AND_NODE) {
//			cur_lb = (*ptr).exact_val;
			cur_ub = ptr->exact_val;
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//				cur_lb += (*sib).lb;
				cur_ub += sib->ub;
			}
		}
		else {
			// OR_NODE node
			std::vector<double> ub_vec {(*ptr).exact_val};
//			lb_vec {(*ptr).exact_val};
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//				lb_vec.push_back( (*sib).lb );
				ub_vec.push_back( sib->ub );
			}
//			cur_lb = logsumexp(lb_vec);
			cur_ub = logsumexp(ub_vec);
		}
//		if ( cur_lb > ptr->lb || cur_ub < ptr->ub ) {
		if ( cur_ub < ptr->ub ) {
//			ptr->lb = cur_lb;
			ptr->ub = cur_ub;
		} else {
			// no need to move upward
			break;
		}
		ptr = exploredTree.parent(ptr);
	}
	// update bounds for the root
//	LB = rt->lb;
	UB = root->ub;
}
void searchsample::findChild(tree<node>::iterator ptr) {
/*
 * find the best child for an UNsolved node
 * and update its info accordingly
 * a solved child will not be considered
 * we assume no removeWorst
 */
	assert(ptr!=NULL); // should not be NULL
	if ( ptr->solved || (exploredTree.begin(ptr) == exploredTree.end(ptr) ) ) return;
//  note that prior here is not the 'actual' priority value, just a relative quantity for comparison.
	double bestPrior = -std::numeric_limits<double>::infinity();
	double prior = 0.0;
	tree<node>::iterator toChild = NULL;
	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// for upper priority, no need to exactly compute the priority to identify which child is the best
		if (sib->solved) continue; // do not consider a solved child
		if ( sib->type == AND_NODE ) {
//			if (sib->toChild != NULL) {
				prior = sib->down_ub + sib->exact_val;
//			} else {
//				prior = sib->ub;
//			}
		}
		else {
			// even work for frontier OR node
			// prt is AND ndoe
			prior =  sib->down_ub + ptr->ub - sib->ub; // since no removeWorst
		}
		if ((toChild == NULL) || (prior > bestPrior) ) {
			toChild = sib;
			bestPrior = prior;
		}
	}
	// if not solved
	assert(toChild != NULL);
	ptr->toChild = toChild;
	// not including exact_val
//		ptr->downPath_lb = toChild->downPath_lb + cur_lb - ptr->exact_val - toChild->lb;
	if (ptr->type == AND_NODE) {
		ptr->down_ub = toChild->down_ub + ptr->ub - ptr->exact_val - toChild->ub;
	} else {
		// if toChild is frontier node, the following equals toChild->ub
		ptr->down_ub = toChild->down_ub + toChild->exact_val;
	}
}
int searchsample::expandBest() {
	/*
	 * forward pass to expand the best node in OPEN
	 */
	if (verbose > 2) {
		std::cout << "expandBest: running..." << std::endl;
	}
//	mex::vector<uint32_t> tuple(nvar, -1); // for safety, initialized to -1
	std::fill(tuple.begin(), tuple.end(), -1); // re-use to accelerate, for safety, initialized to -1
//	auto ptr = exploredTree.begin(); // the root is AND node
	auto ptr = root;
			// forward pass, also update topPath_ub, topPath_lb
	while (ptr->toChild != NULL) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		if (verbose > 2) printNode(*ptr); // debug
		ptr = ptr->toChild;
		if (ptr->type == AND_NODE) {
			tuple[ptr->X] = ptr->val;
		}
	}
	if (verbose > 2) {
		std::cout << "Best node found!" << std::endl;
		printNode(*ptr);
	}

	if (ptr->solved) {
		std::cout << "The best frontier node has been solved, quit!" << std::endl;
		return 3;
	}

	assert(ptr->type == AND_NODE); // without removal, the frontier must be AND

//	if (ptr->type == AND_NODE) {
//		double cur_lb=0.0, cur_ub = 0.0;
	double cur_ub = 0.0;
	std::list<mex::graphModel::vindex> childrenList;
	if (ptr->id > -1) {
		childrenList = ascList[ptr->id];
	} else {
		// if it is the virtual root
		childrenList = rts;
	}
	// if the most promising node corresponds to a leaf node of the pseudo tree,
	// the whole search progress should be terminated, although it's possible that the bound gap remains non-zero.
	// it depends on what priority you use
	// no need to check, because this node should be marked solved
/*	if (childrenList.empty()) {
		std::cout << "The best node: variable " << (*ptr).id << " with value "
				<< (*ptr).val << " has no children, quit!\n";
		printNode(*ptr);
		if (verbose > 1) {
			while ((*ptr).id > -1) {
				printNode(*ptr);
				ptr = exploredTree.parent(ptr);
			}
			printNode(*ptr);
		}
		return 2;
	}*/
	// update bounds later
	bool isProp = false;
//		if (ptr->isReExpand) {
	// if expanded before, reset
	// also no need to propagate bounds

	// should not be necessary to re-calculate again
//	if (ptr->id > -1)
//		ptr->exact_val = wmbUB.heuristicTheta(ptr->X, tuple);
//			isProp = true;
//		}
//	cur_lb = cur_ub = ptr->exact_val;
	cur_ub = ptr->exact_val;
	// initialization
	ptr->solved = true;
	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
		node child;
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
//		mex::Factor dLB(child.X, 0.0);
		//
/*		int nMini = wmbUB.duplicateInBucket(child.X, tuple);
		isProp = isProp || (nMini > 1);
		if (verbose > 2) {
			if (nMini > 1) {
				std::cout
						<< "Bounds should be improved when expanding this node, no. of mini-buckets: "
						<< nMini << std::endl;
			}
		}*/
		if (verbose > 2) {
//				int nMini = wmbUB.duplicateInBucket(child.X, tuple);
			int nMini = wmbUB.getNumberOfMiniBuckets(child.X);
//				isProp = isProp || (nMini > 1);
			if (nMini > 1) {
				std::cout
						<< "Bounds should be improved after expanding this node, no. of mini-buckets: "
						<< nMini << std::endl;
			}
		}
		// initialization
		ptrChild->solved = true;
		for (size_t v = 0; v < (child.X).states(); ++v) {
			node kid;
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
//				kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
			kid.down_ub = wmbUB.heuristicIn(kid.X, tuple);
//				kid.lb = dLB[v] = kid.downPath_lb + kid.exact_val;
			kid.ub = dUB[v] = kid.down_ub + kid.exact_val;
			// whether solved
			kid.solved = checkSolved(kid);
			if (!kid.solved)
				ptrChild->solved = false;

			// delete a solved node only when all its siblings are solved, i.e., its parent is solved
			// it is serve as a book-keeping purpose, essential for two-stage sampling
			// append
			exploredTree.append_child(ptrChild, kid);
			++GSIZE;
		}
		if (!ptrChild->solved)
			ptr->solved = false;
//			ptrChild->lb = dLB.logsumexp();
		ptrChild->ub = dUB.logsumexp();
//			cur_lb += ptrChild->lb;
		cur_ub += ptrChild->ub;
		// if solved, remove all its children
		if (ptrChild->solved) {
			GSIZE -= exploredTree.number_of_children(ptrChild);
			exploredTree.erase_children(ptrChild);
		}
		// find best child
		findChild(ptrChild);
	}
	// more stable
	if ( ptr->ub > cur_ub + epsilon || ptr == root) {
		isProp = true;
	}
	// propagate bounds if necessary
	if (isProp) {
//		if (isProp && !ptr->isReExpand) {
//			ptr->lb = cur_lb;
		ptr->ub = cur_ub;
		updateBounds(ptr);
	}
	if (verbose > 2) {
		std::cout << "the best frontier node after expansion" << std::endl;
		printNode(*ptr);
	}
	// caveat: make sure bounds have been fully updated before pruning
	ptr = pruneSolved(ptr); // lowest UNsolved
	return backwardUpdate(ptr);
}
int searchsample::backwardUpdate(tree<node>::iterator ptr) {
/*
 * backward update after expanding the best
 */
	if (verbose > 2) {
		std::cout << "running backwardUpdate..." << std::endl;
	}

	if (ptr == NULL) {
		std::cout << "backwardUpdate: NULL occurs, problematic or in purpose!" << std::endl;
		return -1;
	}
	// re-calculate the best child since it may point to a solved node
	// check whether the whole tree has been solved or not
//	auto rt = exploredTree.begin();
	if ( root->solved ) {
		std::cout << "The root is solved!\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		return 1;
	}
	while (true) {
		findChild(ptr);
		// go upward
		if ( ptr == root ) {
			break;
		}
		ptr = exploredTree.parent(ptr);
	}
	return 0;
}
void searchsample::runSearch(const unsigned treeSizeLimit) {
	std::cout << "Run into search...\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	double startIter = mex::timeSystem();
	double outFrequency = 0.1; // for output
	bool outNow = false;

//	auto rt = exploredTree.begin();
	if (verbose > 0) {
//		LB = rt->lb;
		UB = root->ub;
		std::cout.precision(10);
		std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}

	while(true){
		if ( root->solved ) {
//			std::cout << "The root is solved, quit!"<<std::endl;
//			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			break;
		}
		if ( GSIZE >= treeSizeLimit ) {
			std::cout << "Reach tree size limit: " << treeSizeLimit << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			break;
		}
		if ( mex::timeSystem()-startProcess > timeBudget ) {
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			break;
		}
		//
		expandBest();
		//
		outNow = verbose>1;
		if ( mex::timeSystem() - startIter > outFrequency  ) {
			outNow = true;
			outFrequency += outFrequency; // exponential
		}

		if (outNow) {
//			LB = rt->lb;
			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}
	}
//	LB = rt->lb;
	UB = root->ub;
	std::cout.precision(10);
	std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
	std::cout << "Tree size: " << GSIZE <<"\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
}
int searchsample::runSearch(const unsigned nnd, const unsigned treeSizeLimit) {
/*
 * nnd: no. of nodes to expand
 *
 */
//	int base = GSIZE;
	double power = 1.0 + std::floor(std::log2(GSIZE+1)); // add

	if (verbose >1) {
		std::cout << "Run " << nnd << " rounds of best-first search...\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	unsigned cnt=0;
	while(cnt < nnd){
		++cnt;
		if ( root->solved ) {
			std::cout << "The root is solved, done!"<<std::endl;
			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//			break;
			return SOLVED;
		}
		if ( GSIZE >= treeSizeLimit ) {
			std::cout << "Reach tree size limit: " << treeSizeLimit << std::endl;
			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//			break;
			return MEMOUT;
		}
		if ( mex::timeSystem()-startProcess > timeBudget ) {
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//			break;
			return TIMEOUT;
		}
		if ( (GSIZE == int(pow(2.0, power))) && (verbose > 0) ) {
			power += 1.0;
			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}
		//
		expandBest();

/*		if (outNow) {
//			LB = rt->lb;
			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}*/
	}
//	LB = rt->lb;
/*	UB = root->ub;
	std::cout.precision(10);
	std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
	std::cout << "Tree size: " << GSIZE <<"\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;*/

//	LB = rt->lb;
	UB = root->ub;
	if (verbose > 1) {
		std::cout.precision(10);
		std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	return 0;
}
double searchsample::logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done) {
/*
 * compute logF(x) over a (possibly partial) configuration of a set of variables
 * If a variable is in the list, then all its ancestors on the pseudo tree must be in the list.
 * if the list is a complete list of all variables, use a fast approach.
 */
	assert(Done.size() < nvar+1); // debug
	if (Done.size() == nvar) return wmbUB.logP(config);
	//
	double logVal = 0.0;
	for (const auto& X : Done) {
		logVal += wmbUB.heuristicTheta(X, config);
	}
	return logVal;
}
tree<searchsample::node>::iterator searchsample::upperSampling(tree<node>::iterator ptr){
/*
 * sample for an OR node based on current upper bounds of its AND children;
 * for each AND child, probability to be picked is proportional to its current upper bound.
 */
	assert( (*ptr).type == OR_NODE );
	auto X = (*ptr).X;
	auto ns = X.states();
	// some children may be solved, and serve as a placeholder
	assert( ns==exploredTree.number_of_children(ptr) );

	mex::Factor fc(X, 0.0);

	int i = 0;
	for ( auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib ) {
		// !!!caveat: factor's sample() function is NOT in log!!!
		// do NOT assume factors are in log!!!
		fc[i] = exp( (*sib).ub );
		++i;
	}
//	fc /= fc.sum();  // we do not have to normalize it to do sampling, actually
	int val = fc.sample(); // NOT in log
	// find the AND child corresponding to the sampled value
	for ( auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib ) {
		if ( (*sib).val == val )
			return sib;
	}
	assert(false); // should not run into this line
	return NULL;
}
std::list<mex::Var> searchsample::getDescList(int id) {
/*
 * get a list of all descendants of a variable in the pseudo tree
 * order in the list matters, ancestor always prior to descendant
 */
	std::list<mex::Var> desc;
	mex::Var X;
	std::stack<mex::Var> stack;
	std::list<mex::graphModel::vindex> childrenList;
	if ( id  > -1 ) {
		childrenList = ascList[id];
	}
	else {
		// if it is the virtual root
		childrenList = rts;
	}
	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
		stack.push(wmbUB.var(*it));
	}
	while(!stack.empty()) {
		X = stack.top();
		stack.pop();
		desc.push_back(X);
		childrenList = ascList[X.label()]; //  label = id
		for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
			stack.push(wmbUB.var(*it));
		}
	}
	return desc;
}
void searchsample::addSample(double logFx, const std::vector<double>& logQxCond) {
/*
 * add one sample to the root
 */
	auto ptr = root;
	double wt = 0.0;
	double logQx = 0.0;

	++NSP;
	++ptr->nsp; // add count first
	// debug
	assert(NSP == ptr->nsp);

	// sum all conditional probabilities
	for(const auto& v : logQxCond) logQx += v;
	//
	wt = logFx - logQx;
	// for OR nodes, only add one count and done
//	if ( (*ptr).type == OR_NODE ) return;
	if ((*ptr).nsp <= 1) {
		(*ptr).logEx = wt;
		(*ptr).logEx2 = 2 * wt;
	} else {
		int nsp = (*ptr).nsp;
		double logEx = (*ptr).logEx;
		double logEx2 = (*ptr).logEx2;
		// to be numerically stable
		// Ex  = (Ex * (samp-1))/samp + dEx/samp;
		(*ptr).logEx = logsumexp( { log(nsp - 1) + logEx, wt }) - log(nsp);
		// Ex2 = (Ex2 * (samp-1))/samp + dEx*dEx/samp;
		(*ptr).logEx2 = logsumexp( { log(nsp - 1) + logEx2, 2 * wt })
				- log(nsp);
	}
	// update necessary quantities
	updateConvEstimate(wt, root->ub);
	updateNormalizedEstimate(wt, root->ub);
	// debug
	if (verbose > 1 && DEBUG) {
		std::cout << "NSP: " << NSP <<  ", Ex: " << ptr->logEx << ", Ex2: " << ptr->logEx2 << "\n";
		std::cout << "convEx: " << convEx << ", sumSqrUB: " << sumSqrUB << ", sumInvSqrUB: " << sumInvSqrUB << "\n";
		std::cout << "normalizedEx: " << normalizedEx << ", normalizedEx2: " << normalizedEx2 << "\n";
		std::cout << "sumInvUB: " << sumInvUB << ", hmUB: " << hmUB << std::endl;
	}
}
std::pair<double, double> searchsample::calcRootEBB(){
	// calculate EBB for the root, since we can use the refined upper bound for the root
	auto ptr = root;
	double lb = -std::numeric_limits<double>::infinity(), ub = std::numeric_limits<double>::infinity();
////	calcNodeBias
	if ( (*ptr).nsp <= 1 ) {
		// deterministic bound, Z <= Zhat + Z^{+}
//		(*ptr).logBias = (*ptr).ub;
		return std::pair<double,double>(lb,ub);
	}

	int nsp = (*ptr).nsp;
//	double confidence = (*ptr).delta;
	double Ex = exp( (*ptr).logEx - (*ptr).ub  ); // normalized
	double Ex2 = exp( (*ptr).logEx2 - 2*(*ptr).ub  ); // normalized

	double var = std::max(Ex2 - Ex*Ex, 0.0); // to be numerically stable
	var *= ( nsp/(nsp-1) ); // should be the unbiased sample variance

	double rng = std::sqrt(2*var*std::log(2.0/DELTA)/nsp) + 7*std::log(2.0/DELTA)/3.0/(nsp-1);
//	(*ptr).logBias = std::log(rng) + (*ptr).ub;
	double logBias = std::log(rng) + (*ptr).ub;
	//
	if ( (*ptr).nsp > 1 ) {
		ub = logsumexp( {logBias,  (*ptr).logEx} );
		if ( (*ptr).logEx > logBias) {
				lb = log( exp( (*ptr).logEx ) - exp( logBias ) );
		}
	}
	return std::pair<double,double>(lb,ub);
}
std::pair<double, double> searchsample::calcRootEBB(double upperbound){
	/*
	 * calculate EBB for root, using given upper bound
	 */
	auto ptr = root;
	assert(upperbound > ptr->ub - epsilon); // should be worse than current bound
	double lb = -std::numeric_limits<double>::infinity(), ub = std::numeric_limits<double>::infinity();
////	calcNodeBias
	int nsp = ptr->nsp;
	if ( nsp <= 1 ) {
		// deterministic bound, Z <= Zhat + Z^{+}
//		(*ptr).logBias = (*ptr).ub;
		return std::pair<double,double>(lb,ub);
	}
//	double confidence = (*ptr).delta;
	double Ex = exp( ptr->logEx - upperbound  ); // normalized
	double Ex2 = exp( ptr->logEx2 - 2*upperbound  ); // normalized

	double var = std::max(Ex2 - Ex*Ex, 0.0); // to be numerically stable
	var *= ( nsp/(nsp-1) ); // should be the unbiased sample variance

	double rng = std::sqrt(2*var*std::log(2.0/DELTA)/nsp) + 7*std::log(2.0/DELTA)/3.0/(nsp-1);
//	ptr->logBias = std::log(rng) + upperbound;
	double logBias = std::log(rng) + upperbound;
	//
	if ( nsp > 1 ) {
		ub = logsumexp( {logBias,  ptr->logEx} );
		if ( ptr->logEx > logBias) {
				lb = log( exp( ptr->logEx ) - exp( logBias ) );
		}
	}
	return std::pair<double,double>(lb,ub);
}
std::pair<double, double> searchsample::calcHoeffding() {
/*
 * Hoeffding's bounds on the convex combination of samples
 * convex weights are proportional to 1/U_i^2
 */
	double bias = ( log( -log(DELTA)/2.0 ) - sumInvSqrUB )/2.0; //
	double lb = -std::numeric_limits<double>::infinity();
//	double ub = log(exp(convEx) + rng);
	double ub = logsumexp( {convEx, bias} );
	if ( convEx > bias) {
//			lb = log( exp( convEx ) - exp( bias ) );
		lb = convEx + log(1 - exp( bias-convEx ) );
	}
	return std::pair<double,double>(lb,ub);
}
std::pair<double, double> searchsample::calcNormalizedEBB() {
/*
 * EBB based on the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
 *
 */
	double lb = -std::numeric_limits<double>::infinity(), ub = std::numeric_limits<double>::infinity();
	if ( NSP <= 1 ) {
		return std::pair<double,double>(lb,ub);
	}

	double zhat = hmUB + normalizedEx; // unbiased estimate of Z
	double var = std::max(exp(normalizedEx2) - exp(normalizedEx)*exp(normalizedEx), 0.0); // to be numerically stable
	var *= NSP/(NSP-1); // should be the unbiased sample variance
	double bias = log( sqrt(2.0*var*log(2.0/DELTA)/NSP) + 7.0*log(2.0/DELTA)/3.0/(NSP-1) ) + hmUB; // store in log
	ub = logsumexp( {bias, zhat} );
	if ( zhat > bias) {
		lb = zhat + log(1 - exp( bias-zhat ) );
	}
	return std::pair<double,double>(lb,ub);
}
void searchsample::updateConvEstimate(double wt, double ub) {
/*
 * update convEx, sumSqrUB, sumInvSqrUB
 * wt: importance weight, i.e., current sample (in log)
 * ub: corresponding global upper bound on the sample (in log)
 * Caveat: NSP should've already been updated
 */
	if (NSP <= 1) {
		// first sample
		// actually, we do not treat this as a special case due to the initialization of those quantities.
		// however, we do this to make code more readable and robust.
		convEx = wt;
		sumSqrUB = 2*ub;
		sumInvSqrUB = -sumSqrUB;
		return;
	}
	double ex = convEx + sumInvSqrUB; // \sum z_i/u_i^2
	ex = logsumexp( {ex, wt-2*ub} );
	sumSqrUB = logsumexp( {sumSqrUB, 2*ub} );
	sumInvSqrUB = logsumexp( {sumInvSqrUB, -2*ub} );
	convEx = ex - sumInvSqrUB;
}
void searchsample::updateNormalizedEstimate(double wt, double ub) {
/*
 * update the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
 * empirical bernstein bound is applicable to this estimate
 * wt: importance weight, i.e., current sample (in log)
 * ub: corresponding global upper bound on the sample (in log)
 * Caveat: NSP should've already been updated
 */
	if (NSP <= 1) {
		// first sample
		normalizedEx = wt - ub;
		normalizedEx2 = 2*(wt-ub);
		sumInvUB = -ub;
		hmUB = ub;
		return;
	}
	normalizedEx = logsumexp( {normalizedEx+log(NSP-1), wt-ub} ) - log(NSP);
	normalizedEx2 = logsumexp( {normalizedEx2+log(NSP-1), 2*(wt-ub)} ) - log(NSP);
	sumInvUB = logsumexp( {sumInvUB, -ub} );
	hmUB = log(NSP) - sumInvUB;
}
void searchsample::twoStageSampling() {
/*
 * two-stage sampling with a fixed tree
 */
	mex::vector<uint32_t> config(nvar, -1); // -1 for safety
	// logQxCond[i] = log( q(x_i | par(x_i)) ), par(x_i) is the set of all ancestors on the pseudo tree
	std::vector<double> logQxCond(nvar, 0.0);
	std::list<mex::Var> Done; // list of all sampled variables
	// sample the path from root to leaf
//	auto rt = exploredTree.begin();
	auto ptr = root, ptrPar = root;
	//values for solved variables, summation part of Rao-Blackwellisation
	double solved_val = 0.0;
	// stage 1 and 2
	std::stack<tree<node>::iterator> stack;
	stack.push(ptr);
	while (!stack.empty()) {
		ptr = stack.top();
		stack.pop();
		if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
			// frontier node
			// any solved node should be a frontier node due to pruning
			// no sampling for a solved frontier
			if (ptr->solved) {
				if (ptr->type == AND_NODE) {
					// For solved AND node, the corresponding variable is IN config.
					solved_val += ptr->ub - ptr->exact_val; // exact_val will be computed logF(config, Done);
				} else {
					// For solved OR node, the corresponding variable is NOT in config.
					solved_val += ptr->ub;
				}
			} else {
				// an unsolved frontier node should be AND node due to the way we expand the best
				assert(ptr->type == AND_NODE);
				// sample from the frontier node according to the mixture proposal
				int ID = -1; // default, root id
				if (ptr->id > -1)
					ID = (ptr->X).label();
				// test, should be the same
				assert(ptr->id == ID);
				// variables for sampling
				auto UnDone = getDescList(ID);
				// debug
				if (DEBUG) {
					std::cout << "ID: " << ID << "\n";
					std::cout << "UnDone: ";
					for (auto X : UnDone) {
						std::cout << X.label() << ", ";
					}
					std::cout
							<< "\n++++++++++++++++++++++++++++++++++++++++++++++++++++"
							<< std::endl;
				}

				wmbUB.conditionalMixtureSampling(ID, config, logQxCond, UnDone); // config, logQxVec, logFx will be updated
//				addSampleToNode(exploredTree.begin(), logFx - logQxVec[nvar]);
						// add a set of sampled variables
				Done.insert(Done.end(), UnDone.begin(), UnDone.end());
			}
		} else {
			// internal node
			assert(!ptr->solved); // debug, internal nodes should not be solved due to pruning
			if (ptr->type == OR_NODE) {
				// not solved,
				// an AND child returned
				ptrPar = ptr;
				ptr = upperSampling(ptr);
				stack.push(ptr);
				//
				config[ptr->X] = ptr->val;
//			loc = priority[ (ptr->X).label() ];
//			logQxVec[loc+1] = ptr->ub - ptrPar->ub;
				logQxCond[ptr->id] = ptr->ub - ptrPar->ub; // ub of an OR parent is equal to sum of all AND children's ub
				// add a sampled variable
				Done.push_back(ptr->X);
			} else {
				// push all OR children to stack
				for (auto sib = exploredTree.begin(ptr);
						sib != exploredTree.end(ptr); ++sib) {
					stack.push(sib);
				}
			}
		}
	}
	// add the sample to the root node
	// TODO add the sample to those internal nodes as well?
	double logFx = logF(config, Done);

	// debug
	if (DEBUG) {
		int nass = 0; // no. of assigned X
		for (auto v : config) {
			if (v >= 0 && v < nvar ) ++nass;
		}
		std::cout << "no. of assigned X in config: " << nass
				<< "\nDone.size: " << Done.size()
				<< "\nlogFx: " << logFx
				<<"\nsolved_val: " << solved_val << std::endl;
		std::cout<< "config (logQxCond): ";
		for (int i=0; i<nvar; ++i) {
			std::cout << config[i] << " (" << logQxCond[i] << "), ";
		}
		std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}

	//
	logFx += solved_val; // add solved
	addSample(logFx, logQxCond);
}
void searchsample::doSampling() {
//	auto rt = exploredTree.begin();
	bool outNow = false;
	double outFrequency = 1;
	std::pair<double, double> prob;

	while(NSP < nSample) {
		twoStageSampling();
		outNow = verbose > 1;
		if (NSP >= outFrequency) {
			outNow = true;
			outFrequency += outFrequency;
		}
		if (outNow) {
			prob = calcRootEBB();
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
//			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		}

		if ( mex::timeSystem()-startProcess > timeBudget ) {
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			break;
		}
	}
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	if (NSP >= nSample) std::cout << "Reach sample size limit: " << nSample << std::endl;
	prob = calcRootEBB();
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
//	std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
	std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
	std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
}
int searchsample::doSampling(int nsp,  bool isInitUB) {
/*
 * draw nsp samples based on current search tree
 */
//	auto rt = exploredTree.begin();
//	bool outNow = false;
//	double outFrequency = 1;
	std::pair<double, double> prob;
	int cnt = 0;
//	double power = 1.0 + std::floor(std::log2(NSP+1)); // add 1 to NSP to avoid log(0)

	while(cnt < nsp) {
		++cnt;
		twoStageSampling();
//		outNow = verbose > 1;
//		if (NSP >= outFrequency) {
//			outNow = true;
//			outFrequency += outFrequency;
//		}
//		if ( (NSP == long(pow(2.0, power)))  || (verbose > 1) ) {
		if ( (NSP == _outFrequency)  || (verbose > 1) ) {
			_outFrequency += _outFrequency;
//			power += 1.0;
//			prob = calcRootEBB();
			if ( isInitUB ) prob = calcRootEBB(initUB);
			else prob = calcRootEBB();
			std::cout.precision(10);
//			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
			std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
//			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
			/***** more ways to aggregate estimates and bounds *****/
			prob = calcHoeffding();
			std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
			prob = calcNormalizedEBB();
			std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
			/*********************************/
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		}

		if ( mex::timeSystem()-startProcess > timeBudget ) {
			if ( isInitUB ) prob = calcRootEBB(initUB);
			else prob = calcRootEBB();
			std::cout.precision(10);
//			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
			std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
//			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
			/***** more ways to aggregate estimates and bounds *****/
			prob = calcHoeffding();
			std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
			prob = calcNormalizedEBB();
			std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
			/*********************************/
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			return TIMEOUT;
		}
	}
//	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//	if (NSP >= nSample) std::cout << "Reach sample size limit: " << nSample << std::endl;
	if (verbose > 1) {
//		prob = calcRootEBB();
//		prob = calcRootEBB(initUB);
		if ( isInitUB ) prob = calcRootEBB(initUB);
		else prob = calcRootEBB();
		std::cout.precision(10);
//		std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
		std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
//		std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
		std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
		std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
		/***** more ways to aggregate estimates and bounds *****/
		prob = calcHoeffding();
		std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
		prob = calcNormalizedEBB();
		std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
		/*********************************/
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	return 0;
}
void searchsample::start() {
//	DEBUG = true;
	// basic setup
//	buildPseudoTree();
//	setRootConfig();
//
//	auto rt = exploredTree.begin();
	// in byte
	int memUnit = sizeof(*root);
	// MemLimit is in megabyte
	auto treeSizeLimit = unsigned(1024*1024*memoryBudget/memUnit);
	std::cout << "Search time limit (sec): " << timeBudget << std::endl;
	std::cout << "Search memory limit (MB): " << memoryBudget << std::endl;
	std::cout << "Tree size limit: " << treeSizeLimit << std::endl;
	std::cout << "Sample size limit: " << nSample << std::endl;

	// run search first
	runSearch(treeSizeLimit);
	//
	if ( root->solved ) {
		std::cout << "Root solved, no sampling" << std::endl;
		return;
	}
	// do sampling
	doSampling();
}
void searchsample::start(unsigned treeSizeLimit) {
//	DEBUG = true;
	// basic setup

//	auto rt = exploredTree.begin();
	// in byte
//	int memUnit = sizeof(*rt);
//	// MemLimit is in megabyte
//	auto treeSizeLimit = unsigned(1024*1024*memoryBudget/memUnit);
	std::cout << "Search time limit (sec): " << timeBudget << std::endl;
	std::cout << "Search memory limit (MB): " << memoryBudget << std::endl;
	std::cout << "Tree size limit: " << treeSizeLimit << std::endl;
	std::cout << "Sample size limit: " << nSample << std::endl;

	// run search first
	runSearch(treeSizeLimit);
	//
	if ( root->solved ) {
		std::cout << "Root solved, no sampling" << std::endl;
		return;
	}
	// do sampling
	doSampling();
}
void searchsample::start(const int nsp, const int nnd) {
/*
 * the strategy is not quite clear yet
 * we may sample after each expansion or several expansions based on the complexity of drawing one sample and the complexity of one expansion?
 * how many samples should we draw in each sampling process
 * nsp: no. of samples to draw in each iteration
 * nnd: no. of nodes to expand in each iteration
 */
	if (_isExactHeur) {
		std::cout << "==== Problem solved during heuristic construction ====" << std::endl;
		return;
	}
	// in byte
	int memUnit = sizeof(*root);
	// MemLimit is in megabyte
	auto treeSizeLimit = unsigned(1024*1024*memoryBudget/memUnit);
	std::cout << "Time limit (sec): " << timeBudget << std::endl;
	std::cout << "Memory limit (MB): " << memoryBudget << std::endl;
	std::cout << "Tree size limit: " << treeSizeLimit << std::endl;
	std::cout << "Sample size limit: " << nSample << std::endl;
	std::cout << "nsp: " << nsp << ", nnd " << nnd <<  std::endl;
//	assert(nsp>0 || nnd>0);
//	const unsigned nrnd = 1; // tentatively
	bool isMemOut = false;
	int status = 0;
	while (true) {
		// -1: memory out
		// -2: time out
		// run search first
		if (!isMemOut && nnd>0) {
			status = runSearch(nnd, treeSizeLimit);
			if (status == MEMOUT) {
				// memory out
				isMemOut = true;
				std::cout << "no search anymore, only sampling" << std::endl;
			}
			if (status == TIMEOUT) break;
			if (status == SOLVED) break;
			// do sampling
			if (nsp>0) status = doSampling(nsp);
		} else {
			// if memory out, keep sampling until time out
			// nsp=0 implies we only do sampling after search, thus we can use the latest UB for EBB
			status = doSampling(std::numeric_limits<int>::max(), nsp > 0);
			if (status == TIMEOUT) break;
		}
	}
}
// EOF

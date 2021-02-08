/*
 * randheur.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: qlou
 */

#include "randheur.h"

//namespace std {
//randheur(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt, int nv, bool isSolved, int ksp);
randheur::randheur(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt, int nv,
		bool isSolved, int ksp, int depth_lookahead):
wmbUB(wmbeub), nvar(nv), ascList(nv), UB(ub), DELTA(delta), verbose(vb), timeBudget(tbgt), memoryBudget(mbgt),
exactHeur(nv, true), tuple(nv,-1), initUB(ub), _isExactHeur(isSolved), _ksp(ksp),
_depth_lookahead(depth_lookahead), _depthVec(nv, -1), _adaMaxDepthVec(nv, -1) {
	startProcess = mex::timeSystem();
	order = wmbeub.getOrder();
	priority = wmbeub.getPriority();
//	nvar = wmbeub.nvar();
	GSIZE = 0;
	LB = -std::numeric_limits<double>::infinity();
//	iBound = wmbeub.getIBound();
	NSP = 0;
	//
	buildPseudoTree();
	setRootConfig();
	setExactHeur();
	setDepth();
	setLookaheadDepth();
}
randheur::~randheur() {
	// TODO Auto-generated destructor stub
	std::cout<<"==== Done stochastic lookahead ===="<<std::endl;
}
void randheur::buildPseudoTree() {
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
	std::cout<<"==== Pseudo tree built ===="<<std::endl;;


	if (verbose > 2 || DEBUG) {
		if (DEBUG) {
			std::cout << "Elimination order:\n";
			for (const auto o : order) {
				std::cout << o << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "no. of nodes without parents: " << rts.size() << "/"
				<< nvar << "\n";
		for (size_t i = 0; i < ascList.size(); ++i) {
			std::cout << "variable " << i << " has no. of children: " << (ascList[i]).size() << "\n";
			if (DEBUG) {
				std::cout << "those children are: ";
				for (const auto j : ascList[i]) {
					std::cout << j << " ";
				}
				std::cout << std::endl;
			}
		}
	}
}
/*void randheur::setPseudoTree() {

 * create a pseudo tree of variables using a tree data structure
 * run after "buildPseudoTree"

	mex::Var dummy; // dummy root
	std::stack<tree<mex::Var>::iterator> stack;
	_pseudotree.insert(_pseudotree.begin(), dummy);
	for (auto it = rts.begin(); it != rts.end(); ++it) {
		auto X = wmbUB.var(*it);
		auto ptr = _pseudotree.append_child(_pseudotree.begin(), X);
		stack.push(ptr);
	}
	while (!stack.empty()) {
		auto par = stack.top();
		stack.pop();
		auto childrenList = ascList[par->label()]; //  label = id
		for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
			auto ptr = _pseudotree.append_child(par, wmbUB.var(*it));
			stack.push(ptr);
		}
	}
	// debug
	assert(_pseudotree.size() == nvar+1);
}*/
void randheur::setRootConfig() {
/*
 * initialize the search tree
 */
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
	std::cout<<"==== Root configuration set ===="<<std::endl;
}
void randheur::setExactHeur() {
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
void randheur::setDepth() {
/*
 * set depth vector via depth-first search
 * 0-based: a variable without a parent is of depth 0
 */
	std::stack<mex::graphModel::vindex> stack;
	std::list<mex::graphModel::vindex> childrenList;
	int height = 0;
	for (auto id : rts) {
		stack.push(id);
		_depthVec[id] = 0; // 0-based
	}
	while(!stack.empty()) {
		auto id = stack.top();
		stack.pop();
		childrenList = ascList[id]; //  label = id
		int depth = _depthVec[id];
		for (auto it : childrenList) {
			stack.push(it);
			_depthVec[it] = depth + 1;
		}
		if (depth > height) height = depth;
	}
	std::cout<<"==== Depth of each variable set ===="<<std::endl;
	std::cout << "Pseudo tree height: " << height + 1 << std::endl;

	if (DEBUG) {
		assert(_depthVec.size() == nvar);
		std::cout << "Depth of each variable:\n";
		for (int i=0; i < _depthVec.size(); ++i) {
			std::cout << i << " (" << _depthVec[i] << "), ";
		}
		std::cout << std::endl;
	}
}
void randheur::setLookaheadDepth() {
/*
 * set max lookahead depth for each variable
 * works for adaptive depth as well
 * strategy for being adaptive: max depth is set to the depth of first descendant that lives in multiple buckets
 * Caveat: run after "setDepth"
 */
	// fixed lookahead depth
	if (_depth_lookahead > 0) {
		for (int i=0; i<nvar; ++i) {
			_adaMaxDepthVec[i] = _depthVec[i] + _depth_lookahead;
		}

		if (DEBUG) {
			assert(_adaMaxDepthVec.size() == nvar);
			std::cout << "max depth of lookahead for each variable:\n";
			for (int i=0; i < _adaMaxDepthVec.size(); ++i) {
				std::cout << i << " (" << _adaMaxDepthVec[i] << "), ";
			}
			std::cout << "\nno. of mini-buckets for each variable:\n";
			for (int i=0; i < nvar; ++i) {
				std::cout << i << " (" << wmbUB.getNumberOfMiniBuckets(wmbUB.var(i)) << "), ";
			}
			std::cout << std::endl;
		}


		return;
	}
	// adaptive lookahead depth
	// current strategy: max depth is set to the depth of first descendant that lives in multiple buckets
	// traverse via breadth-first search
	// Caveat: _adaMaxDepthVec initialized to -1 is important
	auto parents = wmbUB.getPseudotree();
	std::queue<mex::graphModel::vindex> queue;
	for (auto id : rts) {
		queue.push(id);
	}
	while(!queue.empty()) {
		auto id = queue.front(); queue.pop();
		auto X = wmbUB.var(id);
		int depth = _depthVec[id];
		int nMini = wmbUB.getNumberOfMiniBuckets(X);
		// update upward until reach the first ancestor whose "max depth" has been set.
		if (nMini > 1) {
			auto par = parents[id]; // parent
			int parDepth = _adaMaxDepthVec[par]; // -1 for initialization
			while (parDepth < 0) {
				_adaMaxDepthVec[par] = depth;
				par = parents[par];
				if (par < 0 || par >= nvar) break;
				parDepth = _adaMaxDepthVec[par];
			}
		}
		for (auto j : ascList[id]) queue.push(j);
	}
	std::cout<<"==== Lookahead depth set ===="<<std::endl;

	if (DEBUG) {
		assert(_adaMaxDepthVec.size() == nvar);
		std::cout << "max depth of lookahead for each variable:\n";
		for (int i=0; i < _adaMaxDepthVec.size(); ++i) {
			std::cout << i << " (" << _adaMaxDepthVec[i] << "), ";
		}
		std::cout << "\nno. of mini-buckets for each variable:\n";
		for (int i=0; i < nvar; ++i) {
			std::cout << i << " (" << wmbUB.getNumberOfMiniBuckets(wmbUB.var(i)) << "), ";
		}
		std::cout << std::endl;
	}
}
// helper functions
double randheur::logsumexp(std::vector<double>& vec) {
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
double randheur::logsumexp(std::initializer_list<double> vec) {
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
void randheur::printNode(const node& n) {
	std::string str;
	if (n.type == AND_NODE){
		str = "AND node";
	} else{
		str = "OR node";
	}
	std::string sol_str;
	if (n.solved) {
		sol_str = "Solved ";
	} else {
		sol_str = "UnSolved ";
	}
	std::string node_str;
	if (n.toChild == NULL) {
		node_str = "frontier ";
	} else {
		node_str = "internal ";
	}

	std::cout << "~~~~~~~~~~~~~~\n" ;
	std::cout << node_str << sol_str <<  str << ", id " << n.id << ", val " << n.val << "\n";
	std::cout << "ub " << n.ub << ", exact_val " << n.exact_val << ",  est " << n.est << "\n";
	std::cout << "down_ub " << n.down_ub << " > down_est " << n.down_est << "\n";
	std::cout << "nsp: " << n.nsp <<  ", Ex: " << n.logEx << ", Ex2: " << n.logEx2 << std::endl;
//	std::cout << "ub (" << n.ub << ") > est (" << n.est << ")" << std::endl;
	return;
}
bool randheur::checkSolved(const node& n) {
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
tree<randheur::node>::iterator randheur::pruneSolved(tree<node>::iterator ptr) {
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
void randheur::updateBounds(tree<node>::iterator ptr) {
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
void randheur::updateEsts(tree<node>::iterator ptr){
/*
 * update bounds and estimates from the current node all the way back to the root
 * this is a substitute of updateBounds
 * Caveat: ways to update estimates
 * 1. use the new estimate only (current method)
 * 2. combine with the old estimate (how?)
 */
	if (verbose > 2) {
		std::cout << "now update bounds and estimates: " << std::endl;
	}
	ptr = exploredTree.parent(ptr);
	auto top = exploredTree.parent(root);
	double cur_ub = 0.0, cur_est = 0.0;

	while( ptr != top ) {
		// update bounds by re-calculation.
		// more time-consuming than simple update, but more numerically stable
		if (ptr->type == AND_NODE) {
//			cur_lb = (*ptr).exact_val;
			cur_ub = ptr->exact_val;
			cur_est = ptr->exact_val;
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//				cur_lb += (*sib).lb;
				cur_ub += sib->ub;
				cur_est += sib->est;
			}
		}
		else {
			// OR_NODE node
			std::vector<double> ub_vec {(*ptr).exact_val}, est_vec {(*ptr).exact_val};
//			lb_vec {(*ptr).exact_val};
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
//				lb_vec.push_back( (*sib).lb );
				ub_vec.push_back( sib->ub );
				est_vec.push_back( sib->est );
			}
//			cur_lb = logsumexp(lb_vec);
			cur_ub = logsumexp(ub_vec);
			cur_est = logsumexp(est_vec);
		}
/*//		if ( cur_lb > ptr->lb || cur_ub < ptr->ub ) {
		if ( cur_ub < ptr->ub ) {
//			ptr->lb = cur_lb;
			ptr->ub = cur_ub;
		} else {
			// no need to move upward
			break;
		}*/

		ptr->ub = cur_ub;
		ptr->est = cur_est;

		ptr = exploredTree.parent(ptr);
	}
	// update bounds for the root
//	LB = rt->lb;
	UB = root->ub;

}
void randheur::setEst(tree<node>::iterator ptr, mex::vector<uint32_t>& config){
/*
 * generate a fixed number of samples for a newly generated node
 * set estimate to the node
 */
	assert(!ptr->solved); // no need for solved nodes
	if (ptr->type == AND_NODE) {
		// sample from the frontier node according to the mixture proposal
		int ID = -1; // default, root id
		if (ptr->id > -1) {
			ID = ptr->id;
			config[ID] = ptr->val;
		}
		// variables for sampling
		auto UnDone = getDescList(ID);
		// logQxCond[i] = log( q(x_i | par(x_i)) ), par(x_i) is the set of all ancestors on the pseudo tree
		std::vector<double> logQxCond(nvar, 0.0);
		std::vector<double> spls(_nsp);
		for (int i=0; i<_nsp; ++i) {
			wmbUB.conditionalMixtureSampling(ID, config, logQxCond, UnDone); // config, logQxVec, logFx will be updated
			double qx = 0.0;
			for (const auto& v : logQxCond) qx += v;
			auto fx = logF(config, UnDone);
			spls[i] = fx - qx;
		}
		ptr->down_est = logsumexp(spls) - log(_nsp);
		ptr->est = ptr->down_est + ptr->exact_val;
	} else {
		int ncld = exploredTree.number_of_children(ptr);
		int i = 0;
		std::vector<double> cldEsts(ncld);
		for (auto sib=exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			cldEsts[i] = sib->est;
			++i;
		}
		ptr->est = logsumexp(cldEsts);
	}
}


/*void randheur::findChild(tree<node>::iterator ptr) {

 * find the best child for an UNsolved node
 * and update its info accordingly
 * a solved child will not be considered
 * we assume no removeWorst

	assert(ptr!=NULL); // should not be NULL
	if ( ptr->solved || (exploredTree.begin(ptr) == exploredTree.end(ptr) ) ) return;
//  note that prior here is not the 'actual' priority value, just a relative quantity for comparison.
	double bestPrior = -std::numeric_limits<double>::infinity();
	double prior = 0.0;
	tree<node>::iterator toChild = NULL;

	double est = 0.0;
	if (ptr->type == AND_NODE) {
		est = getChildEstimate(ptr);
	}

	auto ubPair = cumUbStack.top();
	cumUbStack.pop();
	auto estPair = cumEstStack.top();
	cumEstStack.pop();

	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// for upper priority, no need to exactly compute the priority to identify which child is the best
		if (sib->solved) continue; // do not consider a solved child
		if ( sib->type == AND_NODE ) {
//			if (sib->toChild != NULL) {
//				prior = sib->down_ub + sib->exact_val;
			// use random priority, difference between ub and estimate
			prior = sib->exact_val + sib->down_ub - sib->down_est;
//			} else {
//				prior = sib->ub;
//			}
		}
		else {
			// prt is AND ndoe
//			prior =  sib->down_ub + ptr->ub - sib->ub;
			// assume no removeWorst
			prior = sib->down_ub + ptr->ub - sib->ub - ( sib->down_est + est - sib->logEx );
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
//		ptr->down_est = toChild->down_est +
	} else {
		// if toChild is frontier node, the following equals toChild->ub
		ptr->down_ub = toChild->down_ub + toChild->exact_val;
	}
}*/
void randheur::findChild(tree<node>::iterator ptr, double cumUB, double cumEst) {
/*
 * find the best child for an UNsolved node and update its info accordingly
 * a solved child will not be considered
 * we assume no removeWorst
 * cumUB, cumEst are top-down equivalents of "topPath_UB", "topPath_Lb"
 */
	assert(ptr!=NULL); // should not be NULL
	if ( ptr->solved || (exploredTree.begin(ptr) == exploredTree.end(ptr) ) ) return;
//  note that "prior" here is the actual priority value instead of a relative quantity in the pure "ub" case.
	double bestPrior = -std::numeric_limits<double>::infinity();
	double prior = 0.0, prior_ub = 0.0, prior_lb = 0.0;
	tree<node>::iterator toChild = NULL;

	// caveat: ptr->est is NOT a simple aggregation of its children
	double est = 0.0, ub = 0.0;
	if (ptr->type == AND_NODE) {
		auto bdpair = getChildBounds(ptr);
		est = bdpair.first;
		ub = bdpair.second;
	}

	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		// not like the upper priority, we have to compute the exact priority
		if (sib->solved) continue; // do not consider a solved child
		if ( sib->type == AND_NODE ) {
//			if (sib->toChild != NULL) {
//				prior = sib->down_ub + sib->exact_val;
			// use random priority, difference between ub and estimate
//			prior = sib->exact_val + sib->down_ub - sib->down_est;

//			} else {
//				prior = sib->ub;
//			}
			prior_ub = cumUB + sib->exact_val + sib->down_ub;
			prior_lb = cumEst + sib->exact_val + sib->down_est;
		}
		else {
			// prt is AND ndoe
//			prior =  sib->down_ub + ptr->ub - sib->ub;
			// assume no removeWorst
//			prior = sib->down_ub + ptr->ub - sib->ub - ( sib->down_est + est - sib->logEx );
			prior_ub = cumUB + sib->down_ub + ub - sib->ub;
			prior_lb = cumEst + sib->down_est + est - sib->est;
		}
		prior = robustLogDiff(prior_ub, prior_lb);

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
//		ptr->down_ub = toChild->down_ub + ptr->ub - ptr->exact_val - toChild->ub;
		ptr->down_ub = toChild->down_ub + ub - toChild->ub;
		ptr->down_est = toChild->down_est + est - toChild->est;
	} else {
		// if toChild is frontier node, the following equals toChild->ub
		ptr->down_ub = toChild->down_ub + toChild->exact_val;
		ptr->down_est = toChild->down_est + toChild->exact_val;
	}
}
void randheur::findChild(tree<node>::iterator ptr) {
/*
 * find the best child for an UNsolved node and update its info accordingly
 * a solved child will not be considered
 * we assume no removeWorst
 * the priority here is: ub_path * ( n.ub - n.est )
 */
	assert(ptr!=NULL); // should not be NULL
	if ( ptr->solved || (exploredTree.begin(ptr) == exploredTree.end(ptr) ) ) return;
//  note that "prior" here is the actual priority value instead of a relative quantity in the pure "ub" case.
	double bestPrior = -std::numeric_limits<double>::infinity();
	double prior = 0.0, prior_ub = 0.0, prior_lb = 0.0;
	tree<node>::iterator toChild = NULL;

	// caveat: ptr->est is NOT a simple aggregation of its children
	double est = 0.0, ub = 0.0;
	if (ptr->type == AND_NODE) {
		auto bdpair = getChildBounds(ptr);
		est = bdpair.first;
		ub = bdpair.second;
	}

	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		// not like the upper priority, we have to compute the exact priority
		if (sib->solved) continue; // do not consider a solved child
		if ( sib->type == AND_NODE ) {
//			if (sib->toChild != NULL) {
//				prior = sib->down_ub + sib->exact_val;
			// use random priority, difference between ub and estimate
//			prior = sib->exact_val + sib->down_ub - sib->down_est;

//			} else {
//				prior = sib->ub;
//			}
//			prior_ub = cumUB + sib->exact_val + sib->down_ub;
//			prior_lb = cumEst + sib->exact_val + sib->down_est;
			prior_ub = sib->exact_val + sib->down_ub;
			prior_lb = sib->exact_val + sib->down_est;
		} else {
			// prt is AND ndoe
//			prior =  sib->down_ub + ptr->ub - sib->ub;
			// assume no removeWorst
//			prior = sib->down_ub + ptr->ub - sib->ub - ( sib->down_est + est - sib->logEx );
//			prior_ub = cumUB + sib->down_ub + ub - sib->ub;
//			prior_lb = cumEst + sib->down_est + est - sib->est;

			prior_ub = sib->down_ub + ub - sib->ub;
//			prior_lb = sib->down_est + est - sib->est;

			// Caveat: not exactly the same as that we used in AAAI-17 for lb
			// here we do not consider "est" from those branches, that's why we need only to consider bottom-up info.
			prior_lb = sib->down_est + ub - sib->ub;

		}
		prior = robustLogDiff(prior_ub, prior_lb);

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
//		ptr->down_ub = toChild->down_ub + ptr->ub - ptr->exact_val - toChild->ub;
		ptr->down_ub = toChild->down_ub + ub - toChild->ub;

		// Caveat: not exactly the same as that we used in AAAI-17 for lb
//		ptr->down_est = toChild->down_est + est - toChild->est;
		ptr->down_est = toChild->down_est + ub - toChild->ub;
	} else {
		// if toChild is a frontier node, the following equals toChild->ub
		ptr->down_ub = toChild->down_ub + toChild->exact_val;
		ptr->down_est = toChild->down_est + toChild->exact_val;
	}
}
int randheur::expandBest() {
	/*
	 * forward pass to expand the best node in OPEN
	 */
	++_nExpansion;
	if (verbose > 2) {
		std::cout << "expandBest: running " << _nExpansion << "-th round..." << std::endl;
	}
//	TODO: may not need to re-initialize?
//	std::fill(tuple.begin(), tuple.end(), -1); // re-use to accelerate, for safety, initialized to -1
	auto ptr = root;
	// forward pass, also update topPath_ub, topPath_lb
//	assert(cumUbStack.empty());
//	assert(cumEstStack.empty());
//	cumEstStack.push(std::make_pair(root, 0.0));
//	cumUbStack.push(std::make_pair(root, 0.0));

	while (ptr->toChild != NULL) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		if (DEBUG) printNode(*ptr); // debug
		ptr = ptr->toChild;
		if (ptr->type == AND_NODE) {
			tuple[ptr->X] = ptr->val;
			// The same as those of its OR parent
//			cumEstStack.push(std::make_pair(ptr, cumEstStack.top().second));
//			cumUbStack.push(std::make_pair(ptr, cumUbStack.top().second));

		} else {
//			auto eubPair = getSibBounds(ptr); // note that parPar->exact_val included here
//			cumEstStack.push(std::make_pair(ptr, eubPair.first + cumEstStack.top().second));
//			cumUbStack.push(std::make_pair(ptr, eubPair.second + cumUbStack.top().second));
		}
	}
	if (verbose > 2 || DEBUG) {
		std::cout << "\nBest node found!" << std::endl;
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
/*	auto& childrenList = rts;// if it is the virtual root
	if (ptr->id > -1) {
		childrenList = ascList[ptr->id];
	}*/
	 auto& childrenList = (ptr->id > -1)? ascList[ptr->id] : rts;
	//	std::vector<double> logQxCond(nvar, 0.0);

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
//	cur_lb = ptr->exact_val;
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
		mex::Factor dEst(child.X, 0.0);
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
		if (verbose > 2 || DEBUG) {
			int nMini = wmbUB.getNumberOfMiniBuckets(child.X);
//				isProp = isProp || (nMini > 1);
			if (nMini > 1) {
				std::cout
						<< "Bounds might be improved after expanding this node, no. of mini-buckets: "
						<< nMini << std::endl;
			}
		}
		// initialization
		ptrChild->solved = true;

//		auto unDone = getDescList(ptrChild->id);
//		std::list<mex::Var> desc, tipVars;
//		getDescVars(ptrChild->id, desc, tipVars);

		for (size_t v = 0; v < (child.X).states(); ++v) {
			node kid;
			kid.type = AND_NODE;
			kid.X = child.X;
			kid.id = child.id;
			kid.val = v;

			tuple[kid.X] = v;
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
			// solved only when all its children are solved
			if (!kid.solved) {
				ptrChild->solved = false;
			}
//			setRandHeur(kid, desc, tipVars, tuple, logQxCond);
			setRandHeur(kid, tuple);
//			dEst[v] = kid.down_est + kid.exact_val; // = kid.est
			dEst[v] = kid.est;
			// delete a solved node only when all its siblings are solved, i.e., its parent is solved
			// it serves as a book-keeping purpose, essential for two-stage sampling
			// append
			exploredTree.append_child(ptrChild, kid);
			++GSIZE;
		}
		if (!ptrChild->solved) ptr->solved = false;
//			ptrChild->lb = dLB.logsumexp();
		ptrChild->ub = dUB.logsumexp();
		ptrChild->logEx = dEst.logsumexp();
//			cur_lb += ptrChild->lb;
		cur_ub += ptrChild->ub;
		// if solved, remove all its children
		if (ptrChild->solved) {
			GSIZE -= exploredTree.number_of_children(ptrChild);
			exploredTree.erase_children(ptrChild);
		}
		// find best child
//		findChild(ptrChild);
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

	// we always have to update estimate
	// update bounds and estimates simultaneously, not necessary at this moment
//	updateEsts(ptr);
	//
	if (verbose > 2) {
		std::cout << "the best frontier node after expansion" << std::endl;
		printNode(*ptr);
	}
	//
//	auto eub_ptr = cumUbStack.top(); cumUbStack.pop();
//	auto est_ptr = cumEstStack.top(); cumEstStack.pop();
//	assert(eub_ptr.first == ptr);
//	assert(est_ptr.first == ptr);

	// now find best child for each child of (unsolved) ptr
	if (!ptr->solved) {
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// assume no re-expansion
//			findChild(sib, eub_ptr.second + ptr->ub - sib->ub, est_ptr.second + ptr->est - sib->est);
			findChild(sib);
		}
	}

	// caveat: make sure bounds have been fully updated before pruning
	ptr = pruneSolved(ptr); // lowest UNsolved

	// finally, backward update
	return backwardUpdate(ptr);
}
int randheur::backwardUpdate(tree<node>::iterator ptr) {
/*
 * backward update after expanding the best
 * at each node, update its best child and local estimate of the subproblem
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
	if ( root->solved ) {
		std::cout << "The root is solved!\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		return 1;
	}
	while (true) {
		//
//		auto eub_ptr = cumUbStack.top(); cumUbStack.pop();
//		auto est_ptr = cumEstStack.top(); cumEstStack.pop();
//		assert(eub_ptr.first == ptr);
//		assert(est_ptr.first == ptr);
//		findChild(ptr, eub_ptr.second, est_ptr.second);
		findChild(ptr);
		// go upward
		if ( ptr == root ) {
			break;
		}
		ptr = exploredTree.parent(ptr);
	}
	return 0;
}
/*void randheur::runSearch(const unsigned treeSizeLimit) {
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
}*/
/*int randheur::runSearch(const unsigned nnd, const unsigned treeSizeLimit) {

 * nnd: no. of nodes to expand
 *

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

//	LB = rt->lb;
	UB = root->ub;
	if (verbose > 1) {
		std::cout.precision(10);
		std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	return 0;
}*/
double randheur::logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done) {
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
double randheur::logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done, const std::list<mex::Var>& tipVars) {
/*
 * this is a generalization of the above "logF".
 * it can handle the case that we only sample until some level not till leaf nodes
 * Done: variables being sampled
 * tipVars: see "getDescVars"
 */
	if (DEBUG) {
		std::cout << "computing logF...\n";
	}
	assert(Done.size() < nvar+1); // debug
	if (Done.size() == nvar) return wmbUB.logP(config);
	// accumulate instantiated factors for all sampled variables
	double logVal = 0.0;
	for (const auto& X : Done) {
		logVal += wmbUB.heuristicTheta(X, config);
		if (DEBUG) {
			std::cout << wmbUB.heuristicTheta(X, config) << ",";
		}
	}
	if (DEBUG) {
		std::cout << "\n";
	}
	// accumulate heuristic values for those (non-leaf) tip vars
	for (const auto& X : tipVars) {
		logVal += wmbUB.heuristicIn(X, config);
		if (DEBUG) {
			std::cout << wmbUB.heuristicIn(X, config) << ",";
		}
	}
	if (DEBUG) {
		std::cout << "\ndone computing logF"<< std::endl;
	}
	return logVal;
}
tree<randheur::node>::iterator randheur::upperSampling(tree<node>::iterator ptr){
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
std::list<mex::Var> randheur::getDescList(int id) {
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
void randheur::getDescVars(int id, std::list<mex::Var>& desc, std::list<mex::Var>& tipVars) {
/*
 * a generalization of getDescList:
 * desc: all descendant variables to some level controlled by "_depth_lookahead"
 * tipVars: variables that do not have children in "desc"
 * note that those leaf variables in the pseudo tree that do not reach the depth limit are NOT included in tipVars
 * since their heuristic values are 0.0, which has no side-effect.
 */
//	std::list<mex::Var> desc;
	assert(desc.empty() && tipVars.empty());
	int max_depth = 0, depth = 0;
	mex::Var X;
	std::stack<mex::Var> stack;
	std::list<mex::graphModel::vindex> childrenList;
	if ( id  > -1 ) {
		childrenList = ascList[id];
		max_depth = _depthVec[id] + _depth_lookahead;
	}
	else {
		// if it is the virtual root
		childrenList = rts;
		max_depth = _depth_lookahead - 1; // _depth_lookahead >= 1, pseudo-root has depth -1
	}
	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
		stack.push(wmbUB.var(*it));
	}
	while(!stack.empty()) {
		X = stack.top();
		stack.pop();
		desc.push_back(X);
		// check depth to find tip vars
		// possibly some leaf variables with small depth not included
		depth = _depthVec[X.label()];
		if (depth == max_depth) {
			tipVars.push_back(X);
			continue;
		}
		//
		childrenList = ascList[X.label()]; //  label = id
		for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
			stack.push(wmbUB.var(*it));
		}
	}
}
void randheur::addSample(double logFx, const std::vector<double>& logQxCond) {
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

	// debug
	assert(wt < root->ub + epsilon);

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
	if (verbose > 1 || DEBUG) {
		std::cout << "NSP: " << NSP <<  ", Ex: " << ptr->logEx << ", Ex2: " << ptr->logEx2 << "\n";
		std::cout << "convEx: " << convEx << ", sumSqrUB: " << sumSqrUB << ", sumInvSqrUB: " << sumInvSqrUB << "\n";
		std::cout << "normalizedEx: " << normalizedEx << ", normalizedEx2: " << normalizedEx2 << "\n";
		std::cout << "sumInvUB: " << sumInvUB << ", hmUB: " << hmUB << std::endl;
	}
}
void randheur::addSampleToNode(node& n, double wt) {
/*
 * add one sample to a node
 */
	++n.nsp; // add count first
	// for OR nodes, only add one count and done
//	if ( (*ptr).type == OR_NODE ) return;
	if (n.nsp <= 1) {
		n.logEx = wt;
		n.logEx2 = 2 * wt;
	} else {
		int nsp = n.nsp;
		double logEx = n.logEx;
		double logEx2 = n.logEx2;
		// to be numerically stable
		// Ex  = (Ex * (samp-1))/samp + dEx/samp;
		n.logEx = logsumexp( { log(nsp - 1) + logEx, wt }) - log(nsp);
		// Ex2 = (Ex2 * (samp-1))/samp + dEx*dEx/samp;
		n.logEx2 = logsumexp( { log(nsp - 1) + logEx2, 2 * wt })
				- log(nsp);
	}
	// update necessary quantities
	// debug
//	if (DEBUG) {
//		std::cout << "NSP: " << NSP <<  ", Ex: " << n.logEx << ", Ex2: " << n.logEx2 << "\n";
//		std::cout << "convEx: " << convEx << ", sumSqrUB: " << sumSqrUB << ", sumInvSqrUB: " << sumInvSqrUB << "\n";
//		std::cout << "normalizedEx: " << normalizedEx << ", normalizedEx2: " << normalizedEx2 << "\n";
//		std::cout << "sumInvUB: " << sumInvUB << ", hmUB: " << hmUB << std::endl;
//		printNode(n);
//	}
}

std::pair<double, double> randheur::calcRootEBB(){
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
std::pair<double, double> randheur::calcRootEBB(double upperbound){
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
std::pair<double, double> randheur::calcHoeffding() {
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
std::pair<double, double> randheur::calcNormalizedEBB() {
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
void randheur::updateConvEstimate(double wt, double ub) {
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
void randheur::updateNormalizedEstimate(double wt, double ub) {
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
void randheur::twoStepSampling() {
/*
 * a customized version of two-step sampling (see the nips-17 DIS paper)
 */
	++_nsp_root;
	//	mex::vector<uint32_t> config(nvar, -1); // -1 for safety
//	auto& config = tuple; // accelerate, should be safe
//	std::vector<double> logQxCond(nvar, 0.0); // TODO: we may pre-locate the memory for this vector to accelerate?
//	std::list<mex::Var> Done; // list of all sampled variables
	// sample the path from root to leaf
	auto ptr = root, ptrPar = root;
	//values for solved variables, summation part of Rao-Blackwellisation
//	double solved_val = 0.0;
	// step 1 and 2
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
/*				if (ptr->type == AND_NODE) {
					// For solved AND node, the corresponding variable is IN config.
					solved_val += ptr->ub - ptr->exact_val; // exact_val will be computed logF(config, Done);
				} else {
					// For solved OR node, the corresponding variable is NOT in config.
					solved_val += ptr->ub;
				}*/
			} else {
				// an unsolved frontier node should be AND node due to the way we expand the best
				assert(ptr->type == AND_NODE);
				// sample from the frontier node according to the mixture proposal
/*				int ID = -1; // default, root id
				if (ptr->id > -1)
					ID = (ptr->X).label();
				// test, should be the same
				assert(ptr->id == ID);*/

				int ID = ptr->id;
				assert(ID > -1); // in stochastic lookahead case, we should not sample from the empty tree

				// variables for sampling
//				auto unDone = getDescList(ID);
//				std::list<mex::Var> unDone, tipVars;
//				getDescVars(ID, unDone, tipVars);
				// debug
				if (DEBUG) {
					std::cout << "ID: " << ID << "\n";
					std::cout << "UnDone: ";
					auto unDone = getDescList(ID);
					for (auto X : unDone) {
						std::cout << X.label() << ", ";
					}
					std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
				}
				// sample those variables in "unDone" conditioned on the path
				// config, logQxCond will be updated
//				wmbUB.conditionalMixtureSampling(ID, config, logQxCond, unDone);
				double wt = wmbUB.conditionalMixtureSampling(ID, tuple, rts, ascList, _depthVec, _adaMaxDepthVec[ID]);
				// add this sample to the corresponding frontier node
/*				double fx = logF(config, unDone, tipVars);
				double qx = 0.0;
				for (const auto& X : unDone) {
					qx += logQxCond[X.label()];
				}*/
//				addSampleToNode(*ptr, fx-qx);
				addSampleToNode(*ptr, wt);
				if (DEBUG) {
					printNode(*ptr);
				}
				// add a set of sampled variables
//				Done.insert(Done.end(), unDone.begin(), unDone.end());
				// update random heuristic for the frontier node and back up info. along the path
				updateRandHeur(ptr);
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
//				config[ptr->X] = ptr->val;
				tuple[ptr->X] = ptr->val;
//				logQxCond[ptr->id] = ptr->ub - ptrPar->ub; // ub of an OR parent is equal to sum of all AND children's ub
				// add a sampled variable
//				Done.push_back(ptr->X);
			} else {
				// push all OR children to stack
				for (auto sib = exploredTree.begin(ptr);
						sib != exploredTree.end(ptr); ++sib) {
					stack.push(sib);
				}
			}
		}
	}
	if (DEBUG) {
	// add the sample to the root
	// not necessary for our current purpose, comment out after debugging to save computational time.
/*		double logFx = logF(config, Done);
		logFx += solved_val; // add solved
		addSample(logFx, logQxCond);*/

	// debug
		printNode(*root);
		int nass = 0; // no. of assigned X
		for (auto v : tuple) {
			if (v >= 0 && v < nvar ) ++nass;
		}
	/*	std::cout << "no. of assigned X in config: " << nass
				<< "\nDone.size: " << Done.size()
				<< "\nlogFx: " << logFx
				<<"\nsolved_val: " << solved_val << std::endl;*/
/*		std::cout<< "config (logQxCond): ";
		for (int i=0; i<nvar; ++i) {
			std::cout << tuple[i] << " (" << logQxCond[i] << "), ";
		}*/
		std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
}
/*void randheur::doSampling() {
//	auto rt = exploredTree.begin();
	bool outNow = false;
	double outFrequency = 1;
	std::pair<double, double> prob;

	while(NSP < nSample) {
		twoStepSampling();
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
}*/
int randheur::doSampling(int nsp) {
/*
 * draw nsp samples from "two-step" sampling
 */
	std::pair<double, double> prob;
	for(int cnt = 0; cnt < nsp; ++cnt) {
		twoStepSampling();
//		if ( (NSP == _outFrequency)  || (verbose > 1) ) {
//			_outFrequency += _outFrequency;
////			power += 1.0;
////			prob = calcRootEBB();
//			if ( isInitUB ) prob = calcRootEBB(initUB);
//			else prob = calcRootEBB();
//			std::cout.precision(10);
////			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
//			std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
////			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
//			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
//			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
//			/***** more ways to aggregate estimates and bounds *****/
//			prob = calcHoeffding();
//			std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
//			prob = calcNormalizedEBB();
//			std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
//			/*********************************/
//			std::cout << "Tree size: " << GSIZE <<"\n";
//			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
//		}
//
//		if ( mex::timeSystem()-startProcess > timeBudget ) {
//			if ( isInitUB ) prob = calcRootEBB(initUB);
//			else prob = calcRootEBB();
//			std::cout.precision(10);
////			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
//			std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
////			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
//			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
//			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
//			/***** more ways to aggregate estimates and bounds *****/
//			prob = calcHoeffding();
//			std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
//			prob = calcNormalizedEBB();
//			std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
//			/*********************************/
//			std::cout << "Tree size: " << GSIZE <<"\n";
//			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
//			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//			return TIMEOUT;
//		}
	}
//	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//	if (NSP >= nSample) std::cout << "Reach sample size limit: " << nSample << std::endl;
/*	if (verbose > 1) {
		if ( isInitUB ) prob = calcRootEBB(initUB);
		else prob = calcRootEBB();
		std::cout.precision(10);
		std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
		std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
		std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
		**** more ways to aggregate estimates and bounds ****
		prob = calcHoeffding();
		std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
		prob = calcNormalizedEBB();
		std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
		*******************************
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}*/
	return 0;
}

/*void randheur::start() {
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
void randheur::start(unsigned treeSizeLimit) {
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
}*/
/*void randheur::start(const int nsp, const int nnd) {

 * the strategy is not quite clear yet
 * we may sample after each expansion or several expansions based on the complexity of drawing one sample and the complexity of one expansion?
 * how many samples should we draw in each sampling process
 * nsp: no. of samples to draw in each iteration
 * nnd: no. of nodes to expand in each iteration

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
}*/
void randheur::setRandHeur(node& n, std::list<mex::Var>& unDone, std::list<mex::Var>& tipVars,
		mex::vector<uint32_t>& config, std::vector<double>& logQxCond) {
/*
 * set random heuristics for newly-generated unSolved frontier nodes
 *
 * TODO: set const for solved nodes
 */
	if (n.solved) {
		n.est = n.ub;
		return;
	}
	assert(n.type == AND_NODE && n.toChild == NULL); // should be AND frontier node only

	// debug
	if (DEBUG) {
		std::cout << "\nrunning SetRandHeur..." << std::endl;
		std::cout << "before adding sample\n";
		printNode(n);
		std::cout << "unDone: ";
		for (auto X : unDone) std::cout << X.label() << " ";
		std::cout << "\n";
		std::cout << "tipVars: ";
		for (auto X : tipVars) std::cout << X.label() << " ";
		std::cout << "\n";
		std::cout << "config: ";
		for (auto x : config) std::cout << x << " ";
		std::cout << "\n";
		std::cout << "logQxCond: ";
		for (auto v : logQxCond) std::cout << v << " ";
		std::cout << std::endl;
	}

	//
	double fx = 0.0, qx = 0.0;
	for(int i=0; i<_ksp; ++i) {
//		std::fill(logQxCond.begin(), logQxCond.end(), 0.0);
		wmbUB.conditionalMixtureSampling(n.id, config, logQxCond, unDone);
//		fx = logF(config, unDone);
		fx = logF(config, unDone, tipVars);
		qx = 0.0;
		for (const auto& X : unDone) {
//			if (DEBUG) {
//				std::cout << X.label() << ", fx " << wmbUB.heuristicTheta(X, config) << "\n";
//			}
			qx += logQxCond[X.label()];
		}
		addSampleToNode(n, fx-qx);

		if (DEBUG) {
			std::cout << "fx (" <<  fx << ") - qx (" << qx << ") = " << fx-qx << "\n";
		}
		// frontier node, no removeWorst
		assert(n.down_ub + epsilon > fx-qx);
	}
	n.down_est = n.logEx; // estimate of the subproblem below this node
	n.est = n.logEx + n.exact_val;

	if (DEBUG) {
		// frontier node
		std::cout << "after adding sample\n";
		printNode(n);
		std::cout << "config: ";
		for (auto x : config) std::cout << x << " ";
		std::cout << "\n";
		std::cout << "logQxCond: ";
		for (auto v : logQxCond) std::cout << v << " ";
		std::cout << std::endl;

//		std::cout << "re-calculate down_ub: " << wmbUB.heuristicIn(n.X, config) << std::endl;
		assert(n.down_ub + epsilon > n.down_est);
	}
}
void randheur::setRandHeur(node& n, mex::vector<uint32_t>& config) {
/*
 * a faster version: avoid creating variables to sample in every expansion
 */
	if (n.solved) {
		n.est = n.ub;
		return;
	}
	assert(n.type == AND_NODE && n.toChild == NULL && n.id > -1); // should be AND frontier node only

	// debug
	if (DEBUG) {
		std::cout << "\nrunning SetRandHeur..." << std::endl;
		std::cout << "before adding sample\n";
		printNode(n);
//		std::cout << "unDone: ";
//		for (auto X : unDone) std::cout << X.label() << " ";
//		std::cout << "\n";
//		std::cout << "tipVars: ";
//		for (auto X : tipVars) std::cout << X.label() << " ";
//		std::cout << "\n";
		std::cout << "config: ";
		for (auto x : config) std::cout << x << " ";
//		std::cout << "\n";
//		std::cout << "logQxCond: ";
//		for (auto v : logQxCond) std::cout << v << " ";
		std::cout << std::endl;
	}
	//
//	double fx = 0.0, qx = 0.0;
	int max_depth = _adaMaxDepthVec[n.id];
	if (DEBUG) {
		std::cout << "max depth for sampling: " << max_depth << std::endl;
	}
	for(int i=0; i<_ksp; ++i) {
		// double wmbe::conditionalMixtureSampling(const int ID, MapType& config, const std::list<vindex>& rts,
//		const std::vector<std::list<vindex>>& ascList, const std::vector<int>& varDepthVec, const int maxDepth)
		double wt = wmbUB.conditionalMixtureSampling(n.id, config, rts, ascList, _depthVec, max_depth);
/*		fx = logF(config, unDone, tipVars);
		qx = 0.0;
		for (const auto& X : unDone) {
			qx += logQxCond[X.label()];
		}*/
//		addSampleToNode(n, fx-qx);
		addSampleToNode(n, wt);
/*		if (DEBUG) {
			std::cout << "fx (" <<  fx << ") - qx (" << qx << ") = " << fx-qx << "\n";
		}
		// frontier node, no removeWorst
		assert(n.down_ub + epsilon > fx-qx);*/
		if (DEBUG) {
			std::cout << "wt (" <<  wt << ")\n";
		}
		assert(n.down_ub + epsilon > wt);
	}
	n.down_est = n.logEx; // estimate of the subproblem below this node
	n.est = n.logEx + n.exact_val;

	if (DEBUG) {
		// frontier node
		std::cout << "after adding sample\n";
		printNode(n);
		std::cout << "config: ";
		for (auto x : config) std::cout << x << " ";
//		std::cout << "\n";
//		std::cout << "logQxCond: ";
//		for (auto v : logQxCond) std::cout << v << " ";
		std::cout << std::endl;
//		std::cout << "re-calculate down_ub: " << wmbUB.heuristicIn(n.X, config) << std::endl;
		assert(n.down_ub + epsilon > n.down_est);
	}
}
void randheur::updateRandHeur(tree<node>::iterator ptr) {
/*
 * update random heuristics and priorities after a global sampling
 * UnSolved frontier node only
 * backwardUpdate required
 */
	if (ptr->solved) return;
	assert(ptr->toChild == NULL);
	ptr->down_est = ptr->logEx;
	ptr->est = ptr->logEx + ptr->exact_val;

	// debug
	if (DEBUG) {
		// frontier node
		std::cout << "running updateRandHeur..." << std::endl;
		printNode(*ptr);
		assert(ptr->ub + epsilon > ptr->est);
	}

	// backwardUpdate required
	backwardUpdate(ptr);
}
std::pair<double,double> randheur::getSibBounds(tree<node>::iterator ptr) {
/*
 * accumulate upper bounds and estimates from sublings (deleted or non-deleted), edge weight also INCLUDED
 * if getAll is true, return all bounds including ptr itself, default, false
 * use the best bounds of its sublings
 * currently only for OR node
 */
	assert( ptr->type == OR_NODE );
	auto ptrPar = exploredTree.parent(ptr);
	double lb = ptrPar->exact_val;
	double ub = lb;
	for (auto sib=exploredTree.begin(ptrPar); sib != exploredTree.end(ptrPar); ++sib) {
		if (sib != ptr ) {
//			lb += sib->lb;
			lb += sib->est;
			ub += sib->ub;
		}
	}
	return std::make_pair(lb,ub);
}
std::pair<double, double> randheur::getChildBounds(tree<node>::iterator ptr) {
/*
 * accumulate estimates and upper bounds from its children
 * Caveat: "ptr->exact_val" NOT included (different from "getSibBounds")
 * needed only for AND_NODE
 */
	assert(ptr->type == AND_NODE);
	//	auto ptrPar = exploredTree.parent(ptr);
	// without removing worst, we have ptr->ub == ub + ptr->exact_val
	double est = 0.0, ub = 0.0;
	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		est += sib->est;
		ub += sib->ub;
	}
	return std::make_pair(est, ub);
}
double randheur::getChildEstimate(tree<node>::iterator ptr) {
/*
 * accumulate estimates of its children,
 * Caveat: "ptr->exact_val" NOT included (different from "getSibBounds")
 * currently only for AND node
 */
	assert( ptr->type == AND_NODE );
//	auto par = exploredTree.parent(ptr);
	double est = 0.0;
	for (auto sib=exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		est += sib->est;
	}
	return est;
}
double randheur::robustLogDiff(double ub, double lb) {
/*
 * compute log (exp( ub ) - exp( lb ) ) assuming large >= small
 */
	double res = lb - ub;
	res = (res < 0.0) ? res : 0.0;
	res = ub + log(1 - exp(res));
	return res;
}
int randheur::runSearch(const int nnd, const unsigned treeSizeLimit) {
/*
 * run search by doing expansion "nnd" times
 *
 */
//	double power = 1.0 + std::floor(std::log2(GSIZE+1)); // add

	if (verbose >1) {
		std::cout << "Running best-first search...\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	for(int cnt=0; cnt<nnd; ++cnt){
		if ( root->solved ) {
//			std::cout << "The root is solved, done!"<<std::endl;
			display();
			return SOLVED;
		}
		if ( GSIZE >= treeSizeLimit ) {
			std::cout << "Reach tree size limit: " << treeSizeLimit << std::endl;
/*			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;*/
			display();
			return MEMOUT;
		}
		if ( mex::timeSystem()-startProcess > timeBudget ) {
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
/*			UB = root->ub;
			std::cout.precision(10);
			std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;*/
			display();
			return TIMEOUT;
		}
/*		if ( GSIZE == int(pow(2.0, power))  &&  verbose > 0 ) {
			power += 1.0;
			display();
		}*/
		if ( _nExpansion >= _increment  &&  verbose > 0) {
			_increment += _increment; // to avoid expensive "log"
			display();
		}
		//
		expandBest();
	}

//	LB = root->lb;
//	UB = root->ub;
	if (verbose > 1) {
//		std::cout.precision(10);
//		std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB<<std::endl;
//		std::cout << "Tree size: " << GSIZE <<"\n";
//		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		display();
	}
	return 0;
}
int randheur::display() {
	UB = root->ub;
	std::cout.precision(10);
	std::cout<<"["<<mex::timeSystem()-startProcess<<"]: "<<std::setw(10)<<LB<<" < ln Z < "<<UB
			<< ", Zhat = " << root->est <<std::endl;
	std::cout << "Tree size: " << GSIZE << ", no. of node expansions: " << _nExpansion
			<< ", no. of full samples: " << _nsp_root;
	std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	return 0;
}
void randheur::start(const int nsp, const int nnd) {
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
//	std::cout << "Sample size limit: " << nSample << std::endl;
	std::cout << "Lookahead depth: " << _depth_lookahead << std::endl;
	std::cout << "nsp " << nsp << ", nnd " << nnd << ", ksp "<< _ksp << std::endl;
	assert(_depth_lookahead >= 0); // 0 for adaptive depth
	assert(nsp>0 || nnd>0);
//	const unsigned nrnd = 1; // tentatively
//	bool isMemOut = false;
	int status = 0;
	while (true) {
		// run search first until memory out
		// currently we do not have memory-limited version
//		if (!isMemOut && nnd>0) {
		status = runSearch(nnd, treeSizeLimit);
		if (status == MEMOUT) {
			std::cout << "memory out, quit!" << std::endl;
			break;
		}
		if (status == TIMEOUT) {
			std::cout << "time out, quit!" << std::endl;
			break;
		}
		if (status == SOLVED) {
			std::cout << "problem solved, quit!" << std::endl;
			break;
		}
		// do sampling
//			if (nsp>0) status = doSampling(nsp);
//		} else {
//			// if memory out, keep sampling until time out
//			// nsp=0 implies we only do sampling after search, thus we can use the latest UB for EBB
//			status = doSampling(std::numeric_limits<int>::max(), nsp > 0);
//			if (status == TIMEOUT) break;
//		}
		// sampling from the root
		doSampling(nsp);
	}
//	std::cout << "Total no. of node expansions: " << _nExpansion << std::endl;
}

//} /* namespace std */

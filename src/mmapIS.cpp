/*
 * mmapIS.cpp
 *
 *  Created on: Dec 1, 2017
 *      Author: qlou
 */

#include "mmapIS.h"

// namespace mex {

mmapIS::mmapIS(mex::wmbe& wmbeub, double ub, double delta, int vb, double tbgt, double mbgt,
		int nv, bool isSolved, int K, const mex::VarSet& maxVars, int nmax):
wmbUB(wmbeub), ascList(nv), UB(ub), DELTA(delta), verbose(vb), timeBudget(tbgt), memoryBudget(mbgt), nvar(nv),
exactHeur(nv, true), tuple(nv,-1), _initUB(ub), _isExactHeur(isSolved), _K(K), _logK(log(K)), _varTypeVec(nv, SUM_VAR),
_nMax(nmax), mapConfig(nv, 0), isLeafVar(nv, false), dfsConfig(nv, -1), _locPair(nv, std::make_pair(nv,nv)),
_mapSpaceSizeVec(nv, -std::numeric_limits<double>::infinity()) {
	startProcess = mex::timeSystem();
	order = wmbeub.getOrder();
	priority = wmbeub.getPriority();
	GSIZE = 0;
	LB = -std::numeric_limits<double>::infinity();
	NSP = 0;
	_mapSpaceSize = 0.0;
	for (auto it = maxVars.begin(); it != maxVars.end(); ++it) {
		_varTypeVec[ (*it).label() ] = MAX_VAR;
		_mapSpaceSize += log( (*it).states() );
	}
	if (DEBUG) std::cout << "mapSpaceSize: " << exp(_mapSpaceSize) << std::endl;

	buildPseudoTree();
	setRootConfig();
	setExactHeur();
	setMapSpaceSize();
//	setDepth();
}
mmapIS::~mmapIS() {
	// TODO Auto-generated destructor stub
	std::cout<<"==== Done Search+Sampling ====\n";
}
void mmapIS::setRootConfig() {
	// initialize root
	// this is a virtual root node, nodes without parents are children of this virtual root node
	node vtRoot; // virtual root
	vtRoot.nodeType = AND_NODE; // virtual AND_NODE node
	vtRoot.varType = MAX_VAR; // a virtual MAX varialbe
	vtRoot.id = -1; // flag
	vtRoot.exact_val = 0.0;
//	vtRoot.lb = LB;
	vtRoot.ub = UB;
	vtRoot.mapSpaceLeft = _mapSpaceSize;
	//
	exploredTree.insert(exploredTree.begin(), vtRoot);
	root = exploredTree.begin();
	GSIZE = 1;
	std::cout<<"==== Root configuration set ====\n";
	if (DEBUG) printNode(vtRoot);
}
void mmapIS::buildPseudoTree() {
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

	/**********generate the pseudo tree structure and store it in one vector**************/
	/* note that _locPair is initialized to "(nvar,nvar)"
	 * a key observation:
	 * 	let "id" be the first child of "par" added to the stack,
	 * 	then _locPair[par].second = _locPair[id].second
	 */
	_pseudotree.clear();
	std::stack< mex::graphModel::vindex > stack;
	for (auto id : rts) {
		// the last child popped out is its first child in the childList
		stack.push(id);
	}
	while(!stack.empty()) {
		auto id = stack.top();
		stack.pop();
		_pseudotree.push_back(wmbUB.var(id));
		_locPair[id].first = _pseudotree.size()-1; // location
		//
		auto childList = ascList[id];
		if (childList.empty()) {
			// a frontier variable, track back
			_locPair[id].second = _locPair[id].first+1;
			auto par = parents[id];
			if (par < 0 || par >= nvar) continue; // this frontier node has no parent
			bool isFirstChild = ( *(ascList[par].begin()) == id ); // whether this is the first child pushed into the stack
			while (isFirstChild){
				_locPair[par].second = _locPair[id].second; // key observation
				id = par;
				par = parents[id];
				if (par < 0 || par >= nvar) break;
				isFirstChild = ( *(ascList[par].begin()) == id );
			}
		} else {
			for (auto i : childList) stack.push(i);
		}
	}
	assert(_pseudotree.size() == nvar);
	std::cout<<"==== Pseudo tree built ====\n";

	/***********************************/
	if (DEBUG) {
		std::cout << "no. of nodes without parents: " << rts.size() << "/" << nvar << "\n";
		for (size_t i = 0; i < ascList.size(); ++i) {
			auto X = wmbUB.var(i);
			std::string str = _varTypeVec[i]==MAX_VAR ? "MAX " : "SUM ";
			std::cout << str << "variable " << i << " with " << X.states() << " states has " << (ascList[i]).size() << " children: [";
			for (auto jj : ascList[i]) std::cout << jj << ", ";
			std::cout << "]" << std::endl;
		}
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "variables in the pseudo tree vector: \n";
		for (auto X : _pseudotree) {
			auto id = X.label();
			auto loc = _locPair[id].first;
			assert(_pseudotree[loc].label() == id);
			std::cout << id << ", ";
		}
		std::cout << "\n a variable and its last descendant locations in the pseudo tree vector: \n";
		for (size_t i = 0; i < nvar; ++ i) {
			std::cout << i << ": (" << _locPair[i].first << ", " << _locPair[i].second << ")\n";
		}
		std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
}
void mmapIS::setMapSpaceSize() {
/*
 * set the MAP space size below each MAP variable
 * Caveat: set "_varTypeVec" before running this function
 */
	//for each leaf MAP variable or SUM variable, the
	std::stack<mex::graphModel::vindex> stackFwd;
	 // stackBwd holds all the variables such that each variable has its descendants on top of it
	auto stackBwd = stackFwd;
	for (auto id : rts) stackFwd.push(id);
	//
	while (!stackFwd.empty()) {
		auto top = stackFwd.top();
		stackFwd.pop();
		stackBwd.push(top);
		if (_varTypeVec[top] == MAX_VAR) _mapSpaceSizeVec[top] = 0.0; // re-initialization
		for (auto id : ascList[top]) {
			// we only need to push MAX children
			if (_varTypeVec[id] == MAX_VAR)
				stackFwd.push(id);
		}
	}
	// compute
	auto parents = wmbUB.getPseudotree();
	while (!stackBwd.empty()) {
		auto top = stackBwd.top();
		auto X = wmbUB.var(top);
		stackBwd.pop();
		auto par = parents[top];
		if ( par >= 0 && par < nvar )  {
			_mapSpaceSizeVec[par] += _mapSpaceSizeVec[top] + log(double(X.states()));
		}
	}
	if (DEBUG) {
		std::cout << "MAP space size below each variable (0 if a SUM variable):\n";
		for (size_t id=0; id<nvar; ++id) {
			std::cout << "["<< id << ", " << exp(_mapSpaceSizeVec[id]) << "], ";
		}
		std::cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
}
/*void mmapIS::setDepth() {

 * set depth vector via depth-first search
 * 0-based: a variable without a parent is of depth 0
 *
 * _descList: for any variable X, all its descendants in the pseudo tree are consecutively located in the list.;
 * i.e., positions between X and any X's descendant are for X's descendants only

	_descList.clear();
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
		// add X to the list when it's popped
		auto X = wmbUB.var(id);
		_descList.push_back(X);
		_descPosVec[id] = _descList.size()-1;
		//
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
		assert(_depthVec.size() == nvar && _descList.size() == nvar);
		std::cout << "Depth of each variable:\n";
		for (int i=0; i < _depthVec.size(); ++i) {
			std::cout << i << " (" << _depthVec[i] << "), ";
		}
		std::cout << std::endl;
	}
}*/
void mmapIS::setExactHeur() {
/*
 * set exactHeur via breadth-first search from leaf nodes of the pseudo tree
 * exactHeur tells whether the heuristics below that node is exact or not
 * note that exactHeur is initialized to be true
 */
	// Caveat: this function should still work for the augmented model
	// because it only checks whether partitioning exists in a bucket

	std::vector<size_t> nVis(nvar,0); // no. of children visited so far
	auto parents = wmbUB.getPseudotree();
	std::queue<mex::graphModel::vindex> queue; // variable id
	for (int i=0; i<nvar; ++i) {
		// add leaf nodes
		if (ascList[i].empty()) {
			queue.push(i);
			isLeafVar[i] = true;
		}
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
	if (verbose > 2 || DEBUG) {
		std::cout<<"~~~~~~~~~~~~~~~Printing exactHeur~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		for (int i=0; i<nvar; ++i) {
			std::cout << "[" << i << "," << exactHeur[i] << "], ";
			assert(nVis[i] == ascList[i].size());
		}
		std::cout<<"\n~~~~~~~~~~~~~~~Done printing~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	}
}
// helper functions
double mmapIS::logsumexp(std::vector<double>& vec) {
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
double mmapIS::logsumexp(std::list<double>& vec) {
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
double mmapIS::logsumexp(std::initializer_list<double> vec) {
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
void mmapIS::printNode(const node& n) {
/*
 * print node info for debugging
 */
	std::string vt_str = (n.varType == MAX_VAR)? "MAX " : "SUM ";
	std::string nt_str = (n.nodeType == AND_NODE)? "AND " : "OR ";
	std::string sol_str = (n.solved)? "Solved " : "Not solved ";
	std::string sumsol_str = (n.sumSolved)? "sumSolved" : "Not sumSolved";
	std::string node_str = (n.toChild == NULL)? "frontier " : "internal ";

	std::cout << "~~~~~~~~~~~~~~\n" ;
	std::cout << vt_str <<  nt_str  << node_str << sol_str << sumsol_str << ", id " << n.id << ", val " << n.val << "\n";
	std::cout << "ub " << n.ub << ", exact_val " << n.exact_val << ", wSum " << n.wSum << "\n";
	std::cout << "ub/K " << n.ub/_K << ", exact_val/K " << n.exact_val/_K << "\n";
	std::cout << "down_ub " << n.down_ub << ", downEst " << n.downEst << "\n";
	std::cout << "nsp " << n.nsp <<  ", mapSpaceLeft (log2) " << n.mapSpaceLeft/log(2) << "\n";
	std::cout << "mapDesc: ";
	for (auto& v : n.mapDesc) std::cout << "(" << v.first << ", " << v.second << "), ";
	std::cout << "\n~~~~~~~~~~~~~~" << std::endl ;
/*	if (n.nodeType == AND_NODE) {
		assert(n.ub + epsilon > n.exact_val + n.downEst );
	} else {
		assert(n.ub + epsilon > n.downEst);
	}*/
}
void mmapIS::checkSolved(node& n) {
/*
 * Check whether the SUM problem below has been solved
 * For SUM nodes, it's equivalent to check its "solvedness"
 * works for both AND_NODE and OR_NODE nodes
 * Caveat: only apply to newly generated nodes
 */
	if (n.id > -1) {
		// if not the virtual root
		n.sumSolved = exactHeur[n.id]; // for leaf node in pseudo tree, this is always true
		n.solved = isLeafVar[n.id]; // leaf variables must be solved no matter MAX or SUM
		//
		if (n.varType == SUM_VAR) n.solved = n.sumSolved; // for SUM nodes, it's equivalent
	}
//
	if (n.solved) assert(n.sumSolved);
//	else {
//		solved = (n.ub - n.lb <= solvedThresh) ;
//	}
//	return solved;
}
tree<mmapIS::node>::iterator mmapIS::pruneSolved(tree<node>::iterator ptr) {
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
		/*		if (ptr->varType == MAX_VAR && ptr->nodeType == OR_NODE) {
		 // for an OR-MAx node, solving its child with highest upper bound is enough to ensure its solvedness
		 tree<node>::iterator largestChild = NULL;
		 double largestUpperBound = -std::numeric_limits<double>::infinity();
		 for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		 if (sib->ub > largestUpperBound) {
		 largestChild = sib;
		 largestUpperBound = sib->ub;
		 }
		 }
		 // solved iff its largest child is solved
		 if (largestChild != NULL) {
		 solved = largestChild->solved;
		 }
		 } else {*/
		// for other cases
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr);
				++sib) {
			// a node is solved only when all its children are solved
			if (!sib->solved) {
				solved = false;
				break;
			}
		}
//		}
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
	if ((*ptr).nodeType == AND_NODE) {
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
		if ( highestSolved->nodeType == AND_NODE ) {
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

	if (verbose > 2 || DEBUG) {
		std::cout << "pruning done: " << GSIZE_pre - GSIZE << " nodes deleted!\n" ;
	}
//	return ptr;
	// return parent if not root
	if ( highestSolved == root ) return root;
	return exploredTree.parent(highestSolved);
}
void mmapIS::updateBounds(tree<node>::iterator ptr) {
	// re-calculate bounds from existing children and exact_val
	// bound update is different for the augmented model
	if (verbose > 2) {
		std::cout << "now update bounds: " << std::endl;
	}
	assert( ptr != NULL );
	ptr = exploredTree.parent(ptr);
	auto top = exploredTree.parent(root);
	double cur_ub = 0.0;

	// bound update is different only for the "AND-MAX parent + OR-SUM child" case
	while( ptr != top ) {
		// update bounds by re-calculation.
		// more time-consuming than simple update, but more numerically stable
		if (ptr->nodeType == AND_NODE) {
			cur_ub = ptr->exact_val;
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if (ptr->varType != sib->varType) {
					// MAX vs SUM
					cur_ub += _K * sib->ub;
				} else {
					cur_ub += sib->ub;
				}
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
//		if ( cur_ub < ptr->ub ) {
		// more stable
		if ( cur_ub + epsilon < ptr->ub ) {
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

tree<mmapIS::node>::iterator mmapIS::updateBoundsPlusPruneSolved(tree<node>::iterator ptr) {
/*
 * Merge bound update and node pruning into one function;
 * act as "updateBound" + "pruneSolved"
 * Must run from the newly expanded node!!!
 * -- make sure this node is NOT a frontier node!!!
 * Caveat: a solved node will be deleted only when all its siblings are solved
 * this is crucial (to Rao-blackwellisation) here
 * return the highest UNsolved node
 *
 * For MAX nodes, it can be pruned only when the SUM problem and the MAX problem below are both solved,
 * -- this is how "solvedness" defined for MAX nodes,
 * which requires bound update in this pruning process
 */
	assert( ptr != NULL && exploredTree.begin(ptr) != exploredTree.end(ptr) );
	if ( DEBUG ) {
		std::cout << "running updateBoundsPlusPruneSolved...\n";
		printNode(*ptr);
	}
	tree<node>::iterator highestSolved = NULL;
	auto lowestUnSolved = ptr;
//	ptr = exploredTree.parent(ptr);
	auto top = exploredTree.parent(root);
	double cur_ub = 0.0;
	// Caveat:
	// bound update is different for the "AND-MAX parent + OR-SUM child" case
	// pruning is different for the OR-MAX case
	while( ptr != top ) {
		// update bounds by re-calculation.
		// time-consuming compared to a simple update by computing the difference, but numerically stable
		bool solved = true;
		bool sumSolved = true;
		if (ptr->nodeType == AND_NODE) {
			cur_ub = ptr->exact_val;
			ptr->mapDesc.clear();
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if (ptr->varType != sib->varType) {
					// MAX vs SUM
					cur_ub += _K * sib->ub;
				} else {
					cur_ub += sib->ub;
				}
				if (!sib->solved) solved = false;
				if (!sib->sumSolved) sumSolved = false;
			}
		}
		else {
			// OR_NODE node
//			std::vector<double> ub_vec {(*ptr).exact_val}; // ptr->exct_val = -inf for OR nodes
			std::list<double> ub_vec; // no initialization, actually safer
			// for the OR-MAX case
			tree<node>::iterator largestChild = NULL; // child with the highest upper bound
			double largestUpperBound = -std::numeric_limits<double>::infinity();
			// "ptr" should have at least one child
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				ub_vec.push_back( sib->ub );
				//
				if (!sib->solved) solved = false;
				if (!sib->sumSolved) sumSolved = false;
				// if OR-MAX node, find its largest child, and check its solvedness
				if (ptr->varType == MAX_VAR && sib->ub >= largestUpperBound) {
					largestChild = sib;
					largestUpperBound = sib->ub;
				}
			}
			// the evaluation order matters: && should perform short-circuit
			if ( (ptr->varType == MAX_VAR) && (largestChild != NULL) && (largestChild->solved) ) {
				solved = true;
				cur_ub = largestChild->ub;
				// update its best descendants right now because we will remove "ptr"'s children soon
/*				ptr->mapDesc = largestChild->mapDesc;
				ptr->mapDesc.push_front(std::make_pair(largestChild->id,largestChild->val));*/
				// faster? reduce one memory allocation step
				ptr->mapDesc.clear();
				ptr->mapDesc.push_back(std::make_pair(largestChild->id,largestChild->val));
				ptr->mapDesc.splice(ptr->mapDesc.end(), largestChild->mapDesc);
				assert(largestChild->mapDesc.empty());
				if (DEBUG) {
					std::cout << "pruning since an OR-MAX node's largest child has been solved!" << std::endl;
				}
			} else {
				cur_ub = logsumexp(ub_vec);
			}
		}
		/**************************************************************************************************/
		ptr->solved = solved;
		ptr->sumSolved = sumSolved;
		if (ptr->solved) {
			ptr->sumSolved = true; // make sure "solved" => "sumSolved". This is useful in the sampling step.
			highestSolved = ptr;
			// since those MAX children will be deleted later, save their MAP configurations in their parent
			// not needed for OR-MAX nodes, already updated
			if (ptr->varType == MAX_VAR && ptr->nodeType == AND_NODE) {
				ptr->mapDesc.clear(); // is it safe? should be. since this node is not a frontier node
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					// debug
					if (sib->varType == SUM_VAR) {
//						assert(sib->mapDesc.empty());
						continue;
					}
//					for (auto pair : sib->mapDesc) ptr->mapDesc.push_back(pair);
					ptr->mapDesc.splice(ptr->mapDesc.end(), sib->mapDesc); // faster, empty sib->mapDesc at the same time
					// debug
					assert(sib->mapDesc.empty());
				}
			}
		}

//		if ( cur_ub < ptr->ub ) {
		// more stable
		if ( cur_ub + epsilon  < ptr->ub ) {
			ptr->ub = cur_ub;
		} else {
			// no need to move upward if no solved + no bound update
			if (!ptr->solved) break;
		}
		ptr = exploredTree.parent(ptr);
	}
	// update bounds for the root
	UB = root->ub;
	// return the input node if no pruning happened
	if (highestSolved == NULL) return lowestUnSolved;
	// finally, delete those solved descendants
	auto GSIZE_pre = GSIZE;
	GSIZE -= exploredTree.size(highestSolved); // size of the subtree rooted at the given node (including that node)
	++GSIZE; // do not count highestSolved
	exploredTree.erase_children(highestSolved);
	//
	if (verbose > 2 || DEBUG) {
		std::cout << "pruning done: " << GSIZE_pre - GSIZE << " nodes deleted!\n" ;
	}
//	return ptr;
	// return parent if not root
	if ( highestSolved == root ) return root;
	return exploredTree.parent(highestSolved);
}
void mmapIS::propMapSpace(tree<node>::iterator ptr) {
/*
 * propagate the remaining MAP space for a given MAX node
 * run after pruning "solved"
 */
	if (ptr->varType == SUM_VAR) {
		// no need to propagate
		return;
	}
//	ptr = exploredTree.parent(ptr); // thus, ptr can not be a frontier node
	double mapSpaceLeft = 0.0;
	auto top = exploredTree.parent(root);
	while (ptr != top) {
		if (ptr->solved) {
			ptr->mapSpaceLeft = 0.0; // exactly one configuration for solved nodes
		} else {
			if (ptr->nodeType == AND_NODE) {
				// AND-MAX
//				ptr->mapSpaceLeft = 0.0;
				mapSpaceLeft = 0.0;
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					if (sib->varType == MAX_VAR) {
						if (sib->solved) sib->mapSpaceLeft = 0.0; // be safe
//						ptr->mapSpaceLeft += sib->mapSpaceLeft;
						mapSpaceLeft += sib->mapSpaceLeft;
					}
				}
			} else {
				// OR-MAX
				std::list<double> spList; // no need to initialize, safer
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					if (sib->solved) sib->mapSpaceLeft = 0.0; // be safe
					spList.push_back(sib->mapSpaceLeft);
				}
//				ptr->mapSpaceLeft = logsumexp(spList);
				mapSpaceLeft = logsumexp(spList);
			}
			// be efficient
			if (mapSpaceLeft + epsilon < ptr->mapSpaceLeft) {
				ptr->mapSpaceLeft = mapSpaceLeft;
			} else {
				// if no improvement, stop here
				return;
			}
		}
		ptr = exploredTree.parent(ptr);
	}
}
void mmapIS::findChild(tree<node>::iterator ptr) {
/*
 * find the best child for an UNsolved node
 * and update its info accordingly
 * a solved child will not be considered because we do not need to expand it in the future
 * Works for both MAX and SUM nodes
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
		if ( ptr->nodeType == OR_NODE ) {
				prior = sib->down_ub + sib->exact_val;
		}
		else {
			// ptr is AND node
			// Caveat for the case AND-MAX parent + OR-SUM child
			if (ptr->varType != sib->varType) {
				prior = _K*sib->down_ub + ptr->ub - _K*sib->ub;
			} else {
				prior =  sib->down_ub + ptr->ub - sib->ub;
			}
		}
		if ( (toChild == NULL) || (prior > bestPrior) ) {
			toChild = sib;
			bestPrior = prior;
		}
	}
	// if not solved
	assert(toChild != NULL);
	ptr->toChild = toChild;
	// not including exact_val
	if (ptr->nodeType == AND_NODE) {
		// Caveat for the case AND-MAX parent + OR-SUM child
		if (ptr->varType != toChild->varType) {
			ptr->down_ub = _K * toChild->down_ub + ptr->ub - ptr->exact_val - _K * toChild->ub;
		} else {
			ptr->down_ub = toChild->down_ub + ptr->ub - ptr->exact_val - toChild->ub;
		}
	} else {
		// if toChild is frontier node, the following equals toChild->ub
		ptr->down_ub = toChild->down_ub + toChild->exact_val;
	}
}
int mmapIS::expandBest() {
	/*
	 * forward pass to expand the best node in OPEN
	 */
	if (verbose > 2 || DEBUG) {
		std::cout << "running expandBest..." << std::endl;
	}
	// need not re-initialize, should be safe
//	std::fill(tuple.begin(), tuple.end(), -1); // re-use to accelerate, for safety, initialized to -1
	auto ptr = root;  // the root is an AND & MAX node
	// forward pass
	while (ptr->toChild != NULL) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
//		if (DEBUG) printNode(*ptr); // debug
		//
		ptr = ptr->toChild;
		if (ptr->nodeType == AND_NODE) {
			tuple[ptr->X] = ptr->val;
		}
	}
	if (verbose > 2 || DEBUG) {
		std::cout << "Best node found!" << std::endl;
		printNode(*ptr);
	}

	if (ptr->solved) {
		std::cout << "The best frontier node has been solved, quit!" << std::endl;
		assert(root->solved);
		return SOLVED;
	} else {
		if (DEBUG && ptr->sumSolved) std::cout << "The SUM problem below the best frontier node has been solved." << std::endl;
	}

	assert(ptr->nodeType == AND_NODE); // without removal, the frontier must be AND

/*
 //	if (ptr->nodeType == AND_NODE) {
//		double cur_lb=0.0, cur_ub = 0.0;
	double cur_ub = 0.0;
	std::list<mex::graphModel::vindex> childrenList = (ptr->id > -1)? ascList[ptr->id] : rts; // save memory "&"
	if (ptr->id > -1) {
		childrenList = ascList[ptr->id];
	} else {
		// if it is the virtual root
		childrenList = rts;
	}
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
	ptr->sumSolved = true;
	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
		node child;
		child.id = *it;
		child.X = wmbUB.var(child.id);
		child.nodeType = OR_NODE;
		child.varType = _varTypeVec[child.id];
		child.exact_val = -std::numeric_limits<double>::infinity(); // for OR_NODE node
		child.mapSpaceLeft = _mapSpaceSizeVec[child.id] + log((child.X).states());
		// append to the explored tree
		auto ptrChild = exploredTree.append_child(ptr, child);
		++GSIZE;
		//
		mex::Factor dUB(child.X, 0.0);
//		mex::Factor dLB(child.X, 0.0);
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
		if (verbose > 2) {
//				int nMini = wmbUB.duplicateInBucket(child.X, tuple);
			int nMini = wmbUB.getNumberOfMiniBuckets(child.X);
//				isProp = isProp || (nMini > 1);
			if (nMini > 1) {
				std::cout
						<< "Bounds are possibly improved after expanding this node, no. of mini-buckets: "
						<< nMini << std::endl;
			}
		}
		// initialization
		ptrChild->solved = true;
		ptrChild->sumSolved = true;
		// to track the largest child for a MAX node
		node largestKid; // child with the highest upper bound
		double largestUpperBound = -std::numeric_limits<double>::infinity();
		//
		for (size_t v = 0; v < (child.X).states(); ++v) {
			node kid;
			kid.varType = child.varType;
			kid.nodeType = AND_NODE;
			kid.X = child.X;
			kid.id = child.id;
			kid.val = v;
			kid.mapSpaceLeft = _mapSpaceSizeVec[kid.id];
			tuple[child.X] = v;
			// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
			// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node

			// both MAX and SUM variables apply.
			// note that factors that only involves MAX variables have already replicated during WMB construction
			kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
			// actually, we can speed up this evaluation by computing each factor only once
			// for a leaf node, leaf_ub, leaf_lb are its ub, lb
//			kid.down_ub = wmbUB.heuristicIn(kid.X, tuple);
			// heuristics are differently calculated for MAX and SUM nodes
			kid.down_ub = wmbUB.heuristicInAug(kid.X, tuple, _K, kid.varType == MAX_VAR);
			// no change, same for both MAX and SUM variables
			kid.ub = dUB[v] = kid.down_ub + kid.exact_val;
			// whether solved
//			kid.solved = checkSolved(kid);
			// those children actually have the same "solvedness" due to the way how "checkSolved" works
			// that's why even for the OR MAX case we do not need to track the largest child here
			checkSolved(kid);
			if (!kid.solved) ptrChild->solved = false;
			if (!kid.sumSolved) ptrChild->sumSolved = false;

			// delete a solved node only when all its siblings are solved, i.e., its parent is solved
			// it is serve as a book-keeping purpose, essential for two-step sampling
			// append
			exploredTree.append_child(ptrChild, kid);
			++GSIZE;
			if (kid.ub >= largestUpperBound ) {
				largestUpperBound = kid.ub;
				largestKid = kid;
			}
		}
		if (!ptrChild->solved) ptr->solved = false;
		if (!ptrChild->sumSolved) ptr->sumSolved = false;
//			ptrChild->lb = dLB.logsumexp();
		ptrChild->ub = dUB.logsumexp();
//			cur_lb += ptrChild->lb;
		// Caveat for the case: AND-MAX parent + OR-SUM child: raised to the K-th power
		if (ptr->varType != ptrChild->varType) {
			cur_ub += _K * ptrChild->ub;
		} else {
			cur_ub += ptrChild->ub;
		}
		// for OR-SUM nodes, if "sumSolved", it must be "solved", remove all its children
		// while for OR-MAX nodes, "sumSolved" does not ensure its "solvedness"
		// remove a node's all children only when it's marked "solved".
		if (ptrChild->solved) {
			// for a SUM node, we can delete all its children
			GSIZE -= exploredTree.number_of_children(ptrChild);
			exploredTree.erase_children(ptrChild);
			if (ptr->varType == MAX_VAR) {
				// for a MAX node, we have to keep the best configuration of MAP descendants
				ptrChild->mapDesc.push_back(std::make_pair(largestKid.id,largestKid.val));
				// we will update "mapSpaceLeft" for its ancestors later on
				ptr->mapSpaceLeft -= ptrChild->mapSpaceLeft;
			}
		}
		// find best child
		findChild(ptrChild);
	}
	// more stable
	if ( ptr->ub > cur_ub + epsilon || ptr == root) {
		isProp = true;
	}

	isProp = (ptr->ub > cur_ub + epsilon) || (ptr == root) || (ptr->solved);

	// propagate bounds if necessary
	if (isProp) {
//		if (isProp && !ptr->isReExpand) {
//			ptr->lb = cur_lb;
		ptr->ub = cur_ub;
//		updateBounds(ptr);
		ptr = updateBoundsPlusPruneSolved(ptr); // do both simultaneously, be cautious!!!
	}
	if (verbose > 2) {
		std::cout << "the best frontier node after expansion" << std::endl;
		printNode(*ptr);
	}

	// propagate the MAP space reduction info
	// TODO: update when necessary?
	propMapSpace(ptr);
*/

//	 caveat: make sure bounds have been fully updated before pruning
//	ptr = pruneSolved(ptr); // lowest UNsolved

	/******************************************************/
	ptr = expandOneFrontierNode(ptr, tuple);
	return backwardUpdate(ptr);
}
tree<mmapIS::node>::iterator mmapIS::expandOneFrontierNode(tree<node>::iterator ptr, mex::vector<uint32_t>& config) {
/*
 * Expand one frontier node:
 *  do necessary bound update and pruning and other information propagation
 */
	assert(ptr != NULL);
	assert(!ptr->solved);
	if (DEBUG) {
		std::cout << "Running expandOneFrontierNode...:\n";
		printNode(*ptr);
	}
	++_nExpansion;
	bool isProp = false;
	double cur_ub = 0.0;
	const auto& childrenList = (ptr->id > -1) ? ascList[ptr->id] : rts; // save memory "&"
	cur_ub = ptr->exact_val;
	// initialization
	ptr->solved = true;
	ptr->sumSolved = true;

	for (auto it = childrenList.cbegin(); it != childrenList.cend(); ++it) {
		node child;
		child.id = *it;
		child.X = wmbUB.var(child.id);
		child.nodeType = OR_NODE;
		child.varType = _varTypeVec[child.id];
		child.exact_val = -std::numeric_limits<double>::infinity(); // for OR_NODE node
		child.mapSpaceLeft = _mapSpaceSizeVec[child.id] + log((child.X).states());
		// append to the explored tree
		auto ptrChild = exploredTree.append_child(ptr, child);
		++GSIZE;
		//
		mex::Factor dUB(child.X, 0.0);
		if (verbose > 2 || DEBUG) {
			int nMini = wmbUB.getNumberOfMiniBuckets(child.X);
			//				isProp = isProp || (nMini > 1);
			if (nMini > 1) {
				std::cout << "Bounds are possibly improved after expanding this node, no. of mini-buckets: "
						<< nMini << std::endl;
			}
		}
		// initialization
		ptrChild->solved = true;
		ptrChild->sumSolved = true;
		// to track the largest child for a MAX node
		node largestKid; // child with the highest upper bound
		double largestUpperBound = -std::numeric_limits<double>::infinity();
		//
		for (size_t v = 0; v < (child.X).states(); ++v) {
			node kid;
			kid.varType = child.varType;
			kid.nodeType = AND_NODE;
			kid.X = child.X;
			kid.id = child.id;
			kid.val = v;
			kid.mapSpaceLeft = _mapSpaceSizeVec[kid.id];
			config[child.X] = v;
			// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
			// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node

			// both MAX and SUM variables apply.
			// note that factors that only involves MAX variables have already replicated during WMB construction
			kid.exact_val = wmbUB.heuristicTheta(kid.X, config);
			// actually, we can speed up this evaluation by computing each factor only once
			// for a leaf node, leaf_ub, leaf_lb are its ub, lb
			//			kid.down_ub = wmbUB.heuristicIn(kid.X, config);
			// Caveat: heuristics are differently calculated for MAX and SUM nodes
			kid.down_ub = wmbUB.heuristicInAug(kid.X, config, _K, kid.varType == MAX_VAR);
			// no change, same for both MAX and SUM variables
			kid.ub = dUB[v] = kid.down_ub + kid.exact_val;
			// whether solved
			//			kid.solved = checkSolved(kid);
			// those children actually have the same "solvedness" due to the way how "checkSolved" works
			// that's why even for the OR MAX case we do not need to track the largest child here
			checkSolved(kid);
			if (!kid.solved)
				ptrChild->solved = false;
			if (!kid.sumSolved)
				ptrChild->sumSolved = false;
			// delete a solved node only when all its siblings are solved, i.e., its parent is solved
			// it is serve as a book-keeping purpose, essential for two-step sampling
			// append
			exploredTree.append_child(ptrChild, kid);
			++GSIZE;
			if (kid.ub >= largestUpperBound) {
				largestUpperBound = kid.ub;
				largestKid = kid;
			}
		}
		if (!ptrChild->solved)
			ptr->solved = false;
		if (!ptrChild->sumSolved)
			ptr->sumSolved = false;
		ptrChild->ub = dUB.logsumexp();
/*		// Caveat for the case: AND-MAX parent + OR-SUM child: raised to the K-th power
		if (ptr->varType != ptrChild->varType) {
			cur_ub += _K * ptrChild->ub;
		} else {
			cur_ub += ptrChild->ub;
		}*/
		// for OR-SUM nodes, if "sumSolved", it must be "solved", remove all its children
		// while for OR-MAX nodes, "sumSolved" does not ensure its "solvedness"
		// remove a node's all children only when it's marked "solved".
		if (ptrChild->solved) {
			// debug
			if (DEBUG) {
				std::cout << "the following newly generated OR node's children will be removed" <<std::endl;
				printNode(*ptr);
			}
			//
			assert(ptrChild->sumSolved);
			// for a SUM node, we can delete all its children
			GSIZE -= exploredTree.number_of_children(ptrChild);
			exploredTree.erase_children(ptrChild);
			if (ptrChild->varType == MAX_VAR) {
				// for a MAX node, we have to keep the best configuration of MAP descendants
				ptrChild->mapDesc.push_back( std::make_pair(largestKid.id, largestKid.val) );
				// we will update "mapSpaceLeft" for its ancestors later on in "propMapSpace"
				ptr->mapSpaceLeft -= ptrChild->mapSpaceLeft; //TODO: remove it because no need to do it here?
				ptrChild->mapSpaceLeft = 0.0; // set to 0.0
				ptrChild->ub = largestKid.ub; // set to its largest child!!!
				if (DEBUG) {
					std::cout << "a just-generated OR-MAX node has been solved!" << std::endl;
				}
			}
		}

		// Caveat for the case: AND-MAX parent + OR-SUM child: raised to the K-th power
		if (ptr->varType != ptrChild->varType) {
			cur_ub += _K * ptrChild->ub;
		} else {
			cur_ub += ptrChild->ub;
		}

		// find best child
		findChild(ptrChild);
	}
	// more stable
	isProp = (ptr->ub > cur_ub + epsilon) || (ptr == root) || (ptr->solved);
	// propagate bounds if necessary
	if (isProp) {
//		ptr->ub = cur_ub; //Cause a bug: do not do it here, do it in updateBoundsPlusPruneSolved
		//		updateBounds(ptr);
		ptr = updateBoundsPlusPruneSolved(ptr); // do both simultaneously, be cautious!!!
	}
	if (verbose > 2 || DEBUG) {
		std::cout << "the best frontier node after expansion" << std::endl;
		printNode(*ptr);
	}
	// propagate the MAP space reduction info
	propMapSpace(ptr); // TODO: update when necessary?
	//
	return ptr;
}
int mmapIS::backwardUpdate(tree<node>::iterator ptr) {
/*
 * backward update after expanding the best
 */
	if (verbose > 2 || DEBUG) {
		std::cout << "running backwardUpdate..." << std::endl;
	}

	assert(ptr != NULL);
/*	if (ptr == NULL) {
		std::cout << "backwardUpdate: NULL occurs, problematic or in purpose!" << std::endl;
		return -1;
	}*/
	// re-calculate the best child since it may point to a solved node
	// check whether the whole tree has been solved or not
//	auto rt = exploredTree.begin();
	if ( root->solved ) {
		std::cout << "The root is solved!\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		// update the best config
		for (auto& pair : root->mapDesc) {
			mapConfig[pair.first] = pair.second;
		}
		return SOLVED;
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
tree<mmapIS::node>::iterator mmapIS::DFS(tree<node>::iterator ptr) {
	/*
	 * Depth-first search: expand one node at a time
	 * Implementing DFS without a stack to be compatible with pruning
	 * input: a frontier node to expand and its associated configuration
	 * output: a frontier node to expand next time
	 * ptr = NULL means no frontier node has been assigned, find the best to expand
	 */
	if (verbose > 2 || DEBUG) {
		std::cout << "Running DFS..." << std::endl;
	}
	/* ptr == NULL means no frontier node has been assigned,
	 * find the best frontier node to expand
	 */
	if (ptr == NULL) {
		ptr = root;  // the root is AND & MAX node
		// forward pass
		while (ptr->toChild != NULL) {
			// note that the virtual root is also AND_NODE with id = -1
			// first node in the loop would be OR_NODE node
			ptr = ptr->toChild;
			if (ptr->nodeType == AND_NODE) {
				dfsConfig[ptr->X] = ptr->val;
			}
		}
		if (verbose > 2 || DEBUG) {
			std::cout << "Best node found!" << std::endl;
			printNode(*ptr);
		}
		if (ptr->solved) {
			std::cout << "The best frontier node has been solved, quit!" << std::endl;
//			return 3;
			return NULL;
		} else {
			if (DEBUG && ptr->sumSolved)
				std::cout << "The SUM problem below of the best frontier node has been solved." << std::endl;
		}
		assert(ptr->nodeType == AND_NODE); // without removal, the frontier must be AND
	}

	if (DEBUG) {
		std::cout << "Best frontier node for DFS to expand:\n";
		printNode(*ptr);
	}

	/*******************************************/
	ptr = expandOneFrontierNode(ptr, dfsConfig); // get the highest unsolved node
	int msg = backwardUpdate(ptr);
	if (msg == SOLVED) {
		return NULL;
	}
	/* get the next to expand
	 * no need to do "dfsConfig[ptr->X] = ptr->val" even if ptr is an AND node for the first node
	 * because it must be an ancestor of the input of "expandOneFrontierNode"
	 */
	while(ptr->toChild != NULL) {
		ptr = ptr->toChild;
		if (ptr->nodeType == AND_NODE) {
			dfsConfig[ptr->X] = ptr->val;
		}
	}
	assert(ptr->nodeType == AND_NODE);
	return ptr; // next node for DFS to expand
}
/*void mmapIS::runSearch(const unsigned treeSizeLimit) {
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
int mmapIS::runSearch(const unsigned nnd, const unsigned long treeSizeLimit) {
/*
 * nnd: no. of nodes to expand
 *
 */
//	int base = GSIZE;
/*	double power = 1.0 + std::floor(std::log2(GSIZE+1)); // add

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
	return 0;*/

	if (verbose >1) {
		std::cout << "Run " << nnd << " rounds of best-first search...\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	for(unsigned cnt=0; cnt<nnd; ++cnt){
		/*****************************/
		expandBest();
		/*****************************/
		if ( root->solved ) {
			display();
			return SOLVED;
		}
		if ( GSIZE >= treeSizeLimit ) {
			std::cout << "Reach tree size limit: " << treeSizeLimit << std::endl;
			display();
			return MEMOUT;
		}
		if ( mex::timeSystem()-startProcess > timeBudget ) {
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			display();
			return TIMEOUT;
		}
		if ( _nExpansion >= _increment ) {
			_increment += _increment; // to avoid expensive "log"
			display();
		}
	}
	if (verbose > 1 || DEBUG) {
		display();
	}
	return 0;
}
int mmapIS::runDFS(const unsigned nnd) {
/*
 * expand a given number of frontier nodes via DFS
 * nnd: no. of nodes to expand
 */
	if (verbose >1) {
		std::cout << "Run " << nnd << " rounds of depth-first search...\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	}
	/*****************************************/
	for(unsigned cnt=0; cnt<nnd; ++cnt){
		dfsFrontier = DFS(dfsFrontier); // DFS expands one node at a time
		// debug
		if (dfsFrontier == NULL) {
			assert(root->solved);
		}
		/***********************/
		if ( root->solved ) {
			display();
			return SOLVED;
		}
		if ( mex::timeSystem()-startProcess > timeBudget ) {
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			display();
			return TIMEOUT;
		}
		if ( _nExpansion >= _increment ) {
			_increment += _increment; // to avoid expensive "log"
			display();
		}
	}
	/***********************/
	if (verbose > 1 || DEBUG) {
		display();
	}
	return 0;
}
double mmapIS::logF(const mex::vector<uint32_t>& config, const std::list<mex::Var>& Done) {
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
tree<mmapIS::node>::iterator mmapIS::upperSampling(tree<node>::iterator ptr){
/*
 * sample for an OR node based on current upper bounds of its AND children;
 * for each AND child, probability to be picked is proportional to its current upper bound.
 */
	if (DEBUG) {
		std::cout << "running upperSampling..." << std::endl;
//		printNode(*ptr);
	}
	assert( (*ptr).nodeType == OR_NODE );
	auto X = (*ptr).X;
	auto ns = X.states();
	if (DEBUG) {
		std::cout << "ns: " << ns << ", nchildren: " << exploredTree.number_of_children(ptr) << std::endl;
		if (ns!=exploredTree.number_of_children(ptr)) {
			std::cout << "printing its children:\n";
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) printNode(*sib);
		}
	}
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
	 // should not run into the following lines
	std::runtime_error("something wrong in upperSampling!");
	return NULL;
}
std::list<mex::Var> mmapIS::getDescList(int id) {
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
void mmapIS::getDesc(int ID, std::list<mex::Var>& unDoneMax, std::list<mex::Var>& unDoneSum, bool maxDescOnly) {
/*
 * get a list of MAX descendants, and a list of SUM descendants,
 * those are ordered lists.
 * maxDescOnly: whether we only need the MAX descendant, default: false
 *
 * TODO: faster implementation?
 */
	unDoneMax.clear(); unDoneSum.clear();
	std::stack<int> stack;
	const auto& childrenList = (ID >= 0)? ascList[ID] : rts;
//	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
//		stack.push(wmbUB.var(*it));
//	}
	for (int k : childrenList) stack.push(k);
//
	while(!stack.empty()) {
		int id = stack.top(); stack.pop();
		auto X = wmbUB.var(id);
		if (_varTypeVec[id] == MAX_VAR) {
			unDoneMax.push_back(X);
		} else {
			if (maxDescOnly) {
				// if we only need those MAX descendants
				continue;
			}
			unDoneSum.push_back(X);
		}
		for (int k : ascList[id]) stack.push(k);
	}
}

std::pair<int,int> mmapIS::getDesc(int ID) {
/*
 * get location of ID's descendants in "__pseudotree"
 * indices from "pair.first" to "pair.second - 1" are its descendants including itself
 */
	if (DEBUG) {
		std::cout << "running getDesc..." << std::endl;
	}
	if (ID < 0) {
		return std::make_pair(-1, nvar);
	}
	return _locPair[ID];
}

/*void mmapIS::addSample(double logFx, const std::vector<double>& logQxCond) {

 * add one sample to the root

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
//	if ( (*ptr).nodeType == OR_NODE ) return;
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
}*/
void mmapIS::addSampleToRoot(double wt) {
/*
 * add one sample to update the overall estimate
 * run after drawing a full sample
 * wt: importance weight for Z_aug
 */
	++NSP;
	// update necessary quantities
//	updateConvEstimate(wt, root->ub);
	updateNormalizedEstimate(wt, root->ub);
	// update accordingly
	sumWeightedMapSpaceSize  = logsumexp({sumWeightedMapSpaceSize, root->mapSpaceLeft-root->ub});
	// debug
	if (DEBUG) {
		std::cout << "after adding one sample: \n";
		std::cout << "NSP: " << NSP << "\n";
		std::cout << "normalizedEx: " << normalizedEx << ", normalizedEx2: " << normalizedEx2 << "\n";
		std::cout << "sumInvUB: " << sumInvUB << ", hmUB: " << hmUB
				<< " sumWeightedMapSpaceSize: " << sumWeightedMapSpaceSize <<  std::endl;
	}
}
double mmapIS::addSampleToNode(tree<node>::iterator ptr, const mex::vector<uint32_t>& config, double est) {
/*
 * Add an estimate created from (1 or K) samples to a node
 * Used for two type of nodes:
 *  1) frontier AND-MAX node (no descendants in the current search tree)
 *  2) OR-SUM node with a MAX parent
 *  est: for MAX nodes, it is an estimate of the underneath MAP problem (not the SUM problem of the augmented model)
 *  	 for SUM nodes, it is an estimate of the underneath SUM problem, gathered from _K samples.
 */
	if (DEBUG) {
		std::cout << "running addSampleToNode...\n";
	}
	assert(!ptr->solved); // should not add to the solved node
	if (ptr->varType == MAX_VAR) {
		assert(exploredTree.begin(ptr) == exploredTree.end(ptr) && ptr->nodeType == AND_NODE); // must be a frontier AND node
		if (ptr->mapDesc.empty()) {
			// if we have not set a sub-config for this MAX node yet, set it now
		/*	std::list<mex::Var> unDoneMax, unDoneSum;
			getDesc(ptr->id, unDoneMax, unDoneSum, true);
			for (auto& X : unDoneMax) {
				auto id  = X.label();
				ptr->mapDesc.push_back( std::make_pair(id, config[id]));
			}*/
			addMapDesc(ptr, config); // to replace the above
			ptr->downEst = est;
			ptr->nsp = 1;
		} else {
			// for MAX nodes, "nsp" is the no. of its current best MAP configuration being sampled
			// "nsp" is useful for those frontier MAX nodes
			assert(ptr->nsp > 0);
			// check whether these two sub-configurations are the same; if so, merge; otherwise, pick the larger one.
			bool isIdentical = true;
			for (const auto& pair : ptr->mapDesc) {
				if (pair.second != config[pair.first]) {
					isIdentical = false;
					break;
				}
			}
			//
			if (isIdentical) {
				// merge
				ptr->downEst = logsumexp( {est,  ptr->downEst + log(ptr->nsp)} ) - log(ptr->nsp+1);
				ptr->nsp += 1;
			} else {
				// compare, pick the larger one.
				// Caveat: a bit unfair if the stored configuration has been sampled multiple times
				// TODO: any better way?
				if (est > ptr->downEst + epsilon) {
					ptr->downEst = est;
					ptr->nsp = 1;
					// reset
					for (auto& pair : ptr->mapDesc) {
						pair.second = config[pair.first];
					}
				}
			}
		}
	} else {
		// OR-SUM node with a MAX parent
		assert(exploredTree.parent(ptr)->varType == MAX_VAR);
		if (ptr->nsp > 0) {
			auto sumInvUB = ptr->wSum - ptr->downEst; // \sum_i 1/U_i = N/HM(U)
			sumInvUB = logsumexp( {sumInvUB, _logK - ptr->ub} );
			ptr->wSum = logsumexp( {ptr->wSum, est + _logK - ptr->ub} ); // est averages K importance weights
			ptr->downEst = ptr->wSum - sumInvUB;
			ptr->nsp += _K;
			// debug
			// HM(U) >= downEst, sumInvUB = N/HM(U)
			assert(log(ptr->nsp)-sumInvUB + epsilon > ptr->downEst);
		} else {
			ptr->wSum = est + _logK - ptr->ub; // est averages K importance weights
			ptr->downEst = est;
			ptr->nsp = _K;
			// debug
			assert(ptr->ub + epsilon > est);
		}
	}
	return ptr->downEst;
}
void mmapIS::addMapDesc(tree<node>::iterator ptr, const mex::vector<uint32_t>& config) {
/*
 * add Max descendants' configuration to "ptr->mapDesc"
 */
	if (DEBUG) {
		std::cout << "running addMapDesc...\n";
	}
//	assert(ptr->mapDesc.empty() && ptr->varType==MAX_VAR);
	const auto& childrenList = (ptr->id > -1) ? ascList[ptr->id] : rts;
	std::stack<uint32_t> stack;
	for (auto id : childrenList) {
		if (_varTypeVec[id] == MAX_VAR) stack.push(id);
	}
	while (!stack.empty()) {
		auto ID = stack.top(); stack.pop();
		ptr->mapDesc.push_back( std::make_pair(ID, config[ID]) );
		for (auto id : ascList[ID]) {
			if (_varTypeVec[id] == MAX_VAR) stack.push(id);
		}
	}
}
/*std::pair<double, double> mmapIS::calcRootEBB(){
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
std::pair<double, double> mmapIS::calcRootEBB(double upperbound){

	 * calculate EBB for root, using given upper bound

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
std::pair<double, double> mmapIS::calcHoeffding() {

 * Hoeffding's bounds on the convex combination of samples
 * convex weights are proportional to 1/U_i^2

	double bias = ( log( -log(DELTA)/2.0 ) - sumInvSqrUB )/2.0; //
	double lb = -std::numeric_limits<double>::infinity();
//	double ub = log(exp(convEx) + rng);
	double ub = logsumexp( {convEx, bias} );
	if ( convEx > bias) {
//			lb = log( exp( convEx ) - exp( bias ) );
		lb = convEx + log(1 - exp( bias-convEx ) );
	}
	return std::pair<double,double>(lb,ub);
}*/
std::pair<double, double> mmapIS::calcNormalizedEBB() {
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
bool mmapIS::calcMapEBB(double& ubZaug, double& lbZaug, double& lbMap) {
/*
 * an extension of "calcNormalizedEBB"
 * 1. EBB for Z_aug based on the normalized estimate: HM(U)*(\sum_i z_i/u_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
 * 2. probabilistic lower bound of the optimal MAP solution
 * If EBB does not kick in, use Markov lower bounds instead: Zhat*delta
 */
	bool isMarkov = false;
	lbZaug = lbMap = -std::numeric_limits<double>::infinity();
	ubZaug = std::numeric_limits<double>::infinity();
	if ( NSP <= 1 ) {
		return isMarkov;
	}

	double zhat = hmUB + normalizedEx; // unbiased estimate of Z
	double var = std::max(exp(normalizedEx2) - exp(normalizedEx)*exp(normalizedEx), 0.0); // to be numerically stable
	var *= NSP/(NSP-1); // should be the unbiased sample variance
	double bias = log( sqrt(2.0*var*log(2.0/DELTA)/NSP) + 7.0*log(2.0/DELTA)/3.0/(NSP-1) ) + hmUB; // store in log
	ubZaug = logsumexp( {bias, zhat} );
	if ( zhat > bias) {
		lbZaug = zhat + log(1 - exp( bias-zhat ) );
		lbMap = ( lbZaug + log(NSP) - hmUB - sumWeightedMapSpaceSize ) / _K;
	} else {
		// Markov lower bound
		isMarkov = true;
		lbZaug = zhat + log(DELTA);
		lbMap = ( lbZaug + log(NSP) - hmUB - sumWeightedMapSpaceSize ) / _K;
	}
	if (DEBUG){
		std::cout << "zhat: " << zhat << " bias: " << bias << " ubZaug: " << ubZaug << std::endl;
		std::cout << "hmUB: " << hmUB << " normalizedEx: " << normalizedEx << std::endl;
	}
	return isMarkov;
}

/*void mmapIS::updateConvEstimate(double wt, double ub) {

 * update convEx, sumSqrUB, sumInvSqrUB
 * wt: importance weight, i.e., current sample (in log)
 * ub: corresponding global upper bound on the sample (in log)
 * Caveat: NSP should've already been updated

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
}*/
void mmapIS::updateNormalizedEstimate(double wt, double ub) {
/*
 * update the normalized estimate: HM(U)*(\sum_i Z_i/U_i)/n, U=(u_1, ..., u_n), HM(U) is harmonic mean of U
 * empirical Bernstein bound is applicable to this estimate
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
/*void mmapIS::twoStepSampling() {

 * two-stage sampling with a fixed tree

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
				if (ptr->nodeType == AND_NODE) {
					// For solved AND node, the corresponding variable is IN config.
					solved_val += ptr->ub - ptr->exact_val; // exact_val will be computed logF(config, Done);
				} else {
					// For solved OR node, the corresponding variable is NOT in config.
					solved_val += ptr->ub;
				}
			} else {
				// an unsolved frontier node should be AND node due to the way we expand the best
				assert(ptr->nodeType == AND_NODE);
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
			if (ptr->nodeType == OR_NODE) {
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
}*/
void mmapIS::twoStepSampling() {
	/*
	 * two-step sampling for the augmented model: each SUM variable has to be sampled K times
	 * Two tasks we accomplish:
	 * 	1). compute the value that contributes to the estimate of Z_aug
	 * 	2). compute the estimate of the sampled MAP solution and update info at each sampled tip MAP node
	 */
	if (DEBUG) {
		std::cout << "running twoStepSampling...\n";
	}
	/***************** sample a full solution tree ******************/
	auto ptr = root;
	//values for solved variables, summation part of Rao-Blackwellisation
//	double solved_val = 0.0;
/*	std::list<mex::Var> unDoneMax;
	std::list<mex::Var> unDoneSum;*/
	double maxConfigEst = 0.0; // estimated value of the sampled MAX configuration
	double augEst = 0.0; // estimated value of Z_aug from this sample
	std::vector<double> estVec(_K, 0.0);
	double logFx, logQx, aveSumEst, prodSumEst;
	// step 1 and 2
	std::stack<tree<node>::iterator> stack;
	stack.push(root);
	while (!stack.empty()) {
		ptr = stack.top();
		stack.pop();
		// debug
		if (DEBUG) {
			std::cout << "node popped out during twoStepSampling:" << std::endl;
			printNode(*ptr);
		}
		/* if we see a SUM node, it must be a child of some MAX node
		 * draw K samples, and add information to that node
		 */
		if (ptr->varType == SUM_VAR) {
			assert(ptr->nodeType == OR_NODE); // because we start from the root which is an AND-MAX node
			if (ptr->solved) {
				// Caveat: it's possible this SUM node is solved
				augEst += _K * ptr->ub;
				maxConfigEst += ptr->ub;
				// debug
				if (DEBUG) {
					std::cout << "solved OR-SUM node in sampling:\n";
				}
			} else {
				for (int j = 0; j < _K; ++j) {
					estVec[j] = twoStepSum(tuple, ptr);
					augEst += estVec[j];
					if (DEBUG) {
						std::cout << "est[j]: " << estVec[j] << "\n";
					}
				}
				auto tmpval = addSampleToNode(ptr, tuple, logsumexp(estVec) - _logK);
				maxConfigEst += tmpval;
				// no need to run the rest
				if (DEBUG) {
					std::cout << "tmpval: " << tmpval << "\n";
				}
			}
			// no need to add its children to the stack
			continue;
		}
		/*
		 * Now we only need to deal with MAX nodes in the sequel
		 */
		if (ptr->nodeType == AND_NODE) {
			maxConfigEst += ptr->exact_val / _K; // since the K-th power is used in creating factors for MAX variables
			augEst += ptr->exact_val;
		}

		/**************************************************************/
		if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
		/*
		 * frontier node, any solved node should be a frontier node due to pruning
		 * 		if "sumSolved", sufficient to just sample MAX descendants
		 *	 	if "solved", no sampling necessary
		*/
			if (ptr->sumSolved) {
				// the summation problem below has been solved but the maximization problem might not
				augEst += ptr->ub; // no matter whether it is MAX solved or not
				if (ptr->nodeType == AND_NODE)
					augEst -= ptr->exact_val; // already counted previously
				//
				if (ptr->solved) {
					// fully solved. no sampling, copy the MAX configuration
					for (auto& pair : ptr->mapDesc)
						tuple[pair.first] = pair.second; //  <id, val> pairs
					// no matter ptr is AND or OR
					maxConfigEst += ptr->ub / _K;
					if (ptr->nodeType == AND_NODE)
						maxConfigEst -= ptr->exact_val / _K; // already counted previously
				} else {
				/* only the summation has been solved. sample all descendants
				 * sample the SUM descendants only once instead of K times because of no variance
				 * this node should've never been expanded, "sumSolved" due to "exactHeur", thus safe to sample once only
				 * TODO: not sample the SUM variables, take heuristics to accelerate
				*/
/*					getDesc(ptr->id, unDoneMax, unDoneSum);
					// sample once for the SUM variables to get the exact value of the sum subproblem since any sample gives the same value
					wmbUB.conditionalMixtureSampling(ptr->id, tuple, 1, unDoneMax, unDoneSum, logFx, logQx, aveSumEst, prodSumEst);*/

					// a faster version
					auto locPair = getDesc(ptr->id); // ptr must be AND node
					wmbUB.conditionalMixtureSampling(ptr->id, tuple, 1, locPair.first, locPair.second,
							_pseudotree, _varTypeVec, MAX_VAR, logFx, logQx, aveSumEst, prodSumEst);

//					maxConfigEst += logFx/_K + aveSumEst;
					maxConfigEst += addSampleToNode(ptr, tuple, logFx/_K + aveSumEst);
				}
			} else {
				// an unsolved (not even sumSolved) frontier MAX node should be an AND node due to the way we expand nodes
				assert(ptr->nodeType == AND_NODE);
				// sample from the frontier node according to the mixture proposal
				// variables for sampling
				//				auto unDone = getDescList(ID);
				//				std::list<mex::Var> unDone, tipVars;
				//				getDescVars(ID, unDone, tipVars);
				// sample those variables in "unDone" conditioned on the path
				// config, logQxCond will be updated
				//				wmbUB.conditionalMixtureSampling(ID, config, logQxCond, unDone);
/*				double wt = wmbUB.conditionalMixtureSampling(ID, tuple, rts, ascList, _depthVec, _adaMaxDepthVec[ID]);
				// add this sample to the corresponding frontier node
				double fx = logF(config, unDone, tipVars);
				double qx = 0.0;
				for (const auto& X : unDone) {
					qx += logQxCond[X.label()];
				}
				//				addSampleToNode(*ptr, fx-qx);
				addSampleToNode(*ptr, wt);
				if (DEBUG) {
					printNode(*ptr);
				}*/
				// add a set of sampled variables
				//				Done.insert(Done.end(), unDone.begin(), unDone.end());

/*				getDesc(ptr->id, unDoneMax, unDoneSum);
				wmbUB.conditionalMixtureSampling(ptr->id, tuple, _K, unDoneMax, unDoneSum, logFx, logQx, aveSumEst, prodSumEst);*/

				// a faster version
				auto locPair = getDesc(ptr->id);
				wmbUB.conditionalMixtureSampling(ptr->id, tuple, _K , locPair.first, locPair.second,
						_pseudotree, _varTypeVec, MAX_VAR, logFx, logQx, aveSumEst, prodSumEst);

				augEst += logFx - logQx + prodSumEst; // sum up importance weights
				maxConfigEst += addSampleToNode(ptr, tuple, logFx/_K + aveSumEst);
			}
		} else {
			// internal node
			assert(!ptr->solved); // debug, internal nodes should not be solved due to pruning
			if (ptr->nodeType == OR_NODE) {
				augEst += ptr->ub; // see below of explanation
				ptr = upperSampling(ptr); //return the sampled AND child
				/* so, we actually compute  "augEst -= ptr->ub - ptrPar->ub" using two steps
				 * where "ptr->ub - ptrPar->ub" is the probability of this AND node being sampled
				 */
				augEst -= ptr->ub;
				tuple[ptr->X] = ptr->val;
				stack.push(ptr);
			} else {
				// push all OR children to stack
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					stack.push(sib);
				}
			}
		}
	}

	/************update the estimate of Z_aug*****************/
	if (DEBUG) {
		std::cout << "augEst/K: " << augEst/_K << ", GSIZE: " << GSIZE << std::endl;
	}
	assert(augEst < root->ub + epsilon); // check boundedness
	addSampleToRoot(augEst);

	/**************update the current best MAP configuration***************/
	/* its estimated value might be changed due to sampling or node expansion
	 * so we have to re-calculate it before comparing it to the current sampled MAP configuration
	 * Caveat: it's possible that some expansions may result in loss of "downEst" and "mapDesc". TODO: Fix it?
	 * TODO: faster by only updating those nodes being just sampled?
	 */
	if (isMapConfigNeverSampled) {
		// for the first time
		isMapConfigNeverSampled = false;
		mapConfig = tuple;
		mapSolEst =  maxConfigEst;
		//
		if (DEBUG){
			std::cout << "current best MAP configuration:\n";
			for (int i=0; i<nvar; ++i) {
				if (_varTypeVec[i]==MAX_VAR) std::cout << "(" << i << ", " << mapConfig[i] << "), ";
			}
			std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		}
		//
		return;
	}
	double mapSolEstBeforeReCalc = mapSolEst;
	if (DEBUG) {
		std::cout << "maxConfigEst: " << maxConfigEst << std::endl;
		std::cout << "MAP solution estimated value:\n";
		std::cout << "before update: " << mapSolEst << "\n";
	}
	mapSolEst = 0.0; // re-set before update
	stack.push(root);
	while (!stack.empty()) {
		ptr = stack.top();
		stack.pop();
		// debug
		if (DEBUG) {
			std::cout << "node popped out during mapSolEst re-calculation:" << std::endl;
			printNode(*ptr);
		}
		//
		if (ptr->varType == SUM_VAR) {
			// OR-SUM node with a MAX parent
//			mapEstVec[ptr->id] = ptr->downEst;
			// debug
			assert(ptr->nodeType == OR_NODE);
			mapSolEst += ptr->downEst;
			continue;
		}
		// Now the case only for MAX nodes
		if (ptr->nodeType == AND_NODE) {
			// AND-MAX nodes
//			mapEstVec[ptr->id] = ptr->exact_val / _K;
			mapSolEst += ptr->exact_val / _K;
			if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
				// a frontier MAX node, may or may not be solved
				if (ptr->solved) {
//					mapEstVec[ptr->id] += (ptr->ub - ptr->exact_val) / _K;
					mapSolEst += (ptr->ub - ptr->exact_val) / _K;
				} else {
//					mapEstVec[ptr->id] += ptr->downEst;
					mapSolEst += ptr->downEst;
				}
				// update anyway since ptr->mapDesc possibly changed
				for (auto& pair : ptr->mapDesc) {
					mapConfig[pair.first] = pair.second;
				}
			}
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) stack.push(sib);
		} else {
			// OR-MAX nodes
			if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
				assert(ptr->solved); // A frontier OR node must be solved!
				mapSolEst += ptr->ub / _K; //
				// update anyway since ptr->mapDesc possibly changed
				for (auto& pair : ptr->mapDesc) {
					mapConfig[pair.first] = pair.second;
				}
			}
			// if internal OR-MAX nodes, add its corresponding child to the stack
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if (mapConfig[sib->id] == sib->val) {
					stack.push(sib);
					break;
				}
			}
		}
	}
	if (DEBUG) {
		std::cout << "after re-calculation: " << mapSolEst << "\n";
		if ( std::isinf(mapSolEst) ) {
			std::cout << "mapSolEst becomes -inf after re-calculation!!!" << std::endl;
		}
	}
	// TODO: better to fix systematically ?
	if ( std::isinf(mapSolEst) ) {
		// mapSolEst = -inf happens due to information loss in node expansion: descendants may not have "down_est" assigned
		mapSolEst = mapSolEstBeforeReCalc; // this is just an ad-hoc solution
	}
	// compare with the just-sampled MAP configuration
	// if the sampled one is better, replace
	if (mapSolEst < maxConfigEst) {
		mapConfig = tuple;
		mapSolEst =  maxConfigEst;
	}
	//
	if (DEBUG) {
		std::cout << "after update: " << mapSolEst << "\n";
		std::cout << "current best MAP configuration:\n";
		for (int i=0; i<nvar; ++i) {
			if (_varTypeVec[i]==MAX_VAR) std::cout << "(" << i << ", " << mapConfig[i] << "), ";
		}
		std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
}
double mmapIS::twoStepSum(mex::vector<uint32_t>& config, tree<node>::iterator ptr) {
/*
* two-step sampling for SUM variables by conditioning on their MAP ancestors
* return the accumulated value
*/
	if (DEBUG) {
		std::cout << "running twoStepSum..." << std::endl;
	}
	// debug
	auto ptr_orig = ptr;
	// must be a SUM node whose parent is a MAX node
	assert(ptr->nodeType == OR_NODE && ptr->varType == SUM_VAR);
	auto ptrPar = exploredTree.parent(ptr);
	assert(ptrPar->varType == MAX_VAR);
	// if solved, do nothing
	if (ptr->solved) return ptr->ub;
	//
	double logFx, logQx, aveSumEst, prodSumEst;
	std::stack<tree<node>::iterator> stack;
	stack.push(ptr);
//	double solved_val = 0.0; // including exact_val
	double local_est = 0.0; // random estimate of the SUM problem rooted at this node
	while (!stack.empty()) {
		ptr = stack.top();
		stack.pop();
		// including exact_val
		if (ptr->nodeType == AND_NODE) {
//			solved_val += ptr->exact_val;
			local_est += ptr->exact_val;
		}
		if (exploredTree.begin(ptr) == exploredTree.end(ptr)) {
			// if this is a frontier node
			if (ptr->solved) {
				//any solved node should be a frontier node; no sampling for a solved node.
				if (ptr->nodeType == AND_NODE) {
					// For solved AND node, the corresponding variable is IN config.
//					solved_val += ptr->ub - ptr->exact_val; // exact_val already added
					local_est += ptr->ub - ptr->exact_val; // "exact_val" already counted
				} else {
					// For solved OR node, the corresponding variable is NOT in config.
//					solved_val += ptr->ub;
					local_est += ptr->ub;
				}
			} else {
				// an unsolved frontier node should be AND node due to the way we expand the best
				assert(ptr->nodeType == AND_NODE);
				// sample from the frontier node according to the mixture proposal

/*				// TODO: be more efficient by directly fetching from memory?
				std::list<mex::Var> unDoneMax;
				std::list<mex::Var> unDoneSum;
				getDesc(ptr->id, unDoneMax, unDoneSum);
				//debug, "unDoneMax" must be empty for a SUM node
				assert(unDoneMax.empty());
				// K = 1 here
				wmbUB.conditionalMixtureSampling(ptr->id, config, 1, unDoneMax, unDoneSum, logFx, logQx, aveSumEst, prodSumEst);*/

				// a faster version
				auto locPair = getDesc(ptr->id);
				wmbUB.conditionalMixtureSampling(ptr->id, config, 1, locPair.first, locPair.second,
						_pseudotree, _varTypeVec, MAX_VAR, logFx, logQx, aveSumEst, prodSumEst);

				//debug, aveSumEst == prodSumEst for K=1
				assert(aveSumEst < prodSumEst + epsilon  && aveSumEst > prodSumEst - epsilon );

				local_est += prodSumEst;
			}
		} else {
			// an internal node
			assert(!ptr->solved); // debug, internal SUM nodes should not be solved due to pruning
			if (ptr->nodeType == OR_NODE) {
				ptrPar = ptr;
				ptr = upperSampling(ptr); // return an AND child
				config[ptr->X] = ptr->val;
				stack.push(ptr);
//				logQxCond[ptr->id] = ptr->ub - ptrPar->ub; // ub of an OR parent is equal to sum of all AND children's ub
				// add a sampled variable
//				Done.push_back(ptr->X);

				// "ptr->ub - ptrPar->ub" is the probability
				local_est -= ptr->ub - ptrPar->ub; // ub of an OR parent is equal to sum of all AND children's ub
			} else {
				// push all OR children to stack
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					stack.push(sib);
				}
			}
		}
	}
	if (DEBUG) {
		std::cout << "local_est: " << local_est << " < local_ub: " << ptr_orig->ub << std::endl;
	}
	assert(local_est < ptr_orig->ub + epsilon);
	return local_est;
}
int mmapIS::doSampling(int nsp,  bool isInitUB) {
/*
 * draw nsp samples based on current search tree
 */
	if (DEBUG) {
		std::cout << "running doSampling...\n";
	}
//	std::pair<double, double> prob;
	int cnt = 0;
//	double power = 1.0 + std::floor(std::log2(NSP+1)); // add 1 to NSP to avoid log(0)
	while(cnt < nsp) {
		++cnt;
		twoStepSampling();
//		outNow = verbose > 1;
//		if (NSP >= outFrequency) {
//			outNow = true;
//			outFrequency += outFrequency;
//		}
//		if ( (NSP == long(pow(2.0, power)))  || (verbose > 1) ) {
		if ( NSP == _outFrequency ) {
			_outFrequency += _outFrequency;
//			power += 1.0;
//			prob = calcRootEBB();
/*			if ( isInitUB ) prob = calcRootEBB(_initUB);
			else prob = calcRootEBB();
			std::cout.precision(10);
//			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
			std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
//			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
			**** more ways to aggregate estimates and bounds ****
			prob = calcHoeffding();
			std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
			prob = calcNormalizedEBB();
			std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
			*******************************
			std::cout << "Tree size: " << GSIZE <<"\n";
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;*/
			display();
		}

		if ( mex::timeSystem()-startProcess > timeBudget ) {
/*			if ( isInitUB ) prob = calcRootEBB(_initUB);
			else prob = calcRootEBB();
			std::cout.precision(10);
//			std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
			std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
//			std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
			std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
			**** more ways to aggregate estimates and bounds ****
			prob = calcHoeffding();
			std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
			prob = calcNormalizedEBB();
			std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
			*******************************
			std::cout << "Tree size: " << GSIZE <<"\n";*/

			display();
			std::cout << "Reach time limit (sec): "<< timeBudget  << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
			return TIMEOUT;
		}
	}
//	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//	if (NSP >= nSample) std::cout << "Reach sample size limit: " << nSample << std::endl;
	if (verbose > 1) {
/*//		prob = calcRootEBB();
//		prob = calcRootEBB(_initUB);
		if ( isInitUB ) prob = calcRootEBB(_initUB);
		else prob = calcRootEBB();
		std::cout.precision(10);
//		std::cout << "[" << mex::timeSystem() - startProcess << "]:  nsp/nSample = " << NSP << "/"  << nSample <<"\n" ;
		std::cout << "[" << mex::timeSystem() - startProcess << "]:  no. of samples: " << NSP <<"\n" ;
//		std::cout << "Deterministic bounds: "<< std::setw(10) << rt->lb << " < " << rt->logEx << " < " << rt->ub << "\n";
		std::cout << "Deterministic bounds: "<< std::setw(10) << " < " << root->logEx << " < " << root->ub << "\n";
		std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << root->logEx << " < " << prob.second << "\n";
		**** more ways to aggregate estimates and bounds ****
		prob = calcHoeffding();
		std::cout << "Hoeffding bounds: " << std::setw(10) << prob.first << " < " << convEx << " < " << prob.second << "\n";
		prob = calcNormalizedEBB();
		std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << hmUB + normalizedEx << " < " << prob.second << "\n";
		*******************************
		std::cout << "Tree size: " << GSIZE <<"\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;*/
		display();
	}
	return 0;
}
int mmapIS::display() {
/*
 * Display function
 * Output the normalized results, e.g., Zhat/_K instead of Zhat
 */
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	UB = root->ub;
//	auto prob = calcNormalizedEBB();
	double ub, lb, lbMap;
	bool isMarkov = calcMapEBB(ub, lb, lbMap);
	std::string str = isMarkov? " (Markov)" : "";
	auto Zhat = hmUB + normalizedEx;
	std::cout.precision(10);
	std::cout<<"["<<mex::timeSystem()-startProcess<<"]: Tree size: " << GSIZE << ", no. of expansions: " << _nExpansion
			<< ", no. of samples: " << NSP << ", MAP space size (log2): " << root->mapSpaceLeft / log(2) << "\n";
	std::cout << "Deterministic bounds: " << std::setw(10)<< LB/_K <<" < ln Z < "<< UB/_K << "\n";
//	std::cout << "Normalized EBB: " << std::setw(10) << prob.first << " < " << Zhat << " < " << prob.second;
	std::cout << "EBB" << str << " for Zaug: " << std::setw(10) << lb/_K << " < " << Zhat/_K << " < " << ub/_K  << "\n";
	//Note that "mapSolEst" may not be larger than "lbMap"
	std::cout << "EBB" << str << " for MAP: " << std::setw(10) << lbMap << " < " << mapSolEst << " < " << ub/_K << "\n";
//	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
	std::cout << "++++++++++++++++++++++++++Printing current MAP configuration++++++++++++++++++++++++++++\n";
	for (int i=0; i<nvar; ++i) {
		if ( _varTypeVec[i] == MAX_VAR ) {
			std::cout << "(" << i << "," << mapConfig[i] << "), ";
		}
	}
	std::cout << "\n++++++++++++++++++++++++++Done printing+++++++++++++++++++++++++++++\n"<<std::endl;
	return 0;
}

void mmapIS::start(const int nsp, const int nnd) {
/*
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
/*	if (DEBUG) {
		treeSizeLimit = 7;
	}*/
	std::cout << "Time limit (sec): " << timeBudget << std::endl;
	std::cout << "Memory limit (MB): " << memoryBudget << std::endl;
	std::cout << "Tree size limit: " << treeSizeLimit << std::endl;
	std::cout << "nsp: " << nsp << ", nnd: " << nnd << ", K: " << _K << std::endl;
	assert(nsp>0 || nnd>0);
	bool isMemOut = false;
	int status = 0;
	while (true) {
		// -1: memory out
		// -2: time out
/*		// run search first
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
		}*/
		if (!isMemOut) {
			status = runSearch(nnd, treeSizeLimit);
		} else {
			status = runDFS(nnd);
		}
		if (status == MEMOUT) {
			std::cout << "memory out, switch to DFS!" << std::endl;
//			// switch to DFS
			isMemOut = true;
//			break;
		}
		if (status == TIMEOUT) {
			std::cout << "time out, quit!" << std::endl;
			break;
		}
		if (status == SOLVED) {
			std::cout << "problem solved, quit!" << std::endl;
			break;
		}
		doSampling(nsp);
	}
}
// EOF

/*
 * mmap.cpp
 *
 *  Created on: Oct 28, 2016
 *      Author: qlou
 */

#include "mmap.h"

mmap::mmap(mex::wmbe& wmbub, int vbs, int nv, int nmax, double ub, double mb, double tb):
wmbUB(wmbub), verbose(vbs), nvar(nv), nMaxVar(nmax), UB(ub), memoryBudget(mb),
timeBudget(tb), ascList(nv), maxConfig(nv, 0), varTypeVec(nv, SUM_VAR), tuple(nv,-1), exactHeur(nv, true) {
	// TODO Auto-generated constructor stub
	startTime = mex::timeSystem();
	buildPseudoTree();
	setRootConfig();
	// initialize varTypeVec
	auto maxVars = wmbUB.getMaxVars();
	for (auto it = maxVars.begin(); it != maxVars.end(); ++it) {
		varTypeVec[ (*it).label() ] = MAX_VAR;
	}
	setExactHeur();
}

mmap::~mmap() {
	// TODO Auto-generated destructor stub
	std::cout<<" ==== Done best-first search for MMAP ====" << std::endl;
}

void mmap::buildPseudoTree() {
	mex::vector<mex::graphModel::vindex> parents = wmbUB.getPseudotree();
//	assert(ascList.size()==nvar); // debug
//	assert((ascList[0]).empty()); // debug
//	assert((ascList[nvar-1]).empty()); // debug
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

void mmap::setRootConfig() {
/*
 * setup for root
 */
	// this is a virtual root node, nodes without parents are children of this virtual root node
	node n; // virtual root is both MAX_NODE and AND_NODE, which is an crucial setup
	n.nodetype = AND_NODE; // virtual root is AND_NODE
	n.vartype = MAX_VAR; // virtual root is MAX_NODE
	n.id = -1; // flag
	n.ub = UB;
	n.exact_val = 0.0;
	// Caveat: must insert first, then set "root = exploredTree.begin()"!!!
	exploredTree.insert(exploredTree.begin(), n);
	root = exploredTree.begin();
	GSIZE = 1;
	std::cout<<"==== Root configuration set ====\n";
}
void mmap::setExactHeur() {
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
void mmap::printNode(const node& n) {
	/*
	 * print information of a node
	 */
//	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	auto nodetype = (n.nodetype==AND_NODE)? "AND" : "OR";
	auto vartype = (n.vartype==MAX_VAR)? "MAX" : "SUM";
	auto isExpand  = n.expanded? "expanded" : "non-expanded";
	auto solved = n.solved? "solved" : "non-solved";
	auto frontierType = (n.frontierType==MAX_VAR)? "MAX" : "SUM";
	auto worstFrontierType = (n.worstFrontierType==MAX_VAR)? "MAX" : "SUM";
	std::cout << isExpand << " " << nodetype << " node, " << vartype
			<<" variable, id "<< n.id << ", value " << n.val << ", " << solved
			<< ", ub " << n.ub << ", exact_val " << n.exact_val
//			<< ", down_ub " << n.down_ub << ", downWorst_ub " << n.downWorst_ub
			<< ", down_ub " << n.down_ub
			<< ", downWorstMax_ub " << n.downWorstMax_ub << std::endl;
	std::cout << "frontierType " << frontierType << ", worstFrontierType " << worstFrontierType << std::endl;
//	std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}
void mmap::printTree() {
	std::cout << "+++++++++++++++++++++++printing the entire tree++++++++++++++++++++++++++++++\n";
	for (auto ptr = exploredTree.begin(); ptr != exploredTree.end(); ++ptr) {
		printNode(*ptr);
	}
	std::cout << "+++++++++++++++++++++++Done printing++++++++++++++++++++++++++++++"<<std::endl;
}
void mmap::printRemovableSet() {
/*
 * print all the removable nodes
 */
	std::cout << "+++++++++++++++++++++++printing the removable set++++++++++++++++++++++++++++++\n";
//
//	std::stack<std::pair<tree<node>::iterator, double>> stack;
//	stack.push(std::make_pair(root, 0.0));
	std::stack<std::pair<tree<node>::iterator, std::vector<double> >> stack;
	std::vector<double> vec(2, 0.0);
	stack.push( std::make_pair(root, vec) );
	while(!stack.empty()){
		auto top = stack.top(); stack.pop();
		auto ptr = top.first;
		vec = top.second;
		auto max_val = vec[0];
		auto val = vec[1];
//		if ( ptr->nodetype == OR_NODE ) {
//			val += getSibUB(ptr);
//		}
		//
		//		if (exploredTree.begin(ptr)==exploredTree.end(ptr)) continue;
		if (ptr->toChild == NULL) {
			// frontier node
			// a solved MAX node is considered as a frontier node
			if (exploredTree.begin(ptr) != exploredTree.end(ptr))
				assert(ptr->solved && ptr->vartype == MAX_VAR);
			continue;
		}

		if ( ptr->nodetype == OR_NODE ) {
			val += getSibUB(ptr);
			if ( ptr->vartype == MAX_VAR )
				max_val += getSibUB(ptr);
		} else {
			val += ptr->exact_val;
			if ( ptr->vartype == MAX_VAR )
				max_val += ptr->exact_val;
		}
		bool isRemovable = true;
//		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib){
//			if (exploredTree.begin(sib) != exploredTree.end(sib)) {
//			if ( sib->toChild != NULL ) {
//				isRemovable = false;
//				if (ptr->nodetype==AND_NODE){
//					stack.push(std::make_pair(sib, val)); // sib is OR node, take val from the AND parent
//				} else {
//					stack.push(std::make_pair(sib, val+sib->exact_val));
//				}
//			}
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// add non-frontier node (solved MAX node is considered as a frontier node) to stack
			if (sib->toChild != NULL) {
				isRemovable = false;
				if (ptr->vartype != sib->vartype) {
					// ptr must be AND-MAX node, sib must be OR-SUM node
					// ptr is the last MAX node for sib
					vec = {max_val - ptr->exact_val + ptr->ub, val};
				} else {
					vec = {max_val, val};
				}
				//
				stack.push(std::make_pair(sib, vec));
			}
		}
//		}
		if (isRemovable) {
			std::cout << "\n+++++++++++++++++++++++parent++++++++++++++++++++++++++++++\n";
			printNode(*ptr);
			auto ptr_max_val = max_val;
			if ( ptr->vartype == MAX_VAR  ) {
				if ( ptr->nodetype == AND_NODE)  {
					ptr_max_val += ptr->ub - ptr->exact_val;
				} else {
					ptr_max_val += ptr->ub;
				}
			}
			// For debugging: max_val of the best frontier node should be equal to the global upper bound!!!
			// because max_val is the upper bound of the solution tree that includes the frontier node
			std::cout << "MAX solution tree ub " << ptr_max_val << std::endl;
			std::cout << "path ub " << ptr->down_ub + val << std::endl;
//			std::cout << "path worst ub " << ptr->downWorst_ub + val << std::endl;
			std::cout << "path worst max ub " << ptr->downWorstMax_ub + val << std::endl;
			std::cout << "+++++++++++++++++++++++chidlren++++++++++++++++++++++++++++++\n";
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib){
//				assert(exploredTree.begin(sib) == exploredTree.end(sib));
				assert( sib->toChild == NULL );
				printNode(*sib);
				auto sib_max_val = max_val;
				auto sib_val = val;
				//
				if (ptr->vartype != sib->vartype) {
					// ptr must be AND-MAX node, sib must be OR-SUM node
					// ptr is the last MAX node for sib
					sib_max_val = max_val - ptr->exact_val + ptr->ub;
				}

				if ( sib->nodetype == OR_NODE ) {
					sib_val += getSibUB(sib);
					if ( sib->vartype == MAX_VAR )
						sib_max_val += getSibUB(sib);
				} else {
					sib_val += sib->exact_val;
					if ( sib->vartype == MAX_VAR )
						sib_max_val += sib->exact_val;
				}

				if (sib->vartype == MAX_VAR) {
					if (sib->nodetype == AND_NODE) {
						sib_max_val += sib->ub - sib->exact_val;
					} else {
						sib_max_val += sib->ub;
					}
				}
				std::cout << "MAX solution tree ub " << sib_max_val << std::endl;
				std::cout << "path ub " << sib->down_ub + sib_val << std::endl;
//				std::cout << "path worst ub " << sib->downWorst_ub + sib_val << std::endl;
				std::cout << "path worst max ub " << sib->downWorstMax_ub + sib_val << std::endl;
//				std::cout << "+++++++++++++++++++++++chidlren++++++++++++++++++++++++++++++\n";
			}
		}
	}
	std::cout << "+++++++++++++++++++++++Done printing++++++++++++++++++++++++++++++"<<std::endl;
}
void mmap::printOPEN() {
/*
 * For debugging only
 * max_val of the best frontier node should be equal to the global upper bound!!!
 * because max_val is the upper bound of the solution tree that includes the frontier node
 */
	std::cout << "+++++++++++++++++++++++printing the OPEN set++++++++++++++++++++++++++++++\n";
	std::stack<std::pair<tree<node>::iterator, std::vector<double> >> stack;
	std::vector<double> vec(2, 0.0);
	stack.push( std::make_pair(root, vec) );
	while(!stack.empty()) {
		auto top = stack.top(); stack.pop();
		auto ptr = top.first;
		vec = top.second;
		auto max_val = vec[0];
		auto val = vec[1];
		//
		if ( ptr->nodetype == OR_NODE ) {
			val += getSibUB(ptr);
			if ( ptr->vartype == MAX_VAR )
				max_val += getSibUB(ptr);
		} else {
			val += ptr->exact_val;
			if ( ptr->vartype == MAX_VAR )
				max_val += ptr->exact_val;
		}

		if (ptr->toChild == NULL) {
			// a solved MAX node is considered as a frontier node
			if ( exploredTree.begin(ptr) != exploredTree.end(ptr) )
				assert( ptr->solved && ptr->vartype == MAX_VAR );
			//
			printNode(*ptr);
			if ( ptr->vartype == MAX_VAR  ) {
				if ( ptr->nodetype == AND_NODE)  {
					max_val += ptr->ub - ptr->exact_val;
				} else {
					max_val += ptr->ub;
				}
			}
			// For debugging: max_val of the best frontier node should be equal to the global upper bound!!!
			// because max_val is the upper bound of the solution tree that includes the frontier node
			std::cout << "MAX solution tree ub " << max_val << std::endl;
			std::cout << "path ub " << ptr->down_ub + val << std::endl;
//			std::cout << "path worst ub " << ptr->downWorst_ub + val << std::endl;
			std::cout << "path worst max ub " << ptr->downWorstMax_ub + val << std::endl;
			std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
		} else {
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if ( ptr->vartype != sib->vartype ) {
					// ptr must be AND-MAX node, sib must be OR-SUM node
					// ptr is the last MAX node for sib
					vec = {max_val - ptr->exact_val + ptr->ub, val};
				} else {
					vec = {max_val, val};
				}
				//
				stack.push(std::make_pair(sib, vec));
			}

		}
	}
	std::cout << "+++++++++++++++++++++++Done printing++++++++++++++++++++++++++++++"<<std::endl;
}
void mmap::printMaxConfig() {
/*
 * print the max configuration
 */
	int count = 0;
	int cntPerLine = 10;
	//
	std::cout << "++++++++++++++++++++++++++Printing current MAX configuration++++++++++++++++++++++++++++"<<std::endl;
	for (int i=0; i<nvar; ++i) {
		if ( varTypeVec[i] == MAX_VAR ) {
			std::cout << "(" << i << "," << maxConfig[i] << "), ";
			++count;
		}
		if (count == cntPerLine) {
			count = 0;
			std::cout << "\n";
		}
	}
	std::cout.precision(10);
	std::cout << "\n" << " MMAP <= " << UB << "\n";
	std::cout << "++++++++++++++++++++++++++Done printing+++++++++++++++++++++++++++++"<<std::endl;
}
void mmap::getMaxConfig() {
/*
 * get the max (possibly partial) configuration based on current explored Tree
 */
	std::fill(maxConfig.begin(), maxConfig.end(), -1); // reset first
	std::stack<tree<node>::iterator> stack; // only store MAX nodes
	auto ptr = root;
	stack.push(root);
	while(!stack.empty()) {
		ptr = stack.top();
		stack.pop();
		//
		assert(ptr != NULL);
		assert(ptr->vartype == MAX_VAR);
		//
		if (verbose > 2) {
			std::cout << "stack size: " << stack.size() << std::endl;
			printNode(*ptr);
		}
		if (ptr->nodetype == AND_NODE) {
//				std::cout << "no. of children (before) : " << exploredTree.number_of_children(ptr) << std::endl;
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if (sib->vartype == MAX_VAR)
					stack.push(sib);
			}
//				std::cout << "no. of children (after): " << exploredTree.number_of_children(ptr) << std::endl;
			if (ptr->id > -1) {
				maxConfig[ptr->X] = ptr->val;
			}
		} else {
			// OR_NODE
			if (ptr->toChild != NULL) {
				// debug
				if (verbose > 2) {
//					printNode(*(ptr->toChild));
//					printNode(*ptr);
//					printNode(*(exploredTree.parent(ptr->toChild)));
					assert(ptr == exploredTree.parent(ptr->toChild));
				}
				assert(!ptr->solved);
				stack.push(ptr->toChild);
			} else {
				// toChild == NULL => frontier node or solved node
				if ( exploredTree.begin(ptr) != exploredTree.end(ptr)  ) {
					// possibly a solved OR-MAX node with exactly one child
					assert(ptr->solved);
					//
					std::cout << "no. of children: "  << exploredTree.number_of_children(ptr) << std::endl;
					assert(  exploredTree.number_of_children(ptr) == 1 );
					stack.push( exploredTree.begin(ptr) );
				}
			}
		}
	}
}
// helper functions
double mmap::logsumexp(std::list<double>& vec) {
	assert(!vec.empty()); // cannot be empty
	double mx = *(std::max_element(vec.begin(), vec.end()));
	double r = 0;
	if (mx <= std::numeric_limits<double>::lowest()) {
		return mx;
	}
	for (double& val : vec) {
		r += exp(val - mx);
	}
	return (log(r) + mx);
}

double mmap::getSibUB(tree<node>::iterator ptr) {
/*
 * accumulate upper bounds of siblings, edge weight NOT included
 * children's best bounds may not correspond to the parent's best bounds
 * only for AND node
 */
	assert( ptr->nodetype == OR_NODE  && ptr != root );
	auto par = exploredTree.parent(ptr);
	double ub = 0.0;

	for (auto sib=exploredTree.begin(par); sib != exploredTree.end(par); ++sib) {
		if (sib != ptr) ub += sib->ub;
	}
	return ub;
}

bool mmap::checkSolved(const node& n) {
	// works for both AND_NODE and OR_NODE nodes
	// caveat: may only apply to frontier nodes!

	// rare case, for a MAX node, even if ub = 0, we may still expand it in order to get a full configuration of MAX variables
//	if ( n.ub <= std::numeric_limits<double>::lowest() ) {
//		return true;
//	}
	bool solved = false;
	if (n.id > -1) {
		// (ascList[n.id]).empty() == true => exactHeur[n.id] == true
		solved = (ascList[n.id]).empty();
		if (n.vartype == SUM_VAR) {
			// currently we only apply this to SUM vars
			solved = exactHeur[n.id]; // for leaf node in pseudo tree, this is always true
		}
	}
	return solved;
}
tree<mmap::node>::iterator mmap::pruneSolved(tree<node>::iterator ptr, tree<node>::iterator ancestor) {
/*
 * Must run from a newly generated node first!!!
 * The highest solved node will be kept in memory as place-holder
 * an AND node or OR-SUM node is solved only when all its children are solved;
 * an OR-MAX node is solved if its best MAX child is solved;
 * all SUM children of a solved node will be removed,
 * but all MAX children corresponding to best sub-configuration of a solved node will be kept.
 * return the highest UNsolved node or specified by "ancestor"
 * "ancestor" must be an ancestor of ptr, thus,   ptr <= ancestor <= root
 * Caveat: even if a descendant is not solved, an ancestor may still be solved!!!
 * -- the reason is that an OR-MAX node only requires its best child to be solved.
 * So, it happens if a descendant is in a sub-optimal branch, which is different from the pure summation case!!!
 */
	if (verbose > 2) {
		std::cout << "~~~~~~~~~~~~~~~~~~~~~Enter pruningSolved~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	}
//	if ( (ptr==NULL) || (!ptr->solved) ) {
	if ( ptr == NULL ) {
		return ptr;
	}
	//
//	auto highestSolved = ptr;
	tree<node>::iterator  highestSolved = NULL;
	auto startPtr = ptr; // cache the start node
	auto par = exploredTree.parent(ancestor);
	int GSIZE_pre = GSIZE;
	int cnt = 0;
	std::list<tree<node>::iterator> rmList; // all removable solved children of a solved node
	//
	while (ptr != par) {
		rmList.clear();
//		if (!ptr->solved) {
		// if labeled as unsolved, check whether solved now
		// if solved, also check its children
		ptr->solved = true; // temporarily set to be true, might be changed later on
		//
		if (ptr->nodetype == AND_NODE || ptr->vartype == SUM_VAR) {
			// an AND node or a SUM node is solved only when all its children are solved
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if (!sib->solved) {
					ptr->solved = false;
					break;
				}
				if (sib->vartype == SUM_VAR)
					rmList.push_back(sib);
			}
		} else {
			// ptr is OR-MAX node
			// an OR-MAX node is solved if its best MAX child is solved;
			double ub = -std::numeric_limits<double>::infinity();
			tree<node>::iterator toChild = NULL;
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				if (toChild == NULL || sib->ub > ub) {
					toChild = sib;
					ub = toChild->ub;
				}
			}
			assert(toChild != NULL); // debug, ptr should have at least one child because it is parent of some node.
			ptr->solved = toChild->solved;
			if (ptr->solved) {
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					if (sib != toChild)
						rmList.push_back(sib); // all non-max AND children should be deleted
				}
			}
		}
//		}
		// check whether solved now
		// Caveat: even if not solved, still go upward because there might be some OR-MAX ancestor solved
		if (ptr->solved) {
			// debug
//			if (verbose > 2) {
//				std::cout << "size of rmList: " << rmList.size() << std::endl;
//			}
			// for a solved node, set toChild to be NULL. However, solved MAX nodes may still have descendants for book keeping
			ptr->toChild = NULL; // Caveat: important to reset toChild here!!!
			ptr->toWorst = NULL;
			for (auto sib : rmList) {
				cnt = exploredTree.size(sib); // the size of the subtree rooted at given node ( including that node )
				GSIZE -= cnt;
				exploredTree.erase(sib);
			}
			highestSolved = ptr;
		}
		ptr = exploredTree.parent(ptr);
	}
	// debug
//	assert(highestSolved->solved);

	if (verbose > 2) {
		std::cout << "pruning done: " << GSIZE_pre - GSIZE << " nodes deleted!" << std::endl;
		std::cout << "~~~~~~~~~~~~~~~~~~~~~exit pruning~~~~~~~~~~~~~~~~~~~~~" <<std::endl;
	}
	//
	if ( highestSolved == NULL ) return startPtr; // none solved, return the start node
	if ( highestSolved == ancestor ) return ancestor;
	// return UNSOLVED parent if not ancestor
	return exploredTree.parent(highestSolved);
}
void mmap::updateBounds(tree<node>::iterator ptr) {
	// re-calculate bounds from children
	// run it from a newly expanded node !!!
	// bounds of ptr must have been updated already
	if (verbose > 2) {
		std::cout << "now update bounds: " << std::endl;
	}
	assert( ptr != NULL );
	ptr = exploredTree.parent(ptr);
	auto top = exploredTree.parent(root);

	double cur_ub = 0.0;
	//
	while( ptr != top ) {
		// update bounds by re-calculation.
		// more time-consuming than simple update, but more numerically stable
		if (ptr->nodetype == AND_NODE) {
			// AND node, no difference between MAX and SUM
			cur_ub = ptr->exact_val;
			// possibly some children have been pruned
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				cur_ub += sib->ub;
			}
		}
		else {
			// OR_NODE node, differentiate MAX and SUM
			// in our case, ptr->exact_val = -inf for OR-SUM,
			// a useful flag for OR-MAX
			std::list<double> ub_vec {ptr->exact_val};
			for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				ub_vec.push_back( sib->ub );
			}
			if ( ptr->vartype == SUM_VAR ) {
				cur_ub = logsumexp(ub_vec);
			} else {
				// ptr: MAX_VAR
				cur_ub =  *(std::max_element(ub_vec.begin(), ub_vec.end()));
			}
		}
		if ( cur_ub < ptr->ub ) {
			ptr->ub = cur_ub;
		} else {
			// no need to move upward
			break;
		}
		ptr = exploredTree.parent(ptr);
	}
	// update bounds for the root
	UB = root->ub;
}
void mmap::findChild(tree<node>::iterator ptr) {
/*
 * Must run pruneSolved before this function!!!
 * find the best and worst child for an UNSOLVED node
 * and update its info accordingly
 * a solved child will not be considered
 * Caveat: for a solved MAX node, it possibly has solved MAX children, but toChild = toWorst = NULL
 * A solved MAX node is considered as a frontier node although it may have descendants
 * Since we keep track of "down" ub to make sure it is the best "down" ub so far,
 * we do not have to cache "prior" any more
 * Priority of a node is defined by its best descendant visited before
 * Caveat: the best descendant and the worst descendant behave differently.
 * Briefly: the best descendant is always on the best branch, the worst may not be on the worst branch
 * Two principles for removing the worst
 * 1) remove the worst among all removable nodes: consistent with our intuition.
 * 2) remove the worst on the worst branch: consistent with the definition of priority. easy to implement, less memory
 * Here we adopt the second principle
 */
	assert(ptr!=NULL); // should not be NULL
	// if solved or no children
	if ( ptr->solved || exploredTree.begin(ptr)==exploredTree.end(ptr) ) {
		ptr->toChild = ptr->toWorst = NULL;
		return;
	}
	// for upper priority, no need to exactly compute the priority to identify which child is the best
	// just a relative quantity for comparison since we only need the order of those children
	tree<node>::iterator toChild = NULL;
	tree<node>::iterator toWorst = NULL;
	// whether ptr is in the removable set, i.e., all its children are frontier nodes (solved children are considered frontier nodes)
	bool isRemovable = true;
	//
	for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
		// never consider a solved child, even it is the only child
		// if a node only has solved children, this node must be labeled solved in pruneSolved already
		// a solved SUM node should have no children
		// a solved MAX node can not be the best child otherwise its parent is already solved
		// any solved node (except for a solved MAX node) should have no descendants
		if (sib->solved) continue;
		//
		if (toChild == NULL) {
			toChild = sib;
		}
		if ( sib->toChild != NULL ) {
			// for a solved node (may have solved descendants), toChild == NULL
			// debug
			 assert(exploredTree.begin(sib) != exploredTree.end(sib));
			//
			isRemovable = false;
			if (toWorst == NULL) toWorst = sib;
		}
		//
		if ( ptr->nodetype == OR_NODE ) {
			// AND children have the same variable type as the parent
			// case 1: Among all AND-MAX children, simply pick the one with highest upper bound
			// case 2: Among all AND-SUM children, pick the one that has the highest-upper-bound frontier descendant
			if ( ptr->vartype == MAX_VAR ) {
				// case 1:
				if ( sib->ub > toChild->ub ) toChild = sib;
				// principle 2), consistent with the priority definition, easy to implement
				if ( exploredTree.begin(sib) != exploredTree.end(sib)
						// equivalent: sib->toChild == NULL, see "continue" above
						&& sib->downWorstMax_ub + sib->exact_val < toWorst->downWorstMax_ub + toWorst->exact_val ) {
//				if ( exploredTree.begin(sib) != exploredTree.end(sib)
//						&& sib->downWorst_ub + sib->exact_val < toWorst->downWorst_ub + toWorst->exact_val ) {
//				if ( exploredTree.begin(sib) != exploredTree.end(sib)  &&  sib->ub < toWorst->ub ) {
					toWorst = sib;
				}
			} else {
				// case 2:
				if ( sib->down_ub + sib->exact_val > toChild->down_ub + toChild->exact_val )
					toChild = sib;
				if ( exploredTree.begin(sib) != exploredTree.end(sib)
//						&& sib->downWorst_ub + sib->exact_val < toWorst->downWorst_ub + toWorst->exact_val )
						&& sib->downWorstMax_ub + sib->exact_val < toWorst->downWorstMax_ub + toWorst->exact_val )
					// for a SUM node, downWorstMax_ub == downWorst_ub
					toWorst = sib;
			}
		} else {
			// ptr->nodetype == AND_NODE
			// OR children may have different variable types
			// OR children's best children belong to one single solution subtree, while worst chidlren may not be the case
			// compare their best frontier descendants: MAX frontier preferred
			// TODO MAX v.s. MAX: may have to specify a priority, but currently we just pick it from left to right
			// SUM v.s. SUM
			if ( sib->frontierType == toChild->frontierType  &&  toChild->frontierType == SUM_VAR ) {
				if ( sib->down_ub + toChild->ub > toChild->down_ub + sib->ub )
					toChild = sib;
			}
			// MAX v.s. SUM
			if ( sib->frontierType == MAX_VAR  &&  toChild->frontierType == SUM_VAR ) {
				// MAX variables are preferred for best child
				toChild = sib;
			}
			// TODO: we may need to specify a priority, but currently we just pick it from left to right, even for worst node
			if ( sib->toChild == NULL ) {
				// debug
				assert( exploredTree.begin(sib) == exploredTree.end(sib) || sib->solved );
				continue; // no children, skip the rest for the worst node
			}
/*			// SUM - SUM
			if ( sib->worstFrontierType == toWorst->worstFrontierType  &&  toWorst->worstFrontierType == SUM_VAR ) {
				if ( sib->downWorst_ub + toWorst->ub < toWorst->downWorst_ub + sib->ub )
					toWorst = sib;
			}
			// MAX - SUM
			if ( sib->worstFrontierType == SUM_VAR && toWorst->worstFrontierType == MAX_VAR ) {
				// SUM variables are preferred for worst child
				toWorst = sib;
			}*/

			// for the worst case, actually, the frontier type does not matter much
			// (only used for tie-breaking, which is quite rare for real values)
			// SUM-SUM, MAX-SUM can be combined here for toWorst
			// note for SUM nodes, downWorstMax_ub == downWorst_ub
//			if ( sib->downWorst_ub + toWorst->ub < toWorst->downWorst_ub + sib->ub )
			if ( sib->downWorstMax_ub + toWorst->ub < toWorst->downWorstMax_ub + sib->ub )
				toWorst = sib;
		}
	}
	// if not solved
	assert(toChild != NULL); // debug
	if (isRemovable) {
		// if ptr is in the removable set
		assert(toWorst == NULL); // debug
		toWorst = toChild;
	}
	// debug
	if (verbose > 2) {
		if (ptr->toWorst != NULL && ptr->toWorst != toWorst) {
			std::cout << "toWorst changed:\n";
			printNode(*ptr);
			std::cout << "previous toWorst:" << std::endl;
			printNode(*(ptr->toWorst));
			std::cout << "current toWorst:" << std::endl;
			printNode(*toWorst);
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}
	}
	//
	ptr->toChild = toChild;
	ptr->toWorst = toWorst;
	// update
	// Caveat: ptr may have different frontierType as its best child due to removeWorst!!!
	// ptr always remembers information about its best frontier descendant ever reached
	// Be careful: some data (e.g., down_ub, downWorstMax_ub) of a newly generated node (ptr) has not been properly assigned yet,
	// they are initialized to be Inf, which is important.
	if (ptr->nodetype == AND_NODE) {
		double sibUB = getSibUB(toChild);
//		ptr->down_ub = std::min(toChild->down_ub + sibUB, ptr->down_ub);
		if ( toChild->down_ub + sibUB <= ptr->down_ub ) {
			// using "<" instead of "<=" ?
			// if toChild's best frontier descendant is better than the best frontier descendant that ptr has ever seen.
			ptr->down_ub = toChild->down_ub + sibUB;
			ptr->frontierType = toChild->frontierType;
		}
		//
//		if (isRemovable) {
			// if ptr is in the removable set
//			ptr->downWorst_ub = ptr->down_ub;
//			ptr->worstFrontierType = ptr->frontierType;
//			ptr->downWorstMax_ub = ptr->ub - ptr->exact_val; // tricky, be careful
//		} else {
			if (toWorst != toChild) sibUB = getSibUB(toWorst); // re-calculate if necessary
			// Caveat: be careful!!!
//			ptr->downWorst_ub = std::min(ptr->down_ub, toWorst->downWorst_ub + sibUB); // take the min
			ptr->worstFrontierType = toWorst->worstFrontierType;
			if ( ptr->vartype != toWorst->vartype ) {
				// MAX v.s. SUM, first time to reach a MAX node
				ptr->downWorstMax_ub = ptr->ub - ptr->exact_val;
			} else {
//				ptr->downWorstMax_ub = toWorst->downWorstMax_ub + sibUB;
				if ( ptr->vartype == MAX_VAR ) {
					// MAX v.s. MAX
					ptr->downWorstMax_ub = std::min(ptr->ub - ptr->exact_val, toWorst->downWorstMax_ub + sibUB);
				} else {
					// SUM v.s.SUM
					ptr->downWorstMax_ub = std::min(ptr->down_ub,  toWorst->downWorstMax_ub + sibUB);
				}
			}
			if (isRemovable) {
				// toChild == toWorst in this case,
				// since toChild is a frontier node, which is newly generated or was once in the removable set
				// thus, we have
				// toChild->worstFrontierType == toChild->frontierType
				// toChild->downWorstMax_ub caches the best
				ptr->worstFrontierType = ptr->frontierType;
			}
//		}
	} else {
		// ptr->nodetype == OR_NODE
		// if toChild is frontier node, the following equals toChild->ub
//		ptr->down_ub = std::min(toChild->down_ub + toChild->exact_val, ptr->down_ub);
		if (toChild->down_ub + toChild->exact_val <= ptr->down_ub) {
			ptr->down_ub = toChild->down_ub + toChild->exact_val;
			ptr->frontierType = toChild->frontierType;
		}
		// for worst
//		ptr->downWorst_ub = std::min( toWorst->downWorst_ub + toWorst->exact_val, ptr->down_ub );
		if (isRemovable) {
//			ptr->downWorst_ub = ptr->down_ub;
			ptr->worstFrontierType = ptr->frontierType;

			if ( ptr->vartype == MAX_VAR ) {
				// MAX-OR
				ptr->downWorstMax_ub = ptr->ub;
			} else {
				ptr->downWorstMax_ub = ptr->down_ub;
			}

		} else {
			// Caveat: be careful!!!
//			ptr->downWorst_ub = toWorst->downWorst_ub + toWorst->exact_val;
//			ptr->downWorst_ub = std::min(ptr->down_ub, toWorst->downWorst_ub + toWorst->exact_val);
			ptr->worstFrontierType = toWorst->worstFrontierType;

			if (ptr->vartype == MAX_VAR) {
				// ptr: OR-MAX
				ptr->downWorstMax_ub = std::min(ptr->ub, toWorst->downWorstMax_ub + toWorst->exact_val);
			} else {
				// ptr: OR-SUM
				ptr->downWorstMax_ub = std::min(ptr->down_ub, toWorst->downWorstMax_ub + toWorst->exact_val);
			}
		}
	}
}
int mmap::expandBest() {
/*
 * Forward pass to expand the best node in OPEN
 */
	if (verbose > 2) {
		std::cout << "expandBest: running..." << std::endl;
		// debug
		printOPEN();
	}
//	double pathUB = 0.0; // cache all branching info and instantiated factors by the path
	std::fill(tuple.begin(), tuple.end(), -1); // re-use to accelerate, for safety, initialized to -1
	auto ptr = root; // the root is AND node
	while (ptr->toChild != NULL) {
		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		if (verbose > 2) printNode(*ptr); // debug
		//
		ptr = ptr->toChild;
		if (ptr->nodetype == AND_NODE) {
			tuple[ptr->X] = ptr->val;
		}
	}
	if (verbose > 2) {
		std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		std::cout << "Best node found!" << std::endl;
		printNode(*ptr);
	}
	if (ptr->solved) {
		std::cout << "The best frontier node has been solved, quit!" << std::endl;
		return 3;
	}
	// with removeWorst, the frontier could be either AND or OR
	if (ptr->nodetype == AND_NODE) {
		std::list<mex::graphModel::vindex> childrenList;
		if (ptr->id > -1) {
			childrenList = ascList[ptr->id];
		} else {
			// if it is the virtual root
			childrenList = rts;
		}
		// if the most promising node corresponds to a leaf node of the pseudo tree
		if (childrenList.empty()) {
			std::cout << "The best node has no children, quit!" << std::endl;
			if (verbose > 2) {
				printNode(*ptr);
				while (ptr->id > -1) {
					printNode(*ptr);
					ptr = exploredTree.parent(ptr);
				}
				printNode(*ptr);
			}
			return 2;
		}
		// update bounds later if necessary
		double cur_ub = 0.0;
		bool isProp = false;
		// if expanded before, regenerate children, no need to propagate bounds
		// Note that exact_val will never be changed after initialization (even for solved nodes)
		// Thus, we do not have to reset it for expanded node like we did before
/*		if (ptr->id > -1 && ptr->expanded) {
			ptr->exact_val = wmbUB.heuristicTheta(ptr->X, tuple);
		}*/
		cur_ub = ptr->exact_val;
		// initialization
		ptr->solved = true;
		for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
			node child;
			child.nodetype = OR_NODE;
			child.id = *it;
			child.X = wmbUB.var(child.id);
			// add variable type
			child.vartype = varTypeVec[child.id];
			child.frontierType = child.worstFrontierType = child.vartype;
			// note that down_ub = inf,
			// for OR_NODE node, exact_val = -inf is an important flag
			child.exact_val = -std::numeric_limits<double>::infinity();
			// append to the explored tree
			auto ptrChild = exploredTree.append_child(ptr, child);
			++GSIZE;
			//
			mex::Factor dUB(child.X, 0.0);
			//
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
				kid.nodetype = AND_NODE;
				kid.X = child.X;
				kid.id = child.id;
				kid.val = v;
				kid.vartype = child.vartype;
				// for frontier nodes, those types should be the same
				kid.frontierType = kid.worstFrontierType = kid.vartype;
				//
				tuple[child.X] = v;
				// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
				// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
				kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
//				kid.downPath_lb = wmbLB.heuristicIn(kid.X, tuple);
				kid.down_ub = wmbUB.heuristicIn(kid.X, tuple);
//				kid.downWorst_ub = kid.down_ub;
				kid.ub = dUB[v] = kid.down_ub + kid.exact_val;
				// both MAX var and SUM var can be assigned, if MAX var, downWorstMax_ub = downWorst_ub
				kid.downWorstMax_ub = kid.down_ub;
				// whether solved
				kid.solved = checkSolved(kid);
				if (!kid.solved) {
					ptrChild->solved = false;
				}
				// delete a solved node only when its parent is solved. book-keeping, simplify implementation
				exploredTree.append_child(ptrChild, kid);
				++GSIZE;
				// debug
				if (verbose > 2) {
					printNode(kid);
				}
			}
			if (!ptrChild->solved) {
				ptr->solved = false;
			}
			// depends on variable type
			if (ptrChild->vartype == MAX_VAR) {
				ptrChild->ub = dUB.max();
			} else {
				ptrChild->ub = dUB.logsumexp();
			}
			// AND node, not depending on variable type
			cur_ub += ptrChild->ub;
			// if solved
			// SUM node, remove all its children
			// MAX node, remove all its children but keep the best MAX child if possible
			// be careful, should be after bound propagation
			// Caveat: an OR-MAX is solved once its best child is solved instead of all its children
			// TODO revise pruneSolved to make it aware of such situation
			pruneSolved(ptrChild, ptrChild); // Must run from a newly generated node!!!
			// find best child and worst child and update parent's info
			findChild(ptrChild);
			// finally
			ptrChild->expanded = true;

			if (verbose > 2) {
				// if ptrChild is solved, down_ub, and others may not be assigned,
				// though this should not affect anything since we only used ub information of a solved node later on.
//				printNode(*ptrChild);
			}
		}

		// more stable
		if ( ptr->ub > cur_ub + epsilon || ptr == root) {
			isProp = true;
		}
		// propagate bounds if necessary
		if (isProp && !ptr->expanded) {
			ptr->ub = cur_ub;
			// ptr must be updated before entering the following function
			updateBounds(ptr);
		}
		// finally
		ptr->expanded = true;
	}
	else {
		// OR node, must be expanded before and not solved
		// not necessary to do updateBounds() and pruneSolved() here
		assert(ptr->expanded); // children must be deleted before for OR node
		auto X = ptr->X;
		// reset exact_val to make it consistent with children's bounds, necessary???
		ptr->exact_val = -std::numeric_limits<double>::infinity();
		for (size_t v = 0; v < X.states(); ++v) {
			tuple[X] = v;
			//
			node kid;
			kid.nodetype = AND_NODE;
			kid.X = X;
			kid.id = ptr->id;
			kid.val = v;
			// AND child
			kid.vartype = ptr->vartype;
			kid.frontierType  = kid.worstFrontierType =  kid.vartype; // for frontier nodes, these two types should be the same
			// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
			// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
			kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
			kid.down_ub = wmbUB.heuristicIn(kid.X, tuple);
//			kid.downWorst_ub = kid.down_ub;
			kid.ub = kid.down_ub + kid.exact_val;
			// both MAX var and SUM var can be assigned, if MAX var, downWorstMax_ub = downWorst_ub
			kid.downWorstMax_ub = kid.down_ub;
			// whether solved
			kid.solved = checkSolved(kid);
			// append
//			auto ptrKid = exploredTree.append_child(ptr, kid);
			exploredTree.append_child(ptr, kid);
			++GSIZE;
		}
	}
	//
	if (verbose > 2) {
//		std::cout << "The best frontier node after expansion" << std::endl;
//		printNode(*ptr);
//		// debug
//		if (  ptr->toChild != NULL) {
//			assert( exploredTree.parent(ptr->toChild) == ptr );
//		}
	}
	// caveat: make sure bounds have been fully updated before pruning
	ptr = pruneSolved(ptr, root); // lowest UNsolved
	return backwardUpdate(ptr);
}
int mmap::removeWorst() {
/*
 * remove the worst node in the removable set of the worst solution subtree
 * priority of a node is defined by its best descendant visited before.
 */
	if (verbose > 2) {
		std::cout << "removeWorst: starting running..." << std::endl;
		//debug
		printRemovableSet();
	}
	auto ptr = root;
	while ( ptr->toWorst != NULL ) {
		if (verbose > 2) {
			// debug
			printNode(*ptr);
		}

		// note that the virtual root is also AND_NODE with id = -1
		// first node in the loop would be OR_NODE node
		ptr = ptr->toWorst;
	}
	if (verbose>2) {
		std::cout << "the worst removable node found!" << std::endl;
		std::cout << "it's best child: " << std::endl;
		printNode(*ptr);
	}
	if ( ptr->solved ) {
		std::cout << "Warning: The worst node is solved!" << std::endl;
	}
	// can be AND or OR node
	if ( ptr == root ) {
		std::cout << "Warning: The worst node is the root, quit!" << std::endl;
		return 3; // rare case when tiny memory
	}
	//
	int gsize = GSIZE;
	assert( ptr->toWorst == ptr->toChild );
	ptr = exploredTree.parent(ptr);
	if (verbose > 2) {
		std::cout << "the worst removable node itself and all its children: " << std::endl;
		printNode(*ptr);
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
			// debug
			printNode(*sib);
			// if sib is a solved max node, it may have multiple descendants as placeholder
//			assert(exploredTree.begin(sib)==exploredTree.end(sib));
		}
	}
//	if (ptr->nodetype == OR_NODE) {
	// size gives you the size of the subtree rooted at given node ( including that node )
	GSIZE -= (exploredTree.size(ptr)-1);
	exploredTree.erase_children(ptr);
	ptr->toChild = ptr->toWorst = NULL;
	// for frontier nodes, we must make sure the following quantities are identical.
//	ptr->downWorst_ub = ptr->down_ub;
	ptr->worstFrontierType = ptr->frontierType;
	// TODO any other quantities?
	//
	if (verbose>2) {
		std::cout << gsize - GSIZE << " nodes removed!" << std::endl;
//		std::cout << "removeWorstNode: print node after children removed"<< std::endl;
//		printNode(*ptr);
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	}
	return backwardUpdate(ptr);
}

int mmap::backwardUpdate(tree<node>::iterator ptr) {
/*
 * Backward update after expanding the best
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
		findChild(ptr);
		// go upward
		if ( ptr == root ) {
			break;
		}
		ptr = exploredTree.parent(ptr);
	}

	// debug
	if (verbose > 2) {
		printRemovableSet();
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		double exact_val = 0.0;
		double max_val = 0.0;
		while ( ptr->toWorst != NULL ) {
			// note that the virtual root is also AND_NODE with id = -1
			// first node in the loop would be OR_NODE node

			if (ptr->nodetype == AND_NODE) {
				exact_val += ptr->exact_val;
			} else {
				max_val += getSibUB(ptr);
			}

			ptr = ptr->toWorst;
			printNode(*ptr);
		}
		if ( ptr->nodetype == OR_NODE ) {
			max_val += getSibUB(ptr);
		}

		std::cout << "exact_val + ptr->ub => " << exact_val << " + " << ptr->ub << " = " << exact_val + ptr->ub << std::endl;
		std::cout << "priority value => " << exact_val + ptr->ub + max_val << std::endl;
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		ptr = exploredTree.parent(ptr);
		std::cout << "the worst removable node itself and all its children: " << std::endl;
		printNode(*ptr);
		for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
				// debug
				printNode(*sib);
		}
	}

	//
	return 0;
}
void mmap::greedyCompletion() {
/*
 * greedy approach to complete current best configuration to obtain a full configuration of all MAX variables
 */
	std::fill(maxConfig.begin(), maxConfig.end(), -1); // reset first
	std::fill(tuple.begin(), tuple.end(), -1); // reset
	double greedyUB = 0.0; // upper bound of the greedily completed configuration
	std::list<mex::graphModel::vindex> childrenList;
	//
	std::stack<tree<node>::iterator> stack; // MAX nodes only
	stack.push(root);
	while(!stack.empty()) {
		auto ptr = stack.top(); stack.pop();
		// ptr must be MAX node
		assert(ptr->vartype == MAX_VAR);
		if (exploredTree.begin(ptr) != exploredTree.end(ptr))
			// has children
			if (ptr->nodetype == AND_NODE) {
				greedyUB += ptr->exact_val;
				if (ptr->id > -1) {
					maxConfig[ptr->X] = ptr->val;
					tuple[ptr->X] = ptr->val;
				}
				for (auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib) {
					if (sib->vartype == MAX_VAR)
						stack.push(sib);
					else
						greedyUB += sib->ub;
				}
			} else {
				// OR_NODE
				if (ptr->toChild != NULL) {
					assert(!ptr->solved); // debug
					stack.push(ptr->toChild);
				} else {
					// solved MAX-OR
					assert(ptr->solved && exploredTree.number_of_children(ptr) == 1); // debug
					stack.push(exploredTree.begin(ptr));
				}
			}
		else {
			// ptr has NO children
			if ( ptr->nodetype == AND_NODE  &&  ptr->id > -1) {
				maxConfig[ptr->X] = ptr->val;
			}
			if ( ptr->solved ) {
				// this node can not have any MAX descendant
				greedyUB += ptr->ub;
				continue;
			}
			// no children
			double gub = greedyUB; // temporary for greedyUB
			std::stack<node> ndStack;
			ndStack.push(*ptr);
			while (!ndStack.empty()) {
				auto nd = ndStack.top(); ndStack.pop();
				if (nd.nodetype == AND_NODE) {
					assert( nd.vartype == MAX_VAR );
					tuple[nd.X] = nd.val;
					if (nd.id > -1) {
						childrenList = ascList[nd.id];
					} else {
						// if it is the virtual root
						childrenList = rts;
					}
					// leaf node of the pseudo tree
					if (childrenList.empty()) {
//						greedyUB += nd.ub;
						gub += nd.ub;
						continue;
					}
//					greedyUB += nd.exact_val;
					gub += nd.exact_val;
					for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
						node child;
						child.nodetype = OR_NODE;
						child.id = *it;
						child.X = wmbUB.var(child.id);
						// add variable type
						child.vartype = varTypeVec[child.id];
//						child.frontierType = child.worstFrontierType = child.vartype;
//						child.exact_val = -std::numeric_limits<double>::infinity();
						ndStack.push(child);
					}
				} else {
					// OR node
					auto X = nd.X;
					std::vector<node> ndVec(X.states());
					double ub = -std::numeric_limits<double>::infinity();
					size_t vv = 0; // be safe to be 0
					mex::Factor dUB(X, 0.0);
					for (size_t v = 0; v < X.states(); ++v) {
						tuple[X] = v;
						//
						node kid;
						kid.nodetype = AND_NODE;
						kid.X = X;
						kid.id = nd.id;
						kid.val = v;
						// AND child
						kid.vartype = nd.vartype;
//						kid.frontierType  = kid.worstFrontierType =  kid.vartype; // for frontier nodes, these two types should be the same
						kid.exact_val = wmbUB.heuristicTheta(kid.X, tuple);
						kid.down_ub = wmbUB.heuristicIn(kid.X, tuple);
						//			kid.downWorst_ub = kid.down_ub;
						kid.ub = dUB[v] = kid.down_ub + kid.exact_val;
//						kid.downWorstMax_ub = kid.down_ub;
//						kid.solved = checkSolved(kid);
						ndVec[v] = kid;
						if ( kid.ub > ub ) {
							ub = kid.ub;
							vv = v;
						}
					}
					if (nd.vartype == SUM_VAR) {
//						greedyUB += dUB.logsumexp();
						gub += dUB.logsumexp();
					} else {
						// only push the best MAX-AND child
						ndStack.push( ndVec[vv] );
					}
				}
			}
			// take the min
			greedyUB = std::min(greedyUB+ptr->ub,  gub);
		}
	}
	std::cout.precision(10);
	// note that it's possible greedyUB > UB due to removeWorst.
	std::cout << "Greedy MMAP <= " << greedyUB << std::endl;
	printGreedyCompletion();
}
void mmap::printGreedyCompletion() {
	// printing function of greedyCompletion
	int count = 0;
	int cntPerLine = nvar; // may be 100?
	std::cout << "++++++++++++++++++++++++++Printing current MAX configuration with greedy completion++++++++++++++++++++++++++++"<<std::endl;
	for (int i=0; i<nvar; ++i) {
		if ( varTypeVec[i] == MAX_VAR ) {
			if ( maxConfig[i] != -1 ) {
				// Caveat: When comparing signed with unsigned, the compiler converts the signed value to unsigned.
				std::cout << "(" << i << "," << maxConfig[i] << "), ";
			}
			else {
				assert( tuple[i] != -1 ); // debug
				std::cout << "(" << i << ",[" << tuple[i] << "]), ";
			}
			++count;
		}
		if (count == cntPerLine) {
			count = 0;
			std::cout << "\n";
		}
	}
	std::cout << "\n++++++++++++++++++++++++++Done printing+++++++++++++++++++++++++++++"<<std::endl;
}
void mmap::start() {
	std::cout << "Search time limit (sec): " << timeBudget << std::endl;
	std::cout << "Search memory limit (MB): " << memoryBudget << std::endl;
//
	int memUnit = sizeof(*root);
	// MemLimit is in megabyte
	auto treeSizeLimit = unsigned(1024 * 1024 * memoryBudget / memUnit);
	// debug
//	treeSizeLimit = 50;
	//
	std::cout << "Tree size limit: " << treeSizeLimit << std::endl;

	double outFrequency = 0.1; // for output
	//
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startTime << "]: " << std::setw(10)
			<< " MMAP < " << UB << "\n";
	std::cout << "Tree size: " << GSIZE << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			<< std::endl;
	//
	int indicator = 0;
	bool isFirstReachMemLimit = false;
	while (true) {
		if (root->solved) {
			std::cout << "The root is solved, quit!" << std::endl;
			break;
		}
		if (mex::timeSystem() - startTime > outFrequency) {
			outFrequency += outFrequency; // exponential
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startTime << "]: "
					<< std::setw(10) << " MMAP < " << UB << "\n";
			std::cout << "Tree size: " << GSIZE << std::endl;
			greedyCompletion();
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
					<< std::endl;
			//
		}
		if (GSIZE < treeSizeLimit) {
			// debug
			if (verbose > 2) {
				std::cout
						<< "~~~~~~~~~~~~~~~~~before expanding the best frontier~~~~~~~~~~~~~~~~~\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
						<< std::setw(10) << " MMAP < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE << "\n";
				std::cout
						<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
						<< std::endl;
				std::cout << "root information:" << std::endl;
				printNode (*root);
				std::cout
						<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
						<< std::endl;
			}
			indicator = expandBest();
			if (indicator != 0) {
				//				std::cout<< "The best frontier has no children, loop terminated!\n";
				break;
			}
			if (verbose > 2) {
				std::cout
						<< "~~~~~~~~~~~~~~~~~after expanding the best frontier~~~~~~~~~~~~~~~~~\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
						<< std::setw(10) << " MMAP < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE << "\n";
				std::cout << "root information:" << std::endl;
				printNode (*root);
				std::cout
						<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
						<< std::endl;
			}
		} else {
			// reach tree size limit (memory limit)
			if (!isFirstReachMemLimit) {
				// only prompt at the first reach memory limit
				std::cout << "First time reach tree size limit: " << treeSizeLimit << " ( memory limit: "<< memoryBudget << " MB )\n";
				std::cout << "[" << mex::timeSystem() - startTime << "]: " << std::setw(10) << " MMAP < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE << std::endl;
				greedyCompletion();
				std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
				isFirstReachMemLimit = true;
			}


			if (verbose > 2) {
				std::cout
						<< "~~~~~~~~~~~~~~~~~before removing the worst~~~~~~~~~~~~~~~~~\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
						<< std::setw(10) << " MMAP < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE << "\n";
				std::cout << "root information:" << std::endl;
				printNode (*root);
				std::cout
						<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
						<< std::endl;
			}

			indicator = removeWorst();
			if (verbose > 2) {
				std::cout
						<< "~~~~~~~~~~~~~~~~~after removing the worst~~~~~~~~~~~~~~~~~\n";
				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startTime << "]: "
						<< std::setw(10) << " MMAP < " << UB << "\n";
				std::cout << "Tree size: " << GSIZE << "\n";
				std::cout
						<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
						<< std::endl;
				std::cout << "root information:" << std::endl;
				printNode (*root);
				std::cout
						<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
						<< std::endl;
			}

		}
		if ( mex::timeSystem() - startTime > timeBudget ) {
			std::cout << "Reach search time limit (sec): "<< timeBudget  << ", quit!" << std::endl;
			break;
		}
	}
	if (verbose > 2) {
		getMaxConfig();
		printMaxConfig();
	}
	std::cout << "[" << mex::timeSystem() - startTime << "]: "
			<< std::setw(10) << " MMAP < " << UB << "\n";
	std::cout << "Tree size: " << GSIZE << std::endl;
	greedyCompletion();
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			<< std::endl;
}
// EOF

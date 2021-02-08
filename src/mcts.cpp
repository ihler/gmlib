/*
 * mcts.cpp
 *
 *  Created on: May 26, 2016
 *      Author: qlou
 */
#include "mcts.h"
//
mcts::mcts(mex::wmbe& wmbeub, mex::wmbe& wmbelb, int ns, double dlt, double ub, double lb):
wmbUB(wmbeub), wmbLB(wmbelb),nSample(ns), DELTA(dlt) {
	order = wmbeub.getOrder();
	priority = wmbeub.getPriority();
	nvar = wmbeub.nvar();
	GSIZE = 0;
	UB = ub;
	LB = lb;
	cur_nSample = 0;
	verbose = 1;
	setRootConfig();
}

mcts::~mcts() {
	std::cout<<" ==== Done MCTS! ====\n";
}
void mcts::setRootConfig() {
	// initialize root
	// this is a virtual root node.
	node rt;
	rt.type = AND_NODE;
	rt.depth = -1;
	rt.delta = DELTA;
	rt.ub = rt.UB = UB;
	rt.lb = rt.LB = LB;
	exploredTree.insert(exploredTree.begin(), rt);
	std::cout<<" ==== Root configuration set! ====\n";
}
void mcts::sample(int num) {
	// sample given number of times
	while(num>0){
		--num;
		sampleOne();
	}
}
void mcts::sampleOne() {
	++cur_nSample;
	// current it only supports sampling from the root
	mex::vector<uint32_t> xhat(nvar);
	std::vector<double> condUB(nvar+1, 0.0);
	std::vector<double> thetaVec(condUB);
	// sample, get xhat, condUB, thetaVec, and logQxVec
	auto logQxVec = wmbUB.sampleMixtureMore( xhat, condUB, thetaVec);
	// add sample to the explored tree
	addSampleToTree(xhat, condUB, logQxVec, thetaVec);
	// normalized weight
//	double Ex = std::exp( logPQ.first - logPQ.second - UB );
}
void mcts::addSampleToTree(const mex::vector<uint32_t>& xhat, const std::vector<double>& condUB,
		const std::vector<double>& logQxVec, const std::vector<double>& thetaVec) {
	// currently only support OR structure
	double curUB = condUB.back(); // upper bound on logZ, should always be the same
	double logQx = logQxVec.back(); // logQxVec[nvar] is log{q(X)}
	double logFx = thetaVec[0]; // exact value of that configuration
	// virtual root
	auto ptr = exploredTree.begin();
	// double check, these two quantities should be the same
	assert( (UB < curUB+1e-5) && (UB > curUB-1e-5) );
	// add the sample to the virtual root
	addSampleToNode(ptr, logFx-logQx);
	calcNodeBias(ptr);
	// if the virtual root has no child yet
	if ( exploredTree.begin(ptr) == exploredTree.end(ptr)  ) {
		node n;
		n.type = OR_NODE;
		n.depth = 0;
		n.X = wmbUB.var(order[nvar - 1]);
		n.delta = DELTA;
		n.ub = n.UB = UB; // caveat: use the initialized value of UB
		n.lb = n.LB = LB;
		exploredTree.append_child(exploredTree.begin(), n);
	}

	ptr = exploredTree.begin(ptr);


	//
//	mex::vector<uint32_t> tuple(nvar);


//	double cost = 0.0;
	for (int i=nvar-1; i>=0; --i) {
		auto X = wmbUB.var( order[i] ); // reverse order
		assert( (X == (*ptr).X)  &&  ( (*ptr).type==OR_NODE) ); // safety guard
		// add AND child node, current variable
//		double ub = condUB[i]; // conditioned on the AND child
//		tuple[X] = xhat[X];
//		cost += wmbUB.heuristicTheta(X, xhat);
		// add AND child node
		ptr = addChild(ptr, X, xhat, logFx - logQxVec[i], thetaVec[i]);
		// add OR child node, associated with next variable
		if ( i > 0 ) {
			ptr = addChild(ptr, wmbUB.var( order[i-1] ), xhat);
		}
	}
}
tree<mcts::node>::iterator mcts::addChild(tree<node>::iterator ptr, mex::Var X, const mex::vector<uint32_t>& config, double wt, double cost) {
	// add one count only if current node is an OR node and its AND node is added here
	if ( (*ptr).type == OR_NODE ) (*ptr).nsp++;
	auto ptrKid = ptr;
	bool isExist = false;
	int val = config[X];
	// check whether corresponding child exists
	for (auto it = exploredTree.begin(ptr); it != exploredTree.end(ptr); ++it) {
		if ((*ptr).type == OR_NODE) {
			if (val == (*it).val) {
				// add the sample to the node and re-compute bias term
				addSampleToNode(it, wt);
				calcNodeBias(it);
				//
				ptrKid = it;
				isExist = true;
				break;
			}
		} else {
			// for AND_NODE, only one child. if no child, won't be in the loop
			assert( exploredTree.number_of_children(ptr) == 1 );
//			(*it).nsp++;
			// no need to add a count to the OR child,
			// an OR child gets one count only when one of its children gets one count
			ptrKid = it;
			isExist = true;
			break;
		}
	}
	// if no corresponding child exists, create one and update necessary information
	if (!isExist) {
		node n;
		if ((*ptr).type == OR_NODE) {
			n.type = AND_NODE;
			n.depth = (*ptr).depth; // same depth as the OR parent
			n.X = (*ptr).X;
			n.val = val;
			n.nsp = 1;
			n.logEx = wt;
			n.logEx2 = 2*wt; // initialization, in log
			n.cost = cost;
					//(*ptr).cost + wmbUB.heuristicTheta(n.X, config);
			// GSIZE counts the no. of AND nodes in the explored tree
			++GSIZE;

		} else {
			// do not add sample to OR_NODE at this moment
			n.type = OR_NODE;
			n.depth = (*ptr).depth + 1;
			n.X = X;
			n.cost = (*ptr).cost; // no cost on AND->OR arc
		}
		//
		setBounds(n, config); // set deterministic bounds
		ptrKid = exploredTree.append_child(ptr, n);
		setConfidence(ptrKid);
		if ( (*ptrKid).type == AND_NODE ) calcNodeBias(ptrKid);

		// debug
		if ((*ptrKid).type == OR_NODE) {
			assert( ! ( (*ptrKid).ub > (*ptr).ub + 1e-5  || (*ptrKid).lb < (*ptr).lb - 1e-5 ) );
//			if ((*ptrKid).ub > (*ptr).ub + 1e-5  || (*ptrKid).lb < (*ptr).lb - 1e-5 ) {
//				//
//				std::cout << "something wrong here!\n";
//				auto top = ptrKid;
//				std::cout << "kid's depth = " << (*top).depth << ", id (val) = "
//						 << " (" << (*top).val
//						<< ") , delta = " << (*top).delta << "\n";
//				std::cout << "cost = " << (*top).cost << ", lb = " << (*top).lb
//						<< ", ub = " << (*top).ub << "\n";
//				top = ptr;
//				std::cout << "par's depth = " << (*top).depth << ", id (val) = "
//						<< " (" << (*top).val
//						<< ") , delta = " << (*top).delta << "\n";
//				std::cout << "cost = " << (*top).cost << ", lb = " << (*top).lb
//						<< ", ub = " << (*top).ub << "\n";
//			}
		}
	}
	// re-compute bias and estimator
	// return
	return ptrKid;
}
void mcts::addSampleToNode(tree<node>::iterator ptr, double wt) {
	// defined for both AND and OR nodes
	++(*ptr).nsp;
	// for OR nodes, only add one count and done
	if ( (*ptr).type == OR_NODE ) return;


	if ( (*ptr).nsp <= 1 ) {
		(*ptr).logEx = wt;
		(*ptr).logEx2 = 2*wt;
		return;
	}

	int nsp = (*ptr).nsp;
	double logEx = (*ptr).logEx;
	double logEx2 = (*ptr).logEx2;
	// to be numerically stable
	// Ex  = (Ex * (samp-1))/samp + dEx/samp;
	(*ptr).logEx = logsumexp(  { std::log(nsp-1) + logEx, wt } ) - std::log(nsp);
	// Ex2 = (Ex2 * (samp-1))/samp + dEx*dEx/samp;
	(*ptr).logEx2 = logsumexp(  { std::log(nsp-1) + logEx2, 2*wt } ) - std::log(nsp);
}
void mcts::calcNodeBias(tree<node>::iterator ptr) {
	// defined for AND and OR nodes

	// if OR nodes, do nothing
	if ((*ptr).type == OR_NODE ) return;

//	assert( (*ptr).type == AND_NODE );
	if ( (*ptr).nsp <= 1 ) {
		// deterministic bound, Z <= Zhat + Zwmb
		(*ptr).logBias = (*ptr).ub;
		return;
	}

	int nsp = (*ptr).nsp;
	double confidence = (*ptr).delta;
	double Ex = std::exp( (*ptr).logEx - (*ptr).ub  ); // normalized
	double Ex2 = std::exp( (*ptr).logEx2 - 2*(*ptr).ub  ); // normalized

	double var = std::max(Ex2 - Ex*Ex, 0.0); // to be numerically stable
	// should be the unbiased sample variance
	var *= ( nsp/(nsp-1) );

	double rng = std::sqrt(2*var*std::log(2.0/confidence)/nsp) + 7*std::log(2.0/confidence)/3.0/(nsp-1);

	(*ptr).logBias = std::log(rng) + (*ptr).ub;
}
std::pair<double, double> mcts::calcEBB(tree<node>::iterator ptr) {
	// currently only defined for AND nodes
	assert( (*ptr).type == AND_NODE );
	// calculate bias term first
	calcNodeBias(ptr);
	//
	double lb = -std::numeric_limits<double>::infinity(), ub = std::numeric_limits<double>::infinity();
	if ( (*ptr).nsp > 1 ) {
		ub = logsumexp( {(*ptr).logBias,  (*ptr).logEx} );
		if ( (*ptr).logEx > (*ptr).logBias) {
				lb = std::log( std::exp( (*ptr).logEx ) - std::exp( (*ptr).logBias ) );
		}
	}
	return std::pair<double,double>(lb,ub);
}
double mcts::setConfidence(tree<node>::iterator ptr) {
	// assign confidence to the node
	// deterministic assignment s.t. sum{ Z_i * log(2/delta_i) } is minimized
	// i.e., delta_i is proportional to Z_i / sum(Z_i).
	double confidence = 0.0;
	auto ptrPar = exploredTree.parent(ptr);

	if ( (*ptr).type == OR_NODE ) {
		confidence = (*ptrPar).delta;
	} else {
		// assignment according to the deterministic upper bounds
		confidence = (*ptrPar).delta * std::exp((*ptr).ub - (*ptrPar).ub);
	}
	(*ptr).delta = confidence;
	return confidence;
}
void mcts::setBounds(node& n, mex::vector<uint32_t> tuple) {
	// set initial deterministic bounds
	// caveat: must be after n.cost is set.
	auto X = n.X;
	if (n.type == AND_NODE) {
		// AND node
		n.ub = n.cost + wmbUB.heuristicIn(X,tuple); // n.cost contains wmbUB.heuristicTheta(X,tuple) for AND nodes
		n.lb = n.cost + wmbLB.heuristicIn(X,tuple); // n.cost contains wmbLB.heuristicTheta(X,tuple) for AND nodes
		// debug
//		assert( wmbUB.heuristicTheta(X,tuple) - 1e-8 <=  wmbLB.heuristicTheta(X,tuple) );
//		assert( wmbUB.heuristicTheta(X,tuple) + 1e-8 >=  wmbLB.heuristicTheta(X,tuple) );
	}
	else {
		// OR node
		mex::Factor dUB(X,0.0), dLB(X,0.0);
		for (size_t v = 0; v < X.states(); ++v) {
			tuple[ X ] = v;
			dUB[v] = wmbUB.heuristicIn(X,tuple) + wmbUB.heuristicTheta(X,tuple); //
			dLB[v] = wmbLB.heuristicIn(X,tuple) + wmbLB.heuristicTheta(X,tuple); //
		}
		n.ub = n.cost + dUB.logsumexp();
		n.lb = n.cost + dLB.logsumexp();
	}
	// set initial bounds
	n.UB = n.ub;
	n.LB = n.lb;
}
// helper functions
double mcts::logsumexp(std::vector<double>& vec) {
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
double mcts::logsumexp(std::initializer_list<double> vec) {
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
int mcts::aggregateBounds(int max_depth) {
	// aggregate bounds from AND nodes with a specified depth
	// all those AND nodes form a cutset
	// traverse the explored tree in BFS, stop at OR node without all states being instantiated
	// do BFS
	int cutsetSize = 0; // the size of the cutset (a subset of AND nodes)
	double meanDepth = 0.0; // mean depth of nodes in the cutset
	double logZhat = -std::numeric_limits<double>::infinity(); // the aggregated estimator of logZ
	double proUB = -std::numeric_limits<double>::infinity(); // probabilistic upper bound on  logZ
	double proLB = -std::numeric_limits<double>::infinity(); // probabilistic lower bound on  logZ
	int curMaxDepth = -1;
	double aggreDelta = 0.0; // should be the same as DELTA, double check

	double nProUb = 0; // count times of taking the probabilistic upper bound over the deterministic upper bound
	double nProLb = 0;

	auto rt = exploredTree.begin(); // root is the virtual AND node with depth -1
	std::queue<tree<node>::iterator> queue;
	queue.push(rt);
	while( !queue.empty() ) {
		auto top = queue.front();
		queue.pop();
		curMaxDepth = std::max( (*top).depth, curMaxDepth );

		if ( (*top).depth > max_depth ) {
			// reach the required depth
			curMaxDepth = max_depth;
			std::cout << "reach max depth: " << max_depth << ", quit!\n";
			break;
		}

		bool doAggregation = false;
		// if an OR node does not have all the states of the associated variable being instantiated,
		// we stop here and trace back to its parent and add its parent to the cutset
		if ( ( (*top).type == OR_NODE )  && ( exploredTree.number_of_children(top) < ((*top).X).states() ) ) {
			// this node may not reach the desired depth
			top = exploredTree.parent(top);
			doAggregation = true;
		}
		// AND node with the desired depth is a member of the cutset.
		if ( ((*top).type == AND_NODE)  && ((*top).depth == max_depth) ) {
			doAggregation = true;
		}
		// an AND leaf node should also be added.
		// its OR parent must have all states instantiated already, o.w. none of its children would be added to the queue
		// note that all leaf nodes are AND nodes
		if ( ((*top).type == AND_NODE)  && (exploredTree.begin(top) == exploredTree.end(top)) ) {
			doAggregation = true;
		}

		if (doAggregation) {
			++cutsetSize;
			aggreDelta += (*top).delta;
			meanDepth = ( meanDepth*(cutsetSize-1) + (*top).depth )/cutsetSize;

			logZhat = logsumexp({ logZhat, (*top).logEx });

			// take the better between the probabilistic bounds and the deterministic bounds
			double ub = logsumexp( {(*top).logBias,  (*top).logEx} );
			if (ub < (*top).ub) nProUb++;

			ub = std::min(ub, (*top).ub);
			proUB = logsumexp( {proUB, ub} );
			// probabilistic lower bound is possibly trivial
			if ( (*top).logEx >  (*top).logBias) {
				double lb = std::log( std::exp( (*top).logEx ) - std::exp( (*top).logBias ) );

				if (lb > (*top).lb) nProLb++;

				lb = std::max(lb, (*top).lb);
				proLB = logsumexp( {proLB, lb} );
			}
			else {
				proLB = logsumexp( { proLB, (*top).lb } );
			}
//				assert ( (*top).type == AND_NODE );
//			// skip the adding children part
//			std::cout << "top's depth = " << (*top).depth << ", id (val) = "<< order[ (*top).X ] << " (" << (*top).val
//					<< ") , delta = " << (*top).delta << "\n";
//			std::cout << "cost = " << (*top).cost << ", lb = " << (*top).lb << ", ub = " << (*top).ub << "\n";


			continue;
		}
		// add children to the queue
		// note that an AND node will be added only if the states of that variable are all instantiated
		for (auto it = exploredTree.begin(top); it != exploredTree.end(top); ++it) {
			queue.push(it);
		}
	}
	assert( (aggreDelta < DELTA + 1e-5) &&  (aggreDelta > DELTA - 1e-5) ); // sanity check
	std::cout << "aggreDelta: " << aggreDelta << " == "  << DELTA << "\n";
	std::cout << "mean depth v.s. max depth reached v.s. max depth desired: " << meanDepth << "/" << curMaxDepth <<"/" << max_depth <<"\n";
	std::cout << "cutset size v.s. tree size: " << cutsetSize <<"/" << GSIZE <<"\n";
	std::cout << "count taking probabilistic lower bound: " << nProLb << "/" << cutsetSize << " = " << nProLb/cutsetSize <<"\n";
	std::cout << "count taking probabilistic upper bound: " << nProUb << "/" << cutsetSize << " = " << nProUb/cutsetSize <<"\n";
	std::cout << "Empirical Bernstein bounds: " << std::setw(10) << proLB << " < " << logZhat << " < " << proUB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	// this is the largest depth we've reached
	return curMaxDepth;
}
int mcts::aggregateBounds(int nsp_thresh, bool flag) {
	// aggregate bounds according to no. of samples on that node
	// pick AND/OR nodes with no. of samples just merely above nsp_thresh

	int cutsetSize = 0; // the size of the cutset (a subset of AND nodes)
	double meanDepth = 0.0; // mean depth of nodes in the cutset
	int maxDepth = -1; // max depth reached
//	double logZhat = -std::numeric_limits<double>::infinity(); // the aggregated estimator of logZ
	double proUB = -std::numeric_limits<double>::infinity(); // probabilistic upper bound on  logZ
	double proLB = -std::numeric_limits<double>::infinity(); // probabilistic lower bound on  logZ
	double aggreDelta = 0.0; // should be the same as DELTA, double check
	double non_aggreDelta = 0.0; // non-used part of confidence

	double meanNsp = 0.0;
	// count times of taking the probabilistic upper bound over the deterministic upper bound
	double nProUb = 0.0, nProLb = 0.0;

	auto rt = exploredTree.begin(); // root is the virtual AND node with depth -1
	std::queue<tree<node>::iterator> queue;
	queue.push(rt);
	while( !queue.empty() ) {
		auto top = queue.front();
		queue.pop();
		//
		maxDepth = std::max( (*top).depth, maxDepth );

		if ( (*top).type == OR_NODE ) {
			auto par = exploredTree.parent(top);
			assert(  (*par).nsp  == (*top).nsp );
		}

		// only aggregate AND nodes just below the threshold, and leaf AND nodes
		if ( (*top).nsp <= nsp_thresh  || ( (*top).type==AND_NODE  && exploredTree.begin(top)==exploredTree.end(top) ) ) {
			// OR_NODE won't be added if its AND parent is already below the threshold
			assert( (*top).type == AND_NODE );
			++ cutsetSize;
			aggreDelta += (*top).delta;
			meanDepth = ( meanDepth*(cutsetSize-1) + (*top).depth )/cutsetSize;
			meanNsp = ( meanNsp*(cutsetSize-1) + (*top).nsp )/cutsetSize;

			// take the better between the probabilistic bounds and the deterministic bounds
			auto ebb = calcEBB(top);
			// probabilistic lower bound is possibly trivial
			double lb = ebb.first, ub = ebb.second;
			if (lb > (*top).lb)
				nProLb++;

			lb = std::max(lb, (*top).lb);
			proLB = logsumexp( { proLB, lb });

//			proLB = logsumexp(  {proLB, (*top).lb} );

			if (ub < (*top).ub)
				nProUb++;

			ub = std::min(ub, (*top).ub);
			proUB = logsumexp( { proUB, ub });

//			proUB = logsumexp( {proUB, (*top).ub} );

			// no need to add its children
			continue;
		}
		//
		double ub = -std::numeric_limits<double>::infinity(), lb = -std::numeric_limits<double>::infinity();
		for (auto it = exploredTree.begin(top); it != exploredTree.end(top); ++it) {
			ub = logsumexp( {ub, (*it).ub} );
			lb = logsumexp( {lb, (*it).lb} );
			queue.push(it);
		}
		// add bounds of those states not visited by any sample so far
		if ( ( (*top).type == OR_NODE )  && ( exploredTree.number_of_children(top) < ((*top).X).states() ) ) {
			double diff = exp( (*top).ub ) - exp(ub);
//			assert( diff >= 0  );
			diff = log( std::max(diff, 0.0) );
//
//			diff = diff>0 ? std::log(diff) : -std::numeric_limits<double>::infinity();
			proUB = logsumexp( {proUB, diff} );
//
			non_aggreDelta += (*top).delta * exp(diff - (*top).ub);
//
			if (lb > std::numeric_limits<double>::lowest() ) {
				diff = std::exp( (*top).lb ) - std::exp(lb);
//				assert( diff >= 0  );
				diff = log( std::max(diff, 0.0) );

				proLB = logsumexp( {proLB, diff} );
			} else {
				proLB = logsumexp( {proLB, (*top).lb} );
			}
		}
	}
	assert( (aggreDelta + non_aggreDelta < DELTA + 1e-5) &&  (aggreDelta + non_aggreDelta > DELTA - 1e-5) ); // sanity check
	std::cout << "aggreDelta + non_aggreDelta: " << aggreDelta << " + " << non_aggreDelta << " = " << aggreDelta + non_aggreDelta << " == "  << DELTA << "\n";
	std::cout << "mean nsp v.s. nsp threshold: " << meanNsp << "/" << nsp_thresh <<"\n";
	std::cout << "mean depth v.s. max depth reached: " << meanDepth << "/" << maxDepth <<"\n";
	std::cout << "cutset size v.s. tree size: " << cutsetSize <<"/" << GSIZE <<"\n";
	std::cout << "count taking probabilistic lower bound: " << nProLb << "/" << cutsetSize << " = " << nProLb/cutsetSize <<"\n";
	std::cout << "count taking probabilistic upper bound: " << nProUb << "/" << cutsetSize << " = " << nProUb/cutsetSize <<"\n";
	std::cout << "Empirical Bernstein bounds: " << std::setw(10) << proLB << " < " << proUB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	return maxDepth;
}
/*
* two-stage sampling
*/
void mcts::updatePathBounds(tree<node>::iterator ptr) {
/*
 * update deterministic bounds from the leaf to the root
 */
	auto ptrPar = ptr;
	while( ptr != exploredTree.begin() ) {
		ptrPar = exploredTree.parent(ptr);
		if ( (*ptr).type == AND_NODE ) {
			std::vector<double> UBvec, LBvec;
			for(auto sib = exploredTree.begin(ptrPar); sib != exploredTree.end(ptrPar); ++sib) {
				UBvec.push_back( (*sib).UB );
				LBvec.push_back( (*sib).LB );
			}
			(*ptrPar).UB = logsumexp(UBvec);
			(*ptrPar).LB = logsumexp(LBvec);
		}
		else {
			// OR node has the same bounds as its parent.
			// caveat: not the same case when we implement the real AND/OR structure
			(*ptrPar).UB = (*ptr).UB;
			(*ptrPar).LB = (*ptr).LB;
		}
		ptr = ptrPar;
	}

}
void mcts::addSampleToPath(tree<node>::iterator ptr, const mex::vector<uint32_t>& config, const std::vector<double>& logQxVec, double logFx) {
/*
 * add one sample to the path defined by itself
 */
	auto rtPar = exploredTree.parent( exploredTree.begin() );
	int id = -1; // variable id
	int loc = nvar; // position in the order
	double wt = 0.0;

	while( ptr != rtPar ) {
		if ( (*ptr).depth<0 ) {
			// virtual root
			loc = nvar;
		} else {
			id = ((*ptr).X).label();
			loc = priority[id];
		}
		wt = logFx - logQxVec[loc];
		addSampleToNode(ptr, wt);
//		calcNodeBias(ptr); // not necessary here, we calculate it when needed
		ptr = exploredTree.parent(ptr);
	}
}
tree<mcts::node>::iterator mcts::expandTree(tree<node>::iterator ptr, const mex::vector<uint32_t>& config,
		const std::vector<double>& logQxVec, double logFx) {
/*
 * expand tree from current leaf node if necessary
 */
	assert( (*ptr).type == AND_NODE );
	assert( exploredTree.begin(ptr) == exploredTree.end(ptr) );

	int loc = nvar; // default, for the virtual root
	if ( (*ptr).depth > -1 ) {
		loc = priority[ ((*ptr).X).label() ];
	// debug
	 	assert( ((*ptr).X).label() == order[loc] );
	}
	--loc; // child order
	if (loc < 0) return ptr; // already reach the first variable in the order
	//
	auto X = wmbUB.var(order[loc]);
	int val = config[X];
	// OR node
	node n;
	n.type = OR_NODE;
	n.depth = (*ptr).depth + 1;
	n.X = X;
	n.cost = (*ptr).cost; // no cost on AND->OR arc
	setBounds(n, config); // set deterministic bounds
	ptr = exploredTree.append_child(ptr, n);
	double wt = logFx - logQxVec[loc];
	addSampleToNode(ptr); // n.nsp will be added here
	// AND nodes
	auto sib = ptr;
	auto tuple = config;
	for (int i=0; i<X.states(); ++i) {
		node kid;
		kid.type = AND_NODE;
		kid.depth = n.depth;
		kid.X = n.X;
		kid.val = i;
		tuple[kid.X] = kid.val;
		kid.cost = n.cost + wmbUB.heuristicTheta(kid.X, tuple);
		setBounds(kid, tuple); // set deterministic bounds
		sib = exploredTree.append_child(ptr, kid);
		if (i == val) {
			addSampleToNode(sib, wt);
		}
		++GSIZE;
	}
	// update path deterministic bounds
	updatePathBounds(ptr);
	//
	return ptr;
}

tree<mcts::node>::iterator mcts::upperBoundBasedSampling(tree<node>::iterator ptr){
/*
 * sample for an OR node based on current upper bounds of its AND children;
 * for each AND child, probability to be picked is proportional to its current upper bound.
 *
 */
	assert( (*ptr).type == OR_NODE );
	auto X = (*ptr).X;
	auto ns = X.states();
	assert( ns==exploredTree.number_of_children(ptr) );

	mex::Factor fc(X, 0.0);

	int i = 0;
	for ( auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib ) {
		// !!!caveat: factor's sample() function is NOT in log!!!
		// do NOT assume factors are in log!!!
		fc[i] = exp( (*sib).UB );

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

void mcts::twoStageSampling(int NSP) {
/*
 * two-stage sampling procedure
 * for internal nodes, using bound based sampling; for leaf nodes, using wmb mixture proposal to sample
 * caveat: for leaf nodes, bounds are their initial WMB bounds
*/
	mex::vector<uint32_t> config(nvar);
	std::vector<double> logQxVec(nvar+1, 0.0);
	double logFx = 0.0;
	// sample the path from root to leaf
	auto ptr = exploredTree.begin();
	auto ptrPar = ptr;
	int loc = nvar;



	// stage 1:
	while ( exploredTree.begin(ptr) != exploredTree.end(ptr) ) {
		// leaf nodes must be AND nodes
		if ( (*ptr).type == OR_NODE  ) {
			// return (randomly) picked AND child
			ptrPar = ptr;
			ptr = upperBoundBasedSampling(ptr); // return pointer to the AND child
			config[ (*ptr).X ] = (*ptr).val;

			loc = priority[ ((*ptr).X).label() ];
			// logQxVec[j] = log q(x_{j-1} | par(x_{j-1}) ) for initialization, see sampleMixtureConditioned for details
			logQxVec[loc+1] = (*ptr).UB - (*ptrPar).UB;
		}
		else {
			// AND node has at most one child
			ptr = exploredTree.begin(ptr);
		}
	}
	assert(  (*ptr).type == AND_NODE );

	// stage 2: sample from the leaf according to the mixture proposal
	int ID = -1; // default, root id
	if ( (*ptr).depth > -1  ) ID = ((*ptr).X).label();
//
	logFx = wmbUB.sampleMixtureConditioned(ID, config, logQxVec); // config, logQxVec, logFx will be updated
	// update information from current node to the root
	addSampleToPath(ptr, config, logQxVec, logFx);
	// expand tree if necessary, any other meaningful criterion?
	if ( (*ptr).nsp > NSP ) expandTree(ptr, config, logQxVec, logFx);
}

std::pair<double, double> mcts::calcRootEBB(){
	// calculate EBB for the root, since we can use the refined upper bound for the root
	auto ptr = exploredTree.begin(); // virtual root

	double lb = -std::numeric_limits<double>::infinity(), ub = std::numeric_limits<double>::infinity();
////	calcNodeBias
	if ( (*ptr).nsp <= 1 ) {
		// deterministic bound, Z <= Zhat + Z^{+}
		(*ptr).logBias = (*ptr).UB;
		return std::pair<double,double>(lb,ub);
	}

	int nsp = (*ptr).nsp;
	double confidence = (*ptr).delta;
	double Ex = std::exp( (*ptr).logEx - (*ptr).UB  ); // normalized
	double Ex2 = std::exp( (*ptr).logEx2 - 2*(*ptr).UB  ); // normalized

	double var = std::max(Ex2 - Ex*Ex, 0.0); // to be numerically stable
	// should be the unbiased sample variance
	var *= ( nsp/(nsp-1) );

	double rng = std::sqrt(2*var*std::log(2.0/confidence)/nsp) + 7*std::log(2.0/confidence)/3.0/(nsp-1);

	(*ptr).logBias = std::log(rng) + (*ptr).UB;
	//

	if ( (*ptr).nsp > 1 ) {
		ub = logsumexp( {(*ptr).logBias,  (*ptr).logEx} );
		if ( (*ptr).logEx > (*ptr).logBias) {
				lb = std::log( std::exp( (*ptr).logEx ) - std::exp( (*ptr).logBias ) );
		}
	}
	return std::pair<double,double>(lb,ub);
}
tree<mcts::node>::iterator mcts::expandTree(tree<node>::iterator ptr) {
	// expand a leaf node, return pointer to its OR child
	assert( (*ptr).type == AND_NODE );
	assert( exploredTree.begin(ptr) == exploredTree.end(ptr) );

	int loc = nvar; // default, for the virtual root
	if ( (*ptr).depth > -1 ) {
		loc = priority[ ((*ptr).X).label() ];
	// debug
//	 	assert( ((*ptr).X).label() == order[loc] );
	}
	--loc; // child order
	if (loc < 0) return NULL; // already reach the first variable in the order
	//
	mex::vector<uint32_t> config(nvar);
	auto pointer = ptr;
	while( pointer != exploredTree.begin() ) {
		if ( (*pointer).type == AND_NODE  ) {
			config[ (*pointer).X ] = (*pointer).val;
		}
		pointer = exploredTree.parent(pointer);
	}

	auto X = wmbUB.var(order[loc]);
	// OR node
	node n;
	n.type = OR_NODE;
	n.depth = (*ptr).depth + 1;
	n.X = X;
	n.cost = (*ptr).cost; // no cost on AND->OR arc
	setBounds(n, config); // set deterministic bounds
	ptr = exploredTree.append_child(ptr, n);
	++GSIZE; // count OR node as well
	// AND nodes
	for (int i=0; i<X.states(); ++i) {
		node kid;
		kid.type = AND_NODE;
		kid.depth = n.depth;
		kid.X = n.X;
		kid.val = i;
		config[kid.X] = kid.val;
		kid.cost = n.cost + wmbUB.heuristicTheta(kid.X, config);
		setBounds(kid, config); // set deterministic bounds
		exploredTree.append_child(ptr, kid);
		++GSIZE;
	}
	// update path deterministic bounds
	updatePathBounds(ptr);
	// return the OR child
	return ptr;
}

void mcts::runGBFS(long treeSizeLimit, double timeLimit) {
/*
 * GBFS
 */
	double startIter = mex::timeSystem();
	// use priority queue to implement for convenience, can be improved later on.
	// order is already defined in pair, compare pair.first, and then pair.second
	std::priority_queue<myPair> queue;
	auto rt = exploredTree.begin(); // virtual root
	auto ptr = rt;

	myPair rtPair;
	rtPair.ptr = rt;
//	rtPair.priority = calcPriority(*rt);
	rtPair.priority = (*rt).ub; // only use upper bound

	queue.push( rtPair );
	int maxDepth = 0;

	while( !queue.empty() ) {
		if ( GSIZE > treeSizeLimit ) {
			std::cout << "reach tree size limit: "<< treeSizeLimit <<std::endl;
			std::cout << "max depth reached: " << maxDepth << std::endl;
			break;
		}
		if ( mex::timeSystem()-startIter > timeLimit  ) {
			std::cout.precision(10);
			std::cout << "reach search time limit (sec): "<< timeLimit <<"\n";
			std::cout << "max depth reached: " << maxDepth << std::endl;
			break;
		}
		auto top = queue.top();
		queue.pop();
		//debug
//		std::cout.precision(10);
//		std::cout << "top's priority = " << top.priority <<std::endl;
		// expand the tree and update bounds
		ptr = expandTree(top.ptr);
		// append top node to the tree
		if (ptr != NULL) {
			maxDepth = std::max(maxDepth, (*ptr).depth);
			for ( auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib ) {
				myPair kidPair;
				kidPair.ptr = sib;
//				kidPair.priority =  calcPriority(*sib);
				kidPair.priority = (*sib).ub; // use upper bound to avoid memory use for lower bound
				queue.push( kidPair );
			}
		}
	}
}
void mcts::twoStageSamplingWithFixedTree() {
/*
 * two-stage sampling procedure with a fixed tree
 * for internal nodes, using bound based sampling; for leaf nodes, using wmb mixture proposal to sample
 * caveat: for leaf nodes, bounds are their initial WMB bounds
*/
	mex::vector<uint32_t> config(nvar);
	std::vector<double> logQxVec(nvar+1, 0.0);
	double logFx = 0.0;
	// sample the path from root to leaf
	auto ptr = exploredTree.begin();
	auto ptrPar = ptr;
	int loc = nvar;

	// stage 1:
	while ( exploredTree.begin(ptr) != exploredTree.end(ptr) ) {
		// leaf nodes must be AND nodes
		if ( (*ptr).type == OR_NODE  ) {
			// return (randomly) picked AND child
			ptrPar = ptr;
			ptr = upperBoundBasedSampling(ptr); // return pointer to the AND child
			config[ (*ptr).X ] = (*ptr).val;

			loc = priority[ ((*ptr).X).label() ];
			// logQxVec[j] = log q(x_{j-1} | par(x_{j-1}) ) for initialization, see sampleMixtureConditioned for details
			logQxVec[loc+1] = (*ptr).UB - (*ptrPar).UB;
		}
		else {
			// AND node has at most one child
			ptr = exploredTree.begin(ptr);
		}
	}
	assert(  (*ptr).type == AND_NODE );

	// stage 2: sample from the leaf according to the mixture proposal
	int ID = -1; // default, root id
	if ( (*ptr).depth > -1  ) ID = ((*ptr).X).label();

/*
	//debug
	loc = nvar;
	if (ID > -1) loc = priority[ ID ];
	std::cout << "before:\n";
	for (int i=0; i<=loc ; ++i  ) {
		std::cout << logQxVec[i] << ", ";
	}
	std::cout << " | " ;
	for (int i=loc+1; i<logQxVec.size() ; ++i  ) {
		std::cout << logQxVec[i] << ", ";
	}
*/

	logFx = wmbUB.sampleMixtureConditioned(ID, config, logQxVec); // config, logQxVec, logFx will be updated


/*
  	 // debug
	std::cout << "\nafter:\n";
	for (int i=0; i<=loc ; ++i  ) {
		std::cout << logQxVec[i] << ", ";
	}
	std::cout << " | " ;
	for (int i=loc+1; i<logQxVec.size() ; ++i  ) {
		std::cout << logQxVec[i] << ", ";
	}
*/



	// we only have to update the virtual root
	//	addSampleToPath(ptr, config, logQxVec, logFx);
	addSampleToNode(exploredTree.begin(), logFx - logQxVec[nvar]);

/*
	// debug
	double tmpval = 0.0;
	if (ID>-1) tmpval = logQxVec[ priority[ID] ];
	else tmpval = logQxVec[nvar];
	tmpval += (*ptr).UB-(*(exploredTree.begin())).UB;
	assert ( logQxVec[nvar] > tmpval - 1e-6  &&  logQxVec[nvar] < tmpval + 1e-6);
	std::cout << "\n logFx - logQxVec[nvar] = "  << logFx - logQxVec[nvar]  << std::endl;

	auto rt = exploredTree.begin();
	std::cout << (*rt).nsp << "/"  << nSample <<"\n" ;
	std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << std::endl;
*/
}
void mcts::startTwoStageSamplingWithFixedTree(double searchTimeLimit, double memBudget, double timeBudget) {
	assert( searchTimeLimit <= timeBudget );
	double startIter = mex::timeSystem(); // just for search
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	// in byte
	double memUnit = sizeof(*(exploredTree.begin()));
//	double wmbMem = wmbUB.memory() + wmbLB.memory() ; // estimated size of wmb object in MB
	double wmbMem = wmbUB.memory(); // only use upper bound

//	long memLimit = 1024*4; // no more than 4 GB
	double searchMem = memBudget - wmbMem;
	if (searchMem < 0) {
		std::cout << "Not enough memory allocated for WMB construction, quit!" << std::endl;
		return;
	}

	auto treeSizeLimit = long(1024*1024*searchMem/memUnit);

	std::cout << "memory budget (MB): " << memBudget << "\n";
	std::cout << "memory used for WMB (MB): " << wmbMem << "\n";
	std::cout << "memory used for search (MB): " << searchMem << "\n";
	std::cout << "node size (B): " << memUnit << "\n";
	std::cout << "tree size limit: " << treeSizeLimit <<"\n";

	// run GBFS
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout << "Running GBFS..." << std::endl;
	runGBFS(treeSizeLimit, searchTimeLimit);
	std::cout << "Done GBFS!\n";

	// debug
//	printTree();

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	//
	auto rt = exploredTree.begin();
	std::cout << "GSIZE (including all AND and OR nodes) = " << GSIZE << "\n";
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << (*rt).LB << " < ln Z < " << (*rt).UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;

	cur_nSample = 0;
	std::pair<double, double> prob;
	long outBase = 1;

	// remaining time for sampling
	double startSampling = mex::timeSystem();
	std::cout << "Running two stage sampling..."<<std::endl;
	while ( cur_nSample < nSample &&  mex::timeSystem()-startIter < timeBudget) {
//		int outBase = int(std::max(1.0, nSample/100.0));
		twoStageSamplingWithFixedTree();
		cur_nSample = (*rt).nsp;

		if (cur_nSample == outBase) {
				// calculate EBB at root, using refined bounds
			prob = calcRootEBB();
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
//			std::cout << "GSIZE = " << GSIZE << "\n";
//			std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << (*rt).logEx << " < " << prob.second << "\n";
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

			outBase += outBase;
//				outBase *= 10;
		}
	}
	std::cout << "Done two stage sampling!\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		// calculate EBB at root, using refined bounds
	prob = calcRootEBB();
	cur_nSample = (*rt).nsp;
	if ( cur_nSample < nSample ) std::cout << "reach time budget(sec): " << timeBudget << "\n";
	else  std::cout << "reach sample size limit: " <<  nSample << "\n";
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
	std::cout << "GSIZE (including all AND and OR nodes) = " << GSIZE << "\n";
	std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
	std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << (*rt).logEx << " < " << prob.second << "\n";
	std::cout << "Real time (sec) for search: " << startSampling - startIter << ", for sampling: " << mex::timeSystem() - startSampling << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
}
void mcts::printTree() {
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~printing the explored tree~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	// traverse the tree by DFS
	for (auto iter = exploredTree.begin(); iter != exploredTree.end(); ++iter) {
		auto n = *iter;
		int id = -1;
		if ( n.depth > -1  ) id = (n.X).label();

		std::cout <<"variable "<< id << " with value " << n.val << ", depth = " << n.depth << "\n";
		std::cout << "ub = " << n.ub << ", " << "UB = " << n.UB <<"\n";
		std::cout << "lb = " << n.lb << ", " << "LB = " << n.LB <<"\n";
	}

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}

void mcts::startTwoStageSampling() {
	double startIter = mex::timeSystem(); // just for search
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	// in byte
	long memUnit = sizeof(*(exploredTree.begin()));
	long memLimit = 1024*4; // no more than 4 GB
	auto treeSizeLimit = long(1024*1024*memLimit/memUnit);

	std::cout << "memLimit (MB): " << memLimit << "\n";
	std::cout << "memUnit (B) " << memUnit << "\n";
	std::cout << "treeSizeLimit = " << treeSizeLimit <<"\n";

	cur_nSample = 0;

	// try different way to increase the tree size
	int NSP = 1;
	int multiplier = 10;


	while ( NSP <= nSample ) {
//		int outBase = int(std::max(1.0, nSample/100.0));
		int outBase = 1;
		auto rt = exploredTree.begin();
		std::pair<double, double> prob;

		while( GSIZE < treeSizeLimit  &&  cur_nSample < nSample) {
			twoStageSampling(NSP);
			cur_nSample = (*rt).nsp;

			if (cur_nSample == outBase) {
				prob = calcEBB(rt);

/*
				//debug
				double test_LB = -std::numeric_limits<double>::infinity();
				double test_UB = test_LB;
				std::queue<tree<node>::iterator> queue;
				queue.push(rt);
				while( !queue.empty() ) {
					auto top = queue.front();
					queue.pop();

					if ( exploredTree.begin(top) == exploredTree.end(top)  ) {
						assert( (*top).type == AND_NODE );
						test_LB = logsumexp(  {test_LB, (*top).lb} );
						test_UB = logsumexp(  {test_UB, (*top).ub} );
					}

					for (auto it = exploredTree.begin(top); it != exploredTree.end(top); ++it)  queue.push(it);
				}
				//debug
*/

				std::cout.precision(10);
				std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
				std::cout << "NSP = " << NSP << ", GSIZE = " << GSIZE << "\n";
//				std::cout << "Test bounds: "<< std::setw(10) << test_LB << " < " << (*rt).logEx << " < " <<  test_UB << "\n";
				std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
				std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << (*rt).logEx << " < " << prob.second << "\n";
				std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

//				outBase += outBase;
				outBase *= 10;
			}
		}

		prob = calcEBB(rt);
		cur_nSample = (*rt).nsp;
		if ( cur_nSample >= nSample ) std::cout << "reach sample limit!\n";
		else  std::cout << "reach memory limit!\n";
		std::cout.precision(10);
		std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
		std::cout << "NSP = " << NSP << ", GSIZE = " << GSIZE << "\n";
		std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
		std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << (*rt).logEx << " < " << prob.second << "\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

		//a hack, to make sure the last run would be the original one-stage sampling
		NSP *=  multiplier;
		if (GSIZE > 0 && NSP > nSample) NSP = nSample;
		// delete the whole tree
		exploredTree.clear();
		// re-initialization
		setRootConfig();
		GSIZE = 0;
		cur_nSample = 0;
	}
}

void mcts::start() {
	double startIter = mex::timeSystem(); // just for search
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	int outBase = 1e3;
	int outMultiplier = 10;
	while (cur_nSample < nSample) {
		sampleOne();
		if (cur_nSample == outBase) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

			outBase *= outMultiplier;
			int nsp_thresh =  5;
			int maxDepth = -1;
			while( nsp_thresh <= cur_nSample) {
				maxDepth = aggregateBounds(nsp_thresh, true);
				if (maxDepth == -1) break;
				nsp_thresh *= 2;
			}
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		}
	}
	std::cout << "[" << mex::timeSystem() - startIter << "]: succeeded!\n";
}

void mcts::genTreeBySampling(long treeSizeLimit, double timeLimit, int NSP) {
/*
 * apply upper-bound-based sampling to generate a fixed tree
 * this strategy may or may not generate a lower upper bound tree compared to GBFS
 * However, it's possible that it generates a more balanced tree, i.e.,
 * for two sibling nodes, the ratio between their upper bounds is closer to the ratio between their actual value.
 * Caveat: we should do pruning to avoid that sampling visits some solved branch again and again
 *
 */
	double startIter = mex::timeSystem();
	int maxDepth = -1;
//	int count =0 ;
//	int expcnt = 0;

	while (true) {
		if (GSIZE > treeSizeLimit) {
			std::cout << "reach tree size limit: " << treeSizeLimit
					<< std::endl;
			std::cout << "max depth reached: " << maxDepth << std::endl;
			break;
		}
		if (mex::timeSystem() - startIter > timeLimit) {
			std::cout.precision(10);
			std::cout << "reach search time limit (sec): " << timeLimit << "\n";
			std::cout << "max depth reached: " << maxDepth << std::endl;
			break;
		}

		// sample the path from root to leaf
		auto ptr = exploredTree.begin();
		auto ptrPar = ptr;

		while (exploredTree.begin(ptr) != exploredTree.end(ptr)) {
			// leaf nodes must be AND nodes
			if ((*ptr).type == OR_NODE) {
				// return (randomly) picked AND child
				ptrPar = ptr;
				ptr = upperBoundBasedSampling(ptr); // return pointer to the AND child
			} else {
				// AND node has at most one child
				ptr = exploredTree.begin(ptr);
			}
		}

		assert((*ptr).type == AND_NODE);
		// sufficient to only change nsp of frontier nodes
		++(*ptr).nsp;
//		++count;
		// expand the tree if necessary
		if ((*ptr).nsp > NSP) {
			ptr = expandTree(ptr);
			if (ptr != NULL) {
//				++expcnt;
				maxDepth = std::max(maxDepth, (*ptr).depth);
			}
		}
	}
	// debug
//	std::cout << "count = " << count << std::endl;
//	std::cout << "expcnt = " << expcnt << std::endl;
}
void mcts::startTwoStageSamplingWithFixedTree(double searchTimeLimit, double searchMemLimit, double timeBudget, int NSP) {
	assert( searchTimeLimit <= timeBudget );
	double startIter = mex::timeSystem(); // just for search
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	// in byte
	int memUnit = sizeof(*(exploredTree.begin()));
//	long memLimit = 1024*4; // no more than 4 GB
	double memLimit = searchMemLimit;
	auto treeSizeLimit = long(1024*1024*memLimit/memUnit);

	std::cout << "search memory Limit (MB): " << memLimit << "\n";
	std::cout << "node size (B): " << memUnit << "\n";
	std::cout << "tree size limit: " << treeSizeLimit <<"\n";
	std::cout << "NSP: " << NSP << "\n";

	// run GBFS
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout << "Sampling to generate the tree..." << std::endl;
	genTreeBySampling(treeSizeLimit, searchTimeLimit, NSP);
	std::cout << "Done tree generation!\n";

	// debug
//	printTree();

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	//
	auto rt = exploredTree.begin();
	//
	(*rt).nsp = 0; // re-set to 0 because it has been used in genTreeBySampling
	//
	std::cout << "GSIZE (including all AND and OR nodes) = " << GSIZE << "\n";
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << (*rt).LB << " < ln Z < " << (*rt).UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;

	cur_nSample = 0;
	std::pair<double, double> prob;
	long outBase = 1;

	// remaining time for sampling
	double startSampling = mex::timeSystem();
	std::cout << "Running two stage sampling..."<<std::endl;
	while ( cur_nSample < nSample &&  mex::timeSystem()-startIter < timeBudget) {
//		int outBase = int(std::max(1.0, nSample/100.0));
		twoStageSamplingWithFixedTree();
		cur_nSample = (*rt).nsp;

		if (cur_nSample == outBase) {
				// calculate EBB at root, using refined bounds
			prob = calcRootEBB();
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
//			std::cout << "GSIZE = " << GSIZE << "\n";
//			std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
			std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << (*rt).logEx << " < " << prob.second << "\n";
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

			outBase += outBase;
//				outBase *= 10;
		}
	}
	std::cout << "Done two stage sampling!\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		// calculate EBB at root, using refined bounds
	prob = calcRootEBB();
	cur_nSample = (*rt).nsp;
	if ( cur_nSample < nSample ) std::cout << "reach time budget(sec): " << timeBudget << "\n";
	else  std::cout << "reach sample size limit: " <<  nSample << "\n";
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: cur_nSample/nSample = " << cur_nSample << "/"  << nSample <<"\n" ;
	std::cout << "GSIZE (including all AND and OR nodes) = " << GSIZE << "\n";
	std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
	std::cout << "Empirical Bernstein bounds: " << std::setw(10) << prob.first << " < " << (*rt).logEx << " < " << prob.second << "\n";
	std::cout << "Real time (sec) for tree generation: " << startSampling - startIter << ", for sampling: " << mex::timeSystem() - startSampling << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
}

tree<mcts::node>::iterator mcts::expandTree(tree<node>::iterator ptr, mex::vector<uint32_t>& config) {
	// expand a leaf node, return pointer to its OR child
	assert( (*ptr).type == AND_NODE );
	assert( exploredTree.begin(ptr) == exploredTree.end(ptr) );

	int loc = nvar; // default, for the virtual root
	if ( (*ptr).depth > -1 ) {
		loc = priority[ ((*ptr).X).label() ];
	// debug
//	 	assert( ((*ptr).X).label() == order[loc] );
	}
	--loc; // child order
	if (loc < 0) return NULL; // already reach the first variable in the order
	//
//	mex::vector<uint32_t> config(nvar);
	auto pointer = ptr;
	while( pointer != exploredTree.begin() ) {
		if ( (*pointer).type == AND_NODE  ) {
			config[ (*pointer).X ] = (*pointer).val;
		}
		pointer = exploredTree.parent(pointer);
	}

	auto X = wmbUB.var(order[loc]);
	// OR node
	node n;
	n.type = OR_NODE;
	n.depth = (*ptr).depth + 1;
	n.X = X;
	n.cost = (*ptr).cost; // no cost on AND->OR arc
	setBounds(n, config); // set deterministic bounds
	ptr = exploredTree.append_child(ptr, n);
	++GSIZE; // count OR node as well
	// AND nodes
	for (int i=0; i<X.states(); ++i) {
		node kid;
		kid.type = AND_NODE;
		kid.depth = n.depth;
		kid.X = n.X;
		kid.val = i;
		config[kid.X] = kid.val;
		kid.cost = n.cost + wmbUB.heuristicTheta(kid.X, config);
		setBounds(kid, config); // set deterministic bounds
		exploredTree.append_child(ptr, kid);
		++GSIZE;
	}
	// update path deterministic bounds
	updatePathBounds(ptr);
	// return the OR child
	return ptr;
}

void mcts::runGBFSWithNewGap(long treeSizeLimit, double timeLimit, int NSP) {
	double startIter = mex::timeSystem();
	// use priority queue to implement for convenience, can be improved later on.
	// order is already defined in pair, compare pair.first, and then pair.second
	mex::vector<uint32_t> tuple(nvar);
	std::vector<double> logQxVec(nvar+1, 0.0);

	std::priority_queue<myPair> queue;
	auto rt = exploredTree.begin(); // virtual root
	auto ptr = rt;

	myPair rtPair;
	rtPair.ptr = rt;
	rtPair.priority = calcPriority(*rt);

	queue.push( rtPair );
	int maxDepth = 0;
	long outBase = 1;

	while( !queue.empty() ) {
		if (GSIZE >= outBase) {
			std::cout.precision(10);
			std::cout << "[" << mex::timeSystem() - startIter << "]: GSIZE = "
					<< GSIZE << "\n";
			std::cout << "Deterministic bounds: " << std::setw(10) << (*rt).LB
					<< " < " << (*rt).UB << "\n";
			std::cout
					<< "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			outBase += outBase;
		}


		if ( GSIZE > treeSizeLimit ) {
			std::cout << "reach tree size limit: "<< treeSizeLimit <<std::endl;
			std::cout << "max depth reached: " << maxDepth << std::endl;
			break;
		}
		if ( mex::timeSystem()-startIter > timeLimit  ) {
			std::cout.precision(10);
			std::cout << "reach search time limit (sec): "<< timeLimit <<"\n";
			std::cout << "max depth reached: " << maxDepth << std::endl;
			break;
		}
		auto top = queue.top();
		queue.pop();

		// expand the tree and update bounds
		ptr = expandTree(top.ptr, tuple);
//		//debug
//		std::cout.precision(10);
//		std::cout << "top's priority = " << top.priority << std::endl;
//		std::cout << "root's bounds: " << (*rt).LB <<  " < " << (*rt).UB << std::endl;

		// append top node to the tree
		if (ptr != NULL) {
			maxDepth = std::max(maxDepth, (*ptr).depth);
			for ( auto sib = exploredTree.begin(ptr); sib != exploredTree.end(ptr); ++sib ) {
				// sample NSP times
				tuple[ (*sib).X ] = (*sib).val;
				(*sib).logEx = (*sib).lb; // for later use
				//
				int ID = ((*ptr).X).label();
				int loc = priority[ID];
				double logFx = 0.0;
				// sample NSP times
				for (int nsp = 0; nsp < NSP; ++nsp) {
					// do not need to re-set logQxVec
					logFx = wmbUB.sampleMixtureConditioned(ID, tuple,
							logQxVec); // config, logQxVec, logFx will be updated

					addSampleToNode(sib, logFx - logQxVec[loc]);
				}

				//debug
//				std::cout << "nsp = " << (*sib).nsp << "\n";
//				std::cout << (*sib).lb << " < " << (*sib).logEx << " < " << (*sib).ub << std::endl;
				//
				myPair kidPair;
				kidPair.ptr = sib;
				kidPair.priority =  (*sib).ub + std::log( 1 - std::exp( (*sib).logEx - (*sib).ub ) );
				queue.push( kidPair );
			}
		}
	}
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: GSIZE = "
			<< GSIZE << "\n";
	std::cout << "Deterministic bounds: " << std::setw(10) << (*rt).LB
			<< " < " << (*rt).UB << "\n";
	std::cout
			<< "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void mcts::startGBFSWithNewGap(double searchTimeLimit, double searchMemLimit, double timeBudget, int NSP) {
	assert( searchTimeLimit <= timeBudget );
	double startIter = mex::timeSystem(); // just for search
	std::cout.precision(10);
	std::cout << "[" << mex::timeSystem() - startIter << "]: "
						<< std::setw(10) << LB << " < ln Z < " << UB << "\n";
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	// in byte
	int memUnit = sizeof(*(exploredTree.begin()));
//	long memLimit = 1024*4; // no more than 4 GB
	double memLimit = searchMemLimit;
	auto treeSizeLimit = long(1024*1024*memLimit/memUnit);
//	//
//	treeSizeLimit = 54;


	std::cout << "search memory Limit (MB): " << memLimit << "\n";
	std::cout << "node size (B): " << memUnit << "\n";
	std::cout << "tree size limit: " << treeSizeLimit <<"\n";
	std::cout << "NSP: " << NSP << "\n";

	// run GBFS
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout << "Sampling to boost GBFS..." << std::endl;
	runGBFSWithNewGap(treeSizeLimit, searchTimeLimit, NSP);
	std::cout << "Done search!\n";
//	std::cout << "GSIZE (including all AND and OR nodes) = " << GSIZE << "\n";
//	std::cout << "Deterministic bounds: "<< std::setw(10) << (*rt).LB << " < " << (*rt).logEx << " < " << (*rt).UB << "\n";
}
//EOF


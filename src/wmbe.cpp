// Last revised by qlou on 09//11/15: set flag for pseudoTree in wmbe::init()
#include <assert.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstring>
#include <map>

#include "wmbe.h"

/* Weighted Mini-Bucket Elimination object source code
*/

namespace mex {

// TODO:
//   "run"-like function:  1 pass fwd, then bwd/fwd etc, timeout? always finish any down pass?
//     (estimate timing and quit "nearest" to stoptime?)
//     keep track of & modify damping?
//
//   forward heuristic computation: need to check Lars / old mbe
//
//   reparameterization: write thetas back into gm factors
//   
//   jglp form: fwd & bwd messages use weight 1; change in-bucket behavior (?)
//
//   incremental builds?  use variational scoring, pass messages down & up from b;
//     look at current memory usage per bucket; ramp up slowly with sbound?
//   
//   match construction; test usefulness of big match sets
//     use sets to control better?  how to remove redundant ones?
//   
//   merge clique operation that preserves jg property 
//     add clique recursively until contained; walk back merging

/*
wmbe& wmbe::operator=(wmbe const& gm) {
	graphModel::operator=((graphModel&)gm);		// copy factors, adjacency, etc.
  _order = gm._order; 
	_priority = gm._priority; 
	_varWeight = gm._varWeight;
	_parents = gm._parents;
	_elim = gm._elim;
	_iBound = gm._iBound; _sBound = gm._sBound; 
  _obj = gm._obj; _dampTheta = gm._dampTheta; _dampWeight = gm._dampWeight;
	
	std::map<nodeID,nodeID> nodeMap;
	for (size_t b=0;b<gm._nodes.size();++b) {
		for (size_t i=0;i<gm._nodes[b].size();++i) {
			nodeID n1 = gm._nodes[b][i];
			nodeID n2 = addNodeBasic( gm.node(n1).clique );
			nodeMap[n1] = n2;
		} // !!! need to fix up parent/children ptrs
	}	
	_factorIn = gm._factorIn;
	for (size_t f=0;f<_factorIn.size();++f) _factorIn[f] = nodeMap[_factorIn[f]];

}
*/

/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
//
void wmbe::init(bool isAndOr) {
  _obj = 0.0;
  if (_order.size()==0) {               // if we need to construct an elimination ordering
    double tic=timeSystem();
     _order=graphModel::order(ordMethod);
		 _priority.clear();
     _parents.clear();                   // (new elim order => need new pseudotree) TODO !!! should do together
     std::cout<<"Order in "<<timeSystem()-tic<<" sec\n";
  }
	if (_priority.size()==0) {
		_priority.resize(_order.size());
		for (size_t i=0;i<_order.size();++i) _priority[_order[i]]=i;
	}
  if (_parents.size()==0) {             // if we need to construct a pseudo-tree
    double tic=timeSystem();
  //  _parents=graphModel::pseudoTree(_order);
    _parents=graphModel::pseudoTree(_order, isAndOr); // qlou: allow you to choose AND/OR or OR, 6/16/2016
    std::cout<<"Pseudo in "<<timeSystem()-tic<<" sec\n";
  }

  // skip preprocessing message passing (do in other ways?)

	_nodes.clear(); _nodes.resize(_order.size());
	_match.clear(); _match.resize(_order.size());
	_heuristic.clear(); _heuristic.resize(_order.size());
  _nodeList.clear();
	_factorIn.resize( nFactors() );
	for (size_t f=0;f<factors().size();++f) {
		nodeID n = addNodeBasic( factor(f).vars() );
		node(n).origFactors.push_back(f);
		_factorIn[f] = n;
	}

	_elim.clear(); _elim.resize(nvar());	// !!! ??? TODO???
}

void wmbe::init() {
  _obj = 0.0;
  if (_order.size()==0) {               // if we need to construct an elimination ordering
    double tic=timeSystem();
     _order=graphModel::order(ordMethod);
		 _priority.clear();
     _parents.clear();                   // (new elim order => need new pseudotree) TODO !!! should do together
     std::cout<<"Order in "<<timeSystem()-tic<<" sec\n";
  }
	if (_priority.size()==0) {
		_priority.resize(_order.size());
		for (size_t i=0;i<_order.size();++i) _priority[_order[i]]=i; 
	}
  if (_parents.size()==0) {             // if we need to construct a pseudo-tree
    double tic=timeSystem();
  //  _parents=graphModel::pseudoTree(_order);
/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
    _parents=graphModel::pseudoTree(_order, false); // qlou: false, just OR, 6/16/2016
    std::cout<<"Pseudo in "<<timeSystem()-tic<<" sec\n";
  }

  // skip preprocessing message passing (do in other ways?)

	_nodes.clear(); _nodes.resize(_order.size());
	_match.clear(); _match.resize(_order.size());
	_heuristic.clear(); _heuristic.resize(_order.size());
  _nodeList.clear();
	_factorIn.resize( nFactors() );
	for (size_t f=0;f<factors().size();++f) {
		nodeID n = addNodeBasic( factor(f).vars() );
		node(n).origFactors.push_back(f);
		_factorIn[f] = n;	
	}

	_elim.clear(); _elim.resize(nvar());	// !!! ??? TODO???
}


wmbe::nodeID wmbe::addNodeBasic(VarSet vs) {
  if (vs.nvar() == 0) vs = VarSet(var(_order[0]));		// if no variables, just pick first eliminated (TODO better fix?)
  //assert( vs.nvar() > 0 );
	
  size_t b = _order.size();     // get bucket for this clique (1st elim variable)
  for (size_t i=0; i<vs.size(); ++i) b = std::min(b,_priority[vs[i]]); 
  Var X = var( _order[b] );
	//std::cout<<"Adding "<<vs<<" to bucket "<<b<<"\n";
	for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
		//if (node(*ni).clique >> vs) std::cout<<"Found "<<node(*ni).clique<<"\n";
		if (node(*ni).clique >> vs) return *ni;		// already contained: return the nodeID
	}
	_nodeList.push_back( jgNode() );
	nodeID n = &(_nodeList.back());
	_nodes[b].push_back(n);
	node(n).clique = vs;
	//std::cout<<"   (added)\n";
	return n;
}


double wmbe::score(Var X, nodeID ni, nodeID nj) {
	double s;
	size_t nv1 = node(ni).clique.size(), nv2 = node(nj).clique.size();
	if (nv1 < nv2) std::swap(nv1,nv2);
	size_t iBound = std::max(std::max(_iBound,nv1-1),nv2-1);
	double sBound = std::max(std::max(_sBound,node(ni).clique.nrStatesDouble()),node(nj).clique.nrStatesDouble());
	VarSet both = node(ni).clique + node(nj).clique;
	if (both.nvar() > iBound+1 || both.nrStatesDouble() > sBound) s = -3;		// doesn't fit : -3
	else if (1) s = double(nv1) + double(nv2)/double(nv1);						// scope-based : larger = higher
	else { }
	return s;
}


void wmbe::build() { //size_t iBound, size_t sBound) {
	typedef detail::nodePair nodePair;
	typedef std::pair<double, nodePair> pairPriority;
	std::multimap<double,nodePair > scores;
	typedef std::map<nodePair, std::multimap<double,nodePair >::iterator> reverseType;
	reverseType reverseScore;

	//std::cout<<"Building with iBound="<<_iBound<<"; sBound="<<_sBound<<"\n";
	for (size_t b=0;b<_order.size();++b) {		// for each bucket b,
		//std::cout<<"Merging bucket "<<b<<"\n";
		Var X = Var(_order[b],0);

		//std::cout<<"Starting with nodes ";
		//for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) std::cout<<node(*ni).clique<<" "; 
		//std::cout<<"\n";

		// Score each pair of mini-buckets; mark self merges as -1
		for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
			for (vector<nodeID>::iterator nj=ni+1;nj!=_nodes[b].end();++nj) {
				double s = score( X, *ni, *nj );
				nodePair np(*ni,*nj);
				reverseScore[np] = scores.insert( pairPriority( s , np ) );
			}
			reverseScore[nodePair(*ni,*ni)] = scores.insert( pairPriority( -1 , nodePair(*ni,*ni) ) );
		}
		// Now, pull pairs in order of their priority and merge them:
		for (;;) {
			std::multimap<double,nodePair>::reverse_iterator top = scores.rbegin();
      if (top->first < 0) break;                            // if can't do any more, quit
      else {
        nodeID ni=top->second.second, nj=top->second.first;
				//std::cout<<"   Merge "<<node(ni).clique<<"+"<<node(nj).clique<<"="<<node(ni).clique+node(nj).clique<<"\n";
				node(ni).clique += node(nj).clique;
				node(ni).theta  += node(nj).theta;
				node(ni).weight += node(nj).weight;
				// TODO fix up parents / mark as invalid / empty?
				node(ni).children.insert( node(ni).children.end(), node(nj).children.begin(), node(nj).children.end() );
				// TODO compute forward message node(ni)?  TODO: set???
				for (set<nodeID>::iterator ci=node(ni).children.begin();ci!=node(ni).children.end();++ci)
					node(*ci).parent = ni;
				node(ni).origFactors.insert( node(ni).origFactors.end(), node(nj).origFactors.begin(), node(nj).origFactors.end() );
				for (vector<findex>::iterator fi=node(ni).origFactors.begin();fi!=node(ni).origFactors.end();++fi)
					_factorIn[*fi] = ni;
				// TODO: need to erase node *nj from _nodes[b] (currently deferred, see next)
				node(nj).clear();
				for (vector<nodeID>::iterator nk=_nodes[b].begin();nk!=_nodes[b].end();++nk) {
					reverseType::iterator ri = reverseScore.find( nodePair(nj,*nk) );
					if (ri != reverseScore.end()) { scores.erase( ri->second ); reverseScore.erase(ri); }
					if (node(*nk).clique.size() && *nk != ni && *nk != nj) {				// TODO: not a great way to do this...
						double s = score( X, ni, *nk );
						nodePair np(ni,*nk);
						scores.erase( reverseScore[np] );
						reverseScore[np] = scores.insert( pairPriority(s,np) );
					}
				}
			}
		}

		//std::cout<<"Stopping with nodes ";
		//for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) std::cout<<node(*ni).clique<<" "; 
		//std::cout<<"\n";

		// Now pass through _nodes[b] and erase blank, add node for parent to others
		for (vector<nodeID>::reverse_iterator ni=_nodes[b].rbegin();ni!=_nodes[b].rend();++ni) {
			//std::cout<<"  Processing clique "<<node(*ni).clique<<" : ";
      if (node(*ni).clique.size() == 0) {
				std::swap(*ni,_nodes[b].back()); _nodes[b].pop_back();
				//std::cout<<"(blank; removed)\n";
			} else {
				if (node(*ni).clique.size() != 1) {
					node(*ni).parent = addNodeBasic( node(*ni).clique - X );
					node(node(*ni).parent).children.push_back(*ni);
				}
				for(size_t k=_parents[X]; k!=vindex(-1); k=_parents[k]) {	// walk up pseudotree (stop at root)
					_heuristic[_priority[k]].push_back(*ni);					// all these buckets have msgFwd pass through them
					if (node(*ni).clique.contains(var(k))) break;			// if we've arrived at parent bucket, we can stop
				}
			}
    }
		//std::cout<<"Finishing with nodes ";
		//for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) std::cout<<node(*ni).clique<<" "; 
		//std::cout<<"\n";

		// Add matching sets: for now, just "all" (TODO)
		//_match[b].push_back( vector<nodeID>(_nodes[b].begin(),_nodes[b].end()) );
		_match[b].push_back( _nodes[b] );		// OK to use reference to nodes list (? TODO)
		// All done with this bucket; move on to the next
		//std::cout<<"\n";
	}
}

/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
// Qi: 1/6/18
void wmbe::setMsgFwdRep(const int K) {
/*
 * Designed for the augmented graphical model
 * Mark those mini-buckets (nodes) if it is a SUM variable and its forward messages go to a MAX variable
 * Also, replicate those factors that only involve MAX variables
 * Caveat: run after "setTheta" but before "msgFowardAug". we assume NO variable change in each mini-bucket (node) after "build"
 */
	_K = K; // assign once
	// top-down, reverse ordering
	for (int b = _nodes.size() - 1; b >= 0; --b) {
		// Caveat: do not use "size_t b", it converts "-1" to a lager positive number when compared with "0"
		Var X = var(_order[b]);
		bool isMaxVar = _maxVars.contains(X);
		for (nodeList::iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
			node(*ni).vid = X.label();
			// if a Max variable, its parent cannot be a SUM variable
			if (isMaxVar) {
				// replicate its thetas.
				// note that a theta merely involves MAX variables
				// iff it is in a MAX variable's mini-bucket due to the constrained elimination ordering
				node(*ni).theta *= K;
				continue;
			}
			// if this mini-bucket has no parent, set to be true
			// for a SUM variable, set "isMsgFwdRep" to be "true" if 1) no parent; 2) parent is a Max variable.
			auto par = node(*ni).parent;
			if (par == NULL) {
				node(*ni).isMsgFwdRep = true;
			} else {
				// debug, remove it
				assert(node(par).vid >=0 && node(par).vid < nvar());
				//
				Var Y = var(node(par).vid);
				node(*ni).isMsgFwdRep = _maxVars.contains(Y);
			}
		}
	}
}

// Use setTheta after init() and build() to write the model factors into the mini-buckets
void wmbe::setTheta() {
	for (size_t f=0;f<_factorIn.size();++f)
		node(_factorIn[f]).theta += log( factor(f) );
}

// Use reparameterize() after setTheta() and inference to replace the model factors with the mini-bucket versions
void wmbe::reparameterize() {
	clearFactors();
	double logZ = 0.0;
	for (size_t b=0;b<_order.size();++b) {    // for each bucket b,
    for (vector<nodeID>::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
			Factor tmp = node(*ni).theta + 0.5 * node(*ni).msgBwd;
			for (vector<nodeID>::const_iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
				tmp += 0.5 * node(*ci).msgFwd;
				tmp -= 0.5 * node(*ci).msgBwd;
			}	
			if (node(*ni).parent != NULL) tmp -= 0.5 * node(*ni).msgFwd;
			double mx = tmp.max();
			logZ += mx;
			tmp  -= mx;
			addFactor( tmp.exp() );
    }
  }
	logZ /= nFactors();  
	double Z = std::exp(logZ);
	for (size_t f=0;f<_factors.size();++f) _factors[f] *= Z;
}


void wmbe::setElimType(Var v, ElimType elimType) {
	// TODO: check if already set to same value & don't change weights?
	_elim[v] = elimType;
	bucketID b = _priority[v];
	vector<nodeID>::iterator ni=_nodes[b].begin();
	double wt = 0.0;
	size_t nNodes = _nodes[b].size();
	switch (elimType) {
		case ElimType::MaxUpper:	wt = 1e-6; break;
		case ElimType::SumUpper:	wt = 1.0 / nNodes;  break;
		case ElimType::SumBethe:	wt = 1.0; break;	// !!! TODO
		case ElimType::SumLower:	wt = 1.0; if (nNodes > 1) { node(*ni).weight = 2.0; ++ni; wt = -1.0/(nNodes-1); } break;
/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
/*		// add new type for MMAP task
		case ElimType::MaxSumUpper:
			if (_maxVars.contains(v)) {
//				setElimType(v, ElimType::MaxUpper);
				// MaxUpper
				wt = 1e-6;
			} else {
//				setElimType(v, ElimType::SumUpper);
				// SumUpper
				wt = 1.0 / nNodes;
			}
			break;*/
	}
	// TODO: more involved ways of setting the weights?
	for (;ni!=_nodes[b].end();++ni) node(*ni).weight = wt;
}


void wmbe::setWeightPositiveFirst(bucketID b, double tot) {
	size_t nNodes = _nodes[b].size();
  node(*_nodes[b].begin()).weight = tot;
	for (vector<nodeID>::iterator ni=++_nodes[b].begin();ni!=_nodes[b].end();++ni)
		node(*ni).weight = 1e-6;
}

void wmbe::setWeightNegativeFirst(bucketID b, double tot) {
	size_t nNodes = _nodes[b].size();
  node(*_nodes[b].begin()).weight = tot;
	for (vector<nodeID>::iterator ni=++_nodes[b].begin();ni!=_nodes[b].end();++ni)
		node(*ni).weight = -1e-6;
}

void wmbe::setWeightPositiveUniform(bucketID b, double tot) {
	size_t nNodes = _nodes[b].size();
	for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni)
		node(*ni).weight = tot / nNodes;
}

void wmbe::setWeightNegativeUniform(bucketID b, double tot, double pos) { // pos default = -1
	if (pos < 0) pos = 2.0*tot;
	size_t nNodes = _nodes[b].size()-1;
	node(_nodes[b][0]).weight = pos;
	for (vector<nodeID>::iterator ni=++_nodes[b].begin();ni!=_nodes[b].end();++ni)
		node(*ni).weight = (tot - pos) / nNodes;
}


void wmbe::dump(std::ostream& os) {
  for (size_t b=0;b<_order.size();++b) {    // for each bucket b,
		os<<"Bucket "<<b<<" (X"<<_order[b]<<"): ";
		for (vector<nodeID>::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
			os<<"["<<*ni<<"] "<<node(*ni).clique<<" ("<<node(*ni).parent<<"), ";
			//os<<"\n   th:"<<node(*ni).theta<<"\n"<<"   mF:"<<node(*ni).msgFwd<<"\n";
		}
		os<<"\n";
	}
}



Factor wmbe::computeNodeBelief(nodeID n) const {
	//std::cout<<" th:"<<node(n).theta<<"\n"<<"  bw:"<<node(n).msgBwd<<"\n";
  Factor bel = node(n).theta + node(n).msgBwd;   // combine theta, bwd, & fwd messages
  for (vector<nodeID>::const_iterator ci=node(n).children.begin();ci!=node(n).children.end();++ci) {
		//std::cout<<"  fw:"<<node(*ci).msgFwd<<"\n";
		bel += node(*ci).msgFwd;
	}
  bel *= (1.0 / node(n).weight);         // take power
  bel -= bel.logsumexp();                // and normalize
	//std::cout<<"  ("<<node(n).weight<<") => "<<bel<<"\n";
  return bel;
}


void wmbe::updateTheta(bucketID b, double damp, const matchList& match, bool updateBeliefs) {
  if (match.size() == 0) return;   	//  (skip if no nodes to match on)
	size_t nNodes = _nodes[b].size();	// total number of nodes in bucket b
  double wTot = 0.0;								// total weight for this match
  VarSet vAll = node(match[0]).clique; 	// find variables to match on
  for (matchList::const_iterator c=match.begin();c!=match.end();++c) {
    wTot += node(*c).weight;
    vAll &= node(*c).clique;
  }

  // Find marginal beliefs, and average marginal belief:
  vector<Factor> delta(nNodes);
  Factor bavg(1.0); bavg.log();					// TODO SPARSE
	size_t j=0;
  for (matchList::const_iterator c=match.begin();c!=match.end();++c,++j) {
    delta[j] = node(*c).belief.logsumexp( node(*c).clique - vAll ); // compute marginal first
    bavg += delta[j] * (node(*c).weight / wTot);             // & weighted average
  }

  // Find desired change in marginal belief, and update belief & theta (with damping)
	j=0;
  for (matchList::const_iterator c=match.begin();c!=match.end();++c,++j) {
		delta[j] = (bavg - delta[j]);
		//( (delta[j]*=-1) += bavg ) *= damp;		// delta = (bavg - delta)*damp
    // //delta[j] *= -1.0;   // !!! TODO cleaner?
    // //delta[j] += bavg;
    // //delta[j] *= damp;          // delta = (bavg - delta)*damp
    // //(bel[j] += delta[j] ) -= bel[j].logsumexp();	// bel = (bel+delta)-lnZ
		node(*c).belief += delta[j]*damp;
    node(*c).belief -= node(*c).belief.logsumexp(); 
		//if (updateBeliefs) {
    //	node(*c).belief += delta[j];
    //	node(*c).belief -= node(*c).belief.logsumexp(); // bel = bel + delta;
		//}
    node(*c).theta += delta[j]*(node(*c).weight*damp); 
    //delta[j] *= node(*c).weight;
    //node(*c).theta += delta[j];         // theta = theta + delta*weight
  }
}


void wmbe::updateWeights(bucketID b, double step) {
	Var X = var(_order[b]);
	size_t nNodes = _nodes[b].size();
  // Nothing to do for "max" (wi = 0) or "bethe" (wi = 1.0), or single buckets
	if (nNodes <= 1 || _elim[X] == ElimType::MaxUpper || _elim[X] == ElimType::SumBethe) return;

  double Htarget = 0.0;
  vector<double> H(nNodes); // conditional entropies of cliques
	nodeList::iterator nPos;
	size_t j = 0;
	for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni,++j) {
		Factor marg = node(*ni).belief.logsumexp(X);
		H[j] = node(*ni).belief.logEntropy() - marg.logEntropy();
		if (_elim[X] == ElimType::SumUpper) Htarget += node(*ni).weight * H[j];	// average entropy for sum+
		else if (node(*ni).weight > 0.0) { Htarget = H[j]; nPos=ni; }						// or find positive entropy for sum-
	}
	double wTotal = 0.0;
	j = 0;		// Take an (exponential) gradient step
	for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni,++j) {
		if (_elim[X] == ElimType::SumUpper) {		// decrease the upper bound
			node(*ni).weight *= std::exp( step * node(*ni).weight * (H[j]-Htarget) );
			wTotal += node(*ni).weight;
		} else if (node(*ni).weight < 0.0) {		// or increase the lower bound
			node(*ni).weight *= std::exp( -step * node(*ni).weight * (H[j]-Htarget) );
			wTotal += node(*ni).weight;
		}
  }					// finally renormalize positive weights (set positive weight for lower bound)
	if (_elim[X] == ElimType::SumUpper) {
		for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni,++j) 
			node(*ni).weight /= wTotal;
  } else {
		node(*nPos).weight = 1.0 - wTotal;
	}
}


// TODO: what about objective? 
// TODO: make more memory efficient; track memory more exactly (see memory())
double wmbe::msgForward(bucketID b, double dampTheta, double stepWeights) {
  double obj = 0.0;
  Var X = var(_order[b]); 
  size_t nNodes = _nodes[b].size();    // !!! "  "
  if (nNodes > 1) {                    // If more than one mini-bucket:
    if (dampTheta > 0.0 || stepWeights > 0.0)  // compute beliefs if needed for matching
			for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni)
				node(*ni).belief = computeNodeBelief(*ni);

    if (dampTheta > 0.0) 
      for (size_t m=0; m<_match[b].size(); ++m) // count over list of matches
        updateTheta(b, dampTheta, _match[b][m], true);	// update, maintain beliefs

    if (stepWeights > 0.0)
      updateWeights(b, stepWeights);
  } // end: matching updates if more than one mini-bucket

  // Compute forward messages: 
  for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
		node(*ni).belief = Factor(0.0);				// clear to save memory
    Factor bel = node(*ni).theta;					// and build "forward" belief
    for (nodeList::iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
			bel += node(*ci).msgFwd;
		}
		// bel *= 1.0/node(*ni).weight;
		// node(*ni).msgFwd = bel.logsumexp(X)*node(*ni).weight;	// use "into" form?
		// if (_calcBeliefs) bel += msgBwd * 1/wt??; **marginalize into variables, factors? **
    node(*ni).msgFwd = (bel*(1.0/node(*ni).weight)).logsumexp(X)*node(*ni).weight; // power-lse 
    // //_msgFwd[j] = bel[j].logsumexpPower(X, 1.0/_weight[n]);  // take power-lse TODO
    if (node(*ni).parent == NULL) obj += node(*ni).msgFwd[0]; // add roots to overall bound
  } // end: forward messages
	return obj;
}

/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
// Qi: 12/2/2017
// forward message passing for an augmented graphical model
// that is augmented from the original one by adding multiple (say, K) identical copies of SUM variables.
// particularly used for "mmapIS".
// every message sent from a SUM variable to a MAX variable will be raised to the K-th power;
// every factor only involves MAX variables should be raised to the K-th power as well.
// Caveat: this may not be compatible with msgBackward. TODO: create a function "msgBackwardAug"?
// In principle, we still treat this as a summation problem instead of MMAP -- easy for implementation
double wmbe::msgForwardAug(bucketID b, double dampTheta, double stepWeights, const int K) {
/*
 * K: no. of copies of SUM variables
 */
	Var X = var(_order[b]);
	bool isMaxVar  = _maxVars.contains(X);
	double obj = 0.0;
//	Var X = var(_order[b]);
	size_t nNodes = _nodes[b].size();    // !!! "  "
// Qi: matching CAN be removed when debugging
	if (nNodes > 1) {                    // If more than one mini-bucket:
		if (dampTheta > 0.0 || stepWeights > 0.0) // compute beliefs if needed for matching
			for (nodeList::iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni)
				node(*ni).belief = computeNodeBelief(*ni);

		if (dampTheta > 0.0)
			for (size_t m = 0; m < _match[b].size(); ++m) // count over list of matches
				updateTheta(b, dampTheta, _match[b][m], true);// update, maintain beliefs

		if (stepWeights > 0.0)
			updateWeights(b, stepWeights);
	} // end: matching updates if more than one mini-bucket

	// Compute forward messages:
	for (nodeList::iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
		node(*ni).belief = Factor(0.0);				// clear to save memory
		Factor bel = node(*ni).theta;			// and build "forward" belief

		// debug only
		assert(node(*ni).vid == X.label());

		for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
			// whether message from this child has to be raised to the K-th power
			if (node(*ci).isMsgFwdRep){
				bel += K * node(*ci).msgFwd; // message from a Sum variable to a Max variable
			} else {
				bel += node(*ci).msgFwd;
			}
			// debug
			if (node(*ci).isMsgFwdRep) {
				Var Y = var(node(*ci).vid);
				assert( isMaxVar &&  !_maxVars.contains(Y));
			} else {
				Var Y = var(node(*ci).vid);
				assert( !isMaxVar ||  _maxVars.contains(Y));
			}
		}
		// bel *= 1.0/node(*ni).weight;
		// node(*ni).msgFwd = bel.logsumexp(X)*node(*ni).weight;	// use "into" form?
		// if (_calcBeliefs) bel += msgBwd * 1/wt??; **marginalize into variables, factors? **
		node(*ni).msgFwd = (bel * (1.0 / node(*ni).weight)).logsumexp(X) * node(*ni).weight; // power-lse
		// //_msgFwd[j] = bel[j].logsumexpPower(X, 1.0/_weight[n]);  // take power-lse TODO
		if (node(*ni).parent == NULL) {
			// to compute the total bound
			if (node(*ni).isMsgFwdRep) {
				obj += K * node(*ni).msgFwd[0];
				// debug
				assert(!isMaxVar);
			} else {
				obj += node(*ni).msgFwd[0]; // add roots to overall bound
				// debug
				assert(isMaxVar);
			}
		}
	} // end: forward messages
	return obj;
}

void wmbe::msgBackward(bucketID b, double dampTheta, double stepWeights) {
	size_t nNodes = _nodes[b].size();
  // TODO: consider matching steps here?
 
  // Compute backward message from n -> c
  for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
    node(*ni).belief = computeNodeBelief(*ni);	// TODO: don't normalize here?
    node(*ni).belief -= node(*ni).belief.max();
    for (nodeList::const_iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
      VarSet vElim = node(*ni).clique - node(*ci).clique;
      node(*ci).msgBwd  = node(*ni).belief.logsumexp(vElim);
			node(*ci).msgBwd *= node(*ci).weight; 
			node(*ci).msgBwd -= node(*ci).msgFwd;
    } // end: loop over children of n
  }   // end: loop over nodes n in bucket i
}


double wmbe::memory() {
  double mem = 0.0, maxTemp = 0.0;
  for (size_t b=0; b<_nodes.size();++b) {
    double temp = 0.0;
    for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
      mem += node(*ni).clique.nrStates() * sizeof(double);                            // parameters (theta)
      mem += (node(*ni).clique - var(_order[b])).nrStates() * sizeof(double) * 2.0;   // msgFwd + msgBwd (?)
      temp += node(*ni).clique.nrStates() * sizeof(double);                           // beliefs for matching
    }
    maxTemp = std::max(maxTemp, temp);
  }
  return (mem + maxTemp)/1000000;
}


/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
unsigned wmbe::getNumberOfMiniBuckets(Var v) const {
	// get the number of mini-buckets for a variable
	return _nodes[_priority[v]].size();
}
// TODO what's the right way to do this?
double wmbe::eliminateFromNow(std::set<Var>& frontier, std::set<Var>& subtree) const {
	// perform exact elimination from a given set of leaf nodes (current frontier)
	// this helps to upper bound the gain from search reaching a given depth (or explored subgraph/subtree), bound of the oracle
	// note: this is a naive implementation, not suitable for large models
	// vars contains the current frontier
	// double obj = 0.0;
	// create a huge factor, naive
	Factor bulk(0.0);
	// compute thetas for each node in the given sub(andor)tree
	for (auto iter = subtree.cbegin(); iter != subtree.cend(); ++iter) {
		Var v = *iter;
		size_t b = _priority[v];	// find bucket to examine
		// double heur = 0.0;
		for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
			// heur += node(*ni).theta.logsumexp();
			bulk += node(*ni).theta;
		}
		// obj += heur;
	}
	// compute messages passing by/into buckets of each node on the frontier
	for (auto iter = frontier.cbegin(); iter != frontier.cend(); ++iter) {
		Var v = *iter;
		size_t b = _priority[v];	// find bucket to examine
		// double heur = 0.0;
		for (nodeList::const_iterator ci=_heuristic[b].begin();ci!=_heuristic[b].end();++ci) {
			// heur += node(*ci).msgFwd.logsumexp();
			bulk += node(*ci).msgFwd;
		}
		// obj += heur;
	}
	// return obj;
	return bulk.logsumexp();
}

/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Factor> wmbe::getCurrentModel(std::set<Var>& frontier, std::set<Var>& subtree) const {
	// TODO ATI check Qi's code; add delta functions for assigned variables, nothing (?) for unassigned missing (?)
	vector<Factor> fs;
	fs.clear();
	for (auto iter = subtree.cbegin(); iter != subtree.cend(); ++iter) {
		Var v = *iter;
		size_t b = _priority[v];	// find bucket to examine
		// double heur = 0.0;
		for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
			// heur += node(*ni).theta.logsumexp();
			fs.push_back(node(*ni).theta);
		}
	}
	for (auto iter = frontier.cbegin(); iter != frontier.cend(); ++iter) {
		Var v = *iter;
		size_t b = _priority[v];	// find bucket to examine
		// double heur = 0.0;
		for (nodeList::const_iterator ci=_heuristic[b].begin();ci!=_heuristic[b].end();++ci) {
			// heur += node(*ci).msgFwd.logsumexp();
			fs.push_back(node(*ci).msgFwd);
		}
		// obj += heur;
	}
	// add factors for those variables not included to make the model consistent with the original one
	// it might avoid some inconsistency issue.
	for (auto iter = _order.cbegin(); iter != _order.cend(); ++iter) {
		Var v = var(*iter);
		if (subtree.find(v) == subtree.end()) {
			// if not in subtree
			Factor dummy(v,0.0);
			fs.push_back(dummy);
		}
	}

	return fs;
}

/*
mex::vector<uint32_t> wmbe::maxSequential() {	// don't need order or isLog for wmb
  vector<uint32_t> vals( nvar() );
  VarSet done;
  for (int b=_nodes.size()-1;b>=0;--b) {
    Factor bel(0.0);
		for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
      VarSet condv = node(*ni).theta.vars() & done;
			size_t condi = 0; if (condv.size()) condi = sub2ind(condv, vals);
      bel += node(*ni).theta.condition( condv , condi );
      for (nodeList::iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
				size_t condi = 0; if (condv.size()) condi = sub2ind(condv, vals);
        bel += node(*ci).msgFwd.condition( condv, condi );
      }
    }
		assert( bel.nvar() == 1 );
		vals[ bel.vars()[0] ] = bel.argmax();		// TODO: randomize a bit
		done |= bel.vars();
  }
	return vals;
}

mex::vector<uint32_t> wmbe::sample() {	// don't need order or isLog for wmb
  vector<uint32_t> vals( nvar() );
  VarSet done;
  for (int b=_nodes.size()-1;b>=0;--b) {
    Factor bel(0.0);
		for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
      VarSet condv = node(*ni).theta.vars() & done;
			size_t condi = 0; if (condv.size()) condi = sub2ind(condv, vals);
      bel += node(*ni).theta.condition( condv , condi );
      for (nodeList::iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
				size_t condi = 0; if (condv.size()) condi = sub2ind(condv, vals);
        bel += node(*ci).msgFwd.condition( condv, condi );
      }
    }
		assert( bel.nvar() == 1 );
		vals[ bel.vars()[0] ] = bel.exp().sample();		// TODO: ?
		done |= bel.vars();
  }
	return vals;
}
*/

/*




void wmbe::reparameterize() {
  // replace factors in gm object with thetas?  include msgFwds?
  // update _obj ? discard current mbe?  ???
}

*/


/*

  // Scoring function for bucket aggregation
  //   Unable to combine => -3; Scope-only => 1.0; otherwise a positive double score
  double score(const vector<Factor>& fin, const Var& VX, size_t i, size_t j, const vector<Factor>& tmp) {
    double err;
    const Factor& F1=fin[i], &F2=fin[j];                          // (useful shorthand)
    size_t iBound = std::max(std::max(_iBound,F1.nvar()-1),F2.nvar()-1);      // always OK to keep same size
    size_t sBound = std::max(std::max(_sBound,F1.nrStates()),F2.nrStates());
    VarSet both = F1.vars()+F2.vars();
    if ( both.nvar()>iBound+1 || both.nrStates()>sBound ) err=-3;  // too large => -3
    //else if (_byScope) err = 1.0/std::log(both.nrStates()+1);          // greedy scope-based (check if useful???)
    else if (_byScope) err = 1.0/(F1.nvar()+F2.nvar());              // greedy scope-based 2 (check if useful???)
    //else if (_byScope) err = 1;                                    // scope-based => constant score
    else {
      Factor F12 = elim(F1*F2,VX), F1F2=tmp[i]*tmp[j];       // otherwise, compute score
      err = F12.distance(F1F2 , distMethod );         // !!! TODO: FIX: factors are log => +?
    }
    return err;
  }


*/


}

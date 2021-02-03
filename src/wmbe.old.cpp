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


/* Old init
(1) build flist with log-factors
(2) associate "vin" to keep track of variables -> factors
(3) keep track of origination info (factor or which bucket)
(4) Loop: for each elim variable (if we have some factors)
  (a) get list of factors with this variable = ids
  (b) multimap< double, (idx,idx) > with scores
  (c) for each pair of factors in ids, score pair & save reverse pointer ; self pairs are -1
  (d) 


TODO: check for factors without any variables?
*/



/* Add clique:
 (1) find 1st bucket for the clique (earliest priority)
 (2) check if any minibuckets contain the clique; if so; stop & return bucket #
 (3) otherwise, add minibucket:
    (a) get empty minibucket, or push-back additional node
    (b) set internals as desired (mainly clique; weight=?)
    (c) set parent to (recursive) addClique( clique - X )
 (4) check if other minibuckets are contained in this clique; if so
    (a) merge with this one: remove from old parents; add children to this one; clear & push to empty
 (5) push back node # onto list for this bucket
 (6) update match structures?

*/

/*
void wbme::setProperties(std::string opt=std::string()) {
  if (opt.length()==0) {
    setProperties("iBound=4,sBound=inf,Order=MinWidth,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=1");
    _byScope = true;
    return;
  }
  std::vector<std::string> strs = mex::split(opt,',');
  for (size_t i=0;i<strs.size();++i) {
    std::vector<std::string> asgn = mex::split(strs[i],'=');
    switch( Property(asgn[0].c_str()) ) {
      //case Property::ElimOp:  elimOp  = ElimOp(asgn[1].c_str()); break;
      case Property::iBound:  setIBound( atol(asgn[1].c_str()) ); break;
      case Property::sBound:  setSBound( atol(asgn[1].c_str()) ); break;
      //case Property::Order:   _order=_gmo.order(graphModel::OrderMethod(asgn[1].c_str())); break;
      case Property::Order:   _order.clear(); _parents.clear();
                              ordMethod=graphModel::OrderMethod(asgn[1].c_str()); break;
      case Property::Distance: distMethod = Factor::Distance(asgn[1].c_str()); _byScope=false; break;
      case Property::DoMplp:  _doMplp  = atol(asgn[1].c_str()); break;
      case Property::DoMatch: _doMatch = atol(asgn[1].c_str()); break;
      case Property::DoFill:  _doFill = atol(asgn[1].c_str()); break;
      case Property::DoJG:    _doJG = atol(asgn[1].c_str()); break;
      case Property::DoHeur:  _doHeur = atol(asgn[1].c_str()); break;
      default: break;
    }
  }
}
*/



void wmbe::init() {
  _obj = 0.0;
  if (_order.size()==0) {               // if we need to construct an elimination ordering
    double tic=timeSystem();
     _order=_gmo.order(ordMethod);
     _parents.clear();                   // (new elim order => need new pseudotree) TODO !!! should do together
     std::cout<<"Order in "<<timeSystem()-tic<<" sec\n";
  }
  if (_parents.size()==0) {             // if we need to construct a pseudo-tree
    double tic=timeSystem();
    _parents=_gmo.pseudoTree(_order);
    std::cout<<"Pseudo in "<<timeSystem()-tic<<" sec\n";
  }
  // TODO: skip preprocessing message passing (do in other ways?)

  
  // TODO !!! TODO !!! TODO !!! TODO !!! TODO !!! TODO !!! TODO !!! TODO !!! TODO !!!  ////////
}



wmbe::nodeID wmbe::helperAddCliqueBasic(VarSet clique) {
	assert( clique.nvar() > 0 );
	size_t i = _order.size();			// get (priority of) 1st eliminated variable
	for (size_t v=0; v<clique.size(); ++v) i = std::min(i,_priority[clique[v]]);
	size_t x = _order[i];
	for (size_t nn=0;nn<_nodes[x].size();++nn)	// check if it's already contained:
		if (_clique[_nodes[x][nn]] >= clique) return _nodes[x][nn];	
	size_t newNode;
	if (_nodesVacant.size()) {		// if we already have one available just use it
 		newNode = _nodesVacant.back();
		_nodesVacant.pop_back();
	} else {											// otherwise need to increase storage a bit
		newNode = _clique.size();
		_clique.push_back(VarSet());
		_theta.push_back(Factor(0.0));
		_parent.push_back(-1);
		_children.push_back(vector<findex>());
		_weight.push_back(0.0);
		_msgFwd.push_back(Factor());
		_msgBwd.push_back(Factor());
	}
	_clique[newNode] = clique;		// Fill in blank node info and set clique 
	_parent[newNode] = 0;
	_children[newNode].clear();
	_weight[newNode] = 0.0;
	_theta[newNode] = _msgFwd[newNode] = _msgBwd[newNode] = Factor(0.0);
	return newNode;
}

void wmbe::helperAddCliqueParent(Var x, nodeID n) {
	if (_clique[n].size() > 1 && _parent[n] == 0) {
		nodeID par = helperAddCliqueBasic( _clique[n] - x );
		_parent[n] = par;
	}
	if (_theta[n].numel() > 1) {
		// or if incoming messages non-empty???
  	// set downward message if non-empty? 
	}
  // add clique[n] - xi using addBasic
  // set parent of node to returned value
}

wmbe::nodeID wmbe::helperAddClique(VarSet clique) {
  // addCliqueBasic( newclique )
  // recurse on clique - xi
  // check for subsumed cliques and merge with this one
}





Factor wmbe::computeNodeBelief(nodeID n) const {
  Factor bel = _theta[n] + _msgBwd[n];   // combine theta, bwd, & fwd messages
  for (vector<findex>::const_iterator c=_children[n].begin();c!=_children[n].end();++c) bel+=_msgFwd[*c];
  bel *= (1.0/_weight[n]);               // take power
  bel -= bel.logsumexp();                // and normalize
  return bel;
}

void wmbe::updateTheta(bucketID i, double damp, const vector<findex>& match, vector<Factor>& bel) {
  if (match.size() == 0) return;   //  (skip if empty!)
	typedef vector<findex> matchList;
	size_t nNodes = _nodes[i].size();
  double wTot = 0.0;
  VarSet vAll = _clique[match[0]]; 
  for (matchList::const_iterator c=match.begin();c!=match.end();++c) {
    wTot += _weight[*c];
    vAll &= _clique[*c];
  }
  std::map<size_t,size_t> lookup;  // find reverse index: nodes[j]=n <=> lookup[n]=j 
  for (size_t j=0; j<_nodes[i].size(); ++j) lookup[ _nodes[i][j] ] = j; 

  // Find marginal beliefs, and average marginal belief:
  vector<Factor> delta(_nodes[i].size());  // _nodes[i].size()
  Factor bavg(0.0);
  for (matchList::const_iterator c=match.begin();c!=match.end();++c) {
    size_t j = lookup[*c]; //!!!
    delta[j] = bel[j].logsumexp( bel[j].vars() - vAll ); // compute marginal first
    bavg += delta[j] * (_weight[*c] / wTot);             // & weighted average
  }
  // Find desired change in marginal belief, and update belief & theta (with damping)
  for (matchList::const_iterator c=match.begin();c!=match.end();++c) {
    size_t j = lookup[*c]; //!!!
		( (delta[j]*=-1) += bavg ) *= damp;		// delta = (bavg - delta)*damp
    //delta[j] *= -1.0;   // !!! TODO cleaner?
    //delta[j] += bavg;
    //delta[j] *= damp;          // delta = (bavg - delta)*damp
    //(bel[j] += delta[j] ) -= bel[j].logsumexp();	// bel = (bel+delta)-lnZ
    bel[j]   += delta[j];
    bel[j]   -= bel[j].logsumexp(); // bel = bel + delta;
    delta[j] *= _weight[*c];
    _theta[*c] += delta[j];         // theta = theta + delta*weight
  }
}

void wmbe::updateWeights(bucketID i, vector<Factor>& bel) {
  // TODO: figure out if max+, sum+, or sum- for this variable; check weights?
  double Havg = 0.0, wTot = 0.0;
	size_t nNodes = _nodes[i].size();
  vector<double> H(nNodes); // conditional entropies of cliques
  // for each node:
  //   compute conditional entropy
  // TODO::: FINISH !!!!
}


// TODO: what about objective? 
void wmbe::msgForward(bucketID i, double dampTheta, double stepWeights) {
  typedef vector<findex> matchList;
  typedef vector<findex> nodeList;
  double obj = 0.0;
  // Iterate through the variables to eliminate, in order:
//  for (size_t i=0; i<_order.size(); ++i) {
    Var X = Var(_order[i] , 0); // TODO?
    size_t nNodes = _nodes[i].size();       // !!! "  "
    vector<Factor> bel(nNodes);          // Storage for beliefs in this bucket
    if (nNodes > 1) {                    // If more than one mini-bucket:
      if (dampTheta > 0.0 || stepWeights > 0.0)  // compute beliefs if needed for matching
        for (size_t j=0; j<nNodes; ++j) bel[j] = computeNodeBelief( _nodes[i][j] );

      if (dampTheta > 0.0) 
        for (size_t m=0; m<_match[i].size(); ++m) // count over list of matches
          updateTheta(i, dampTheta, _match[i][m], bel);

      if (stepWeights > 0.0)
        updateWeights(i, bel);

      } // end: matching updates if more than one mini-bucket

      // Compute forward messages: 
      for (size_t j=0;j<nNodes;++j) {
        size_t n = _nodes[i][j];
        bel[j] = _theta[n];
        for (nodeList::iterator c=_children[n].begin();c!=_children[n].end();++c) bel[j]+=_msgFwd[*c];
        _msgFwd[j] = (bel[j]*(1.0/_weight[n])).logsumexp(X)*_weight[n];  // take power-lse 
        //_msgFwd[j] = bel[j].logsumexpPower(X, 1.0/_weight[n]);  // take power-lse TODO
        bel[j] = Factor(); // clear memory when done
        if (_parent[n] == -1) obj += _msgFwd[j][0];  // add roots to overall bound
      } // end: forward messages

//    } // end: loop over elimination order 
}

void wmbe::msgBackward(bucketID i, double dampTheta, double stepWeights) {
  typedef vector<findex> nodeList;
	size_t nNodes = _nodes[i].size();
  vector<Factor> bel(nNodes);          // Storage for beliefs in this bucket
  // TODO: consider matching steps here?
 
  // Compute backward message from n -> c
  for (size_t j=0;j<nNodes;++j) {
    size_t n = _nodes[i][j];
    bel[j] = computeNodeBelief( _nodes[i][j] );
    bel[j] -= bel[j].max();
    bel[j] *= (1.0 / _weight[n]);
    for (nodeList::const_iterator c=_children[n].begin();c!=_children[n].end();++c) {
      VarSet vElim = _clique[n] - _clique[*c]; 
      _msgBwd[*c] = bel[j].logsumexp(vElim) * _weight[*c] - _msgFwd[*c];
    } // end: loop over children of n
  }   // end: loop over nodes n in bucket i
}



void wmbe::reparameterize() {
  // replace factors in gm object with thetas?  include msgFwds?
  // update _obj ? discard current mbe?  ???
}




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

  // Simulated version for memory checking (scope only!)
  double scoreSim(const vector<VarSet>& fin, const Var& VX, size_t i, size_t j) {
    double err;
    const VarSet& F1=fin[i], &F2=fin[j];                          // (useful shorthand)
    size_t iBound = std::max(std::max(_iBound,F1.nvar()-1),F2.nvar()-1);      // always OK to keep same size
    size_t sBound = std::max(std::max(_sBound,F1.nrStates()),F2.nrStates());
    VarSet both = F1+F2;
    if ( both.nvar()>iBound+1 || both.nrStates()>sBound ) err=-3;  // too large => -3
    //else if (_byScope) err = 1.0/std::log(both.nrStates()+1);      // greedy scope-based 1 (check if useful???)
    else if (_byScope) err = 1.0/(F1.nvar()+F2.nvar());              // greedy scope-based 2 (check if useful???)
    //else if (_byScope) err = 1;                                    // scope-based => constant score
    else {
      throw std::runtime_error("Cannot simulate dynamic scoring procedure");
    }
    return err;
  }


  void init() {
    _logZ = 0.0;
    if (_order.size()==0) {               // if we need to construct an elimination ordering
      double tic=timeSystem();
      _order=_gmo.order(ordMethod);
      _parents.clear();                   // (new elim order => need new pseudotree) !!! should do together
      std::cout<<"Order in "<<timeSystem()-tic<<" sec\n";
    }
    if (_parents.size()==0) {             // if we need to construct a pseudo-tree
      double tic=timeSystem();
      _parents=_gmo.pseudoTree(_order);
      std::cout<<"Pseudo in "<<timeSystem()-tic<<" sec\n";
    }

    if (_doMplp) {
      mplp _mplp(_gmo.factors());
      char niter[16]; sprintf(niter,"StopIter=%d",_doMplp);
      _mplp.setProperties("Schedule=Fixed,Update=Var,StopObj=-1.0,StopMsg=-1.0");
      _mplp.setProperties(niter);
      _mplp.init();
      _mplp.run();
      //std::cout<<"After Mplp: "<<_mplp.ub()<<"\n";
      _gmo = graphModel(_mplp.beliefs());
    }

    if (_doJG && _doHeur) throw std::runtime_error("MBE: Incompatible options doJG,doHeur");
    std::cout<<"Use scope only? "<<_byScope<<"\n";

    vector<Factor> fin(_gmo.factors());
    for (size_t i=0;i<_gmo.nFactors();++i) fin[i].log();
    vector<double> Norm(_gmo.nFactors(),0.0);
    vector<flist>  vin;

    for (size_t i=0;i<_gmo.nvar();++i) vin.push_back(_gmo.withVariable(var(i)));
    atElim.clear(); atElim.resize(_gmo.nvar());
    atElimNorm.clear(); atElimNorm.resize(_gmo.nvar(),0.0);
    vector<Factor> tmp; if (!_byScope) tmp.resize(_gmo.nFactors());  // temporary factors used in score heuristics
    vector<flist> Orig(_gmo.nFactors());              // origination info: which original factors are
    for (size_t i=0;i<Orig.size();++i) Orig[i]|=i;    //   included for the first time, and which newly
    vector<flist> New(_gmo.nFactors());                //   created clusters feed into this cluster


    //// Eliminate each variable in the sequence given: /////////////////////////////////
    for (VarOrder::const_iterator x=_order.begin();x!=_order.end();++x) {
      //std::cout<<"Eliminating "<<*x<<"\n";
      Var VX = var(*x);
      if (*x >= vin.size() || vin[*x].size()==0) continue;    // check that we have some factors over this variable

      flist ids = vin[*x];                 // list of factor IDs contained in this bucket

      //// Select allocation into buckets ///////////////////////////////////////
      typedef flist::const_iterator flistIt;
      typedef std::pair<double,sPair> _INS;
      std::multimap<double,sPair > scores;
      std::map<sPair,std::multimap<double,sPair>::iterator> reverseScore;

      //std::cout<<"Initial table sizes: "; for (flistIt i=ids.begin();i!=ids.end();++i) std::cout<<fin[*i].numel()<<" "; std::cout<<"\n";

      //// Populate list of pairwise scores for aggregation //////////////////////
      if (!_byScope) for (flistIt i=ids.begin();i!=ids.end();++i) tmp[*i] = elim(fin[*i],VX);
      for (flistIt i=ids.begin();i!=ids.end();++i) {
        for (flistIt j=ids.begin(); j!=i; ++j) {
          double err = score(fin,VX,*i,*j,tmp); sPair sp(*i,*j);
          reverseScore[sp]=scores.insert(_INS(err,sp));       // save score 
        }
        reverseScore[sPair(*i,*i)]=scores.insert(_INS(-1,sPair(*i,*i)));         // mark self index at -1
      }

      //// Run through until no more pairs can be aggregated: ////////////////////
      //   Find the best pair (ii,jj) according to the scoring heuristic and join
      //   them as jj; then remove ii and re-score all pairs with jj
      for(;;) {
        std::multimap<double,sPair>::reverse_iterator top = scores.rbegin();
         //multimap<double,_IDX>::reverse_iterator  last=scores.lower_bound(top->first);  // break ties randomly !!!
        //std::advance(last, randi(std::distance(top,last)));
        //std::cout<<top->first<<" "<<top->second.first<<" "<<top->second.second<<"\n";

        if (top->first < 0) break;                            // if can't do any more, quit
        else {
          size_t ii=top->second.first, jj=top->second.second;
          //std::cout<<"Joining "<<ii<<","<<jj<<"; size "<<(fin[ii].vars()+fin[jj].vars()).nrStates()<<"\n";
//!!!          if (fin[ii].vars()>>fin[jj].vars()) { fin[jj].swap(fin[ii]); }
          fin[jj] += fin[ii];                             // combine into j
          Norm[jj] += Norm[ii];
          ////double mx = fin[jj].max(); fin[jj]/=mx; mx=std::log(mx); _logZ+=mx; Norm[jj]+=mx;
          erase(vin,ii,fin[ii].vars()); fin[ii]=Factor(0.0);  //   & remove i

          Orig[jj] |= Orig[ii]; Orig[ii].clear();      // keep track of list of original factors in this cluster
          New[jj]  |= New[ii];  New[ii].clear();       //   list of new "message" clusters incoming to this cluster

          if (!_byScope) tmp[jj] = elim(fin[jj],VX);       // update partially eliminated entry for j

          for (flistIt k=ids.begin();k!=ids.end();++k) {     // removing entry i => remove (i,k) for all k
            scores.erase(reverseScore[sPair(ii,*k)]);
          }
          ids /= ii;

          for (flistIt k=ids.begin();k!=ids.end();++k) {  // updated j; rescore all pairs (j,k) 
            if (*k==jj) continue;
            double err = score(fin,VX,jj,*k,tmp); sPair sp(jj,*k);
            scores.erase(reverseScore[sp]);                        // change score (i,j)
            reverseScore[sp]=scores.insert(_INS(err,sp));  //
          }
        }
      }

      //std::cout<<"End table sizes: "; for (flistIt i=ids.begin();i!=ids.end();++i) std::cout<<fin[*i].numel()<<" "; std::cout<<"\n";

      //// Perform any matching? /////////////////////////////////////////////////
      //    "Matching" here is: compute the largest overlap of all buckets, and ensure that the
      //    moments on that subset of variables are identical in all buckets.
      //    !!! also, add extra variables if we can afford them?
      //size_t beta=0;
      if (_doMatch && ids.size()>1) {

        // Possibly we should first extend the scopes of the buckets to be larger
        if (_doFill) {
          VarSet all=fin[ids[0]].vars(); for (size_t i=1;i<ids.size();i++) all|=fin[ids[i]].vars();
          //reorder all?
          for (size_t i=0;i<ids.size();i++) {
            VarSet vsi = fin[ids[i]].vars(), orig=fin[ids[i]].vars();
            for (VarSet::const_iterator vj=all.begin();vj!=all.end();++vj) {
              VarSet test = vsi + *vj;
              if (test.nvar() <= _iBound+1 && test.nrStates() <= _sBound) vsi=test;
            }
            if (vsi.nvar()>orig.nvar()) fin[ids[i]] += Factor(vsi-orig,0.0);
          }
        }

        vector<Factor> ftmp(ids.size());            // compute geometric mean
        VarSet var = fin[ids[0]].vars();            // on all mutual variables
        for (size_t i=1;i<ids.size();i++) var &= fin[ids[i]].vars();
        Factor fmatch(var,0.0);
        for (size_t i=0;i<ids.size();i++) {
          ftmp[i] = marg(fin[ids[i]],var);
          fmatch += ftmp[i];
        }
        fmatch *= (1.0/ids.size());                  // and match each bucket to it
        for (size_t i=0;i<ids.size();i++) fin[ids[i]] += fmatch - ftmp[i];

        //beta = addFactor( Factor(fmatch.vars(),1.0) );  // add node to new cluster graph
        //atElim[*x] |= beta;
      }


      //// Weight heuristic? /////////
      //  Currently, we just take the first bucket; !!! add heuristics from matlab code
      size_t select = ids[0];

      //// Eliminate individually within buckets /////////////////////////////////
      //   currently does not use weights except 0/1; !!! add sumPower alternatives from matlab code
      vector<findex> alphas;
      //std::cout<<"Table sizes: ";
      for (flistIt i=ids.begin();i!=ids.end();++i) {
        // 
        // Create new cluster alpha over this set of variables; save function parameters also
        findex alpha = findex(-1), alpha2 = findex(-1);
        if (_doJG) { alpha=addFactor(fin[*i]); alphas.push_back(alpha); }

        //std::cout<<fin[*i].numel()<<" ";

        //if (*i==select) fin[*i] = elim     (fin[*i],VX);
        //else            fin[*i] = elimBound(fin[*i],VX);  
        if (*i==select) { Factor FE = elim     (fin[*i],VX); fin[*i].swap(FE); }
        else            { Factor FE = elimBound(fin[*i],VX); fin[*i].swap(FE); }

        if (_doJG) _factors[alpha] -= fin[*i];

        if (_doHeur) alpha2 = addFactor(fin[*i]);          /// !!! change to adding a message not a factor?

        //if (_doJG && _doMatch && ids.size()>1) addEdge(alpha,beta);
        if (_doJG) for (size_t j=0;j<alphas.size()-1;++j) addEdge(alpha,alphas[j]);
        if (_doJG) for (flistIt j=New[*i].begin();j!=New[*i].end();++j) addEdge(*j,alpha);

        Orig[*i].clear(); New[*i].clear(); New[*i]|=alpha;  // now incoming nodes to *i is just alpha

        size_t k=_parents[*x];                        //  mark next bucket and intermediates with msg for heuristic calc
        for (; k!=vindex(-1) && !fin[*i].vars().contains(var(k)); k=_parents[k]) { atElim[k]|=alpha2; atElimNorm[k]+=Norm[*i]; }
        if (k!=vindex(-1)) { atElim[k] |= alpha2; atElimNorm[k]+=Norm[*i]; }  // need check?

        insert(vin,*i,fin[*i].vars());              // recompute and update adjacency

        // Shortcut return if model is inconsistent (value = -inf)
        if (VX.states()>0 && fin[*i].max() == -infty()) { _logZ=_bound=-infty(); _dObj=_dObjCurr=0.0; _iter=0; return; }
      }
      //std::cout<<"\n";

    }
    /// end for: variable elim order /////////////////////////////////////////////////////////

    Factor F(0.0);
    for (size_t i=0;i<fin.size();++i) F+=fin[i];
    assert( F.nvar() == 0 );
    _logZ += F[0];

    _bound=_logZ; _dObj=_dObjCurr=infty(); _iter=0;
    //std::cout<<"Bound "<<_logZ<<"\n";
  }





*/


}

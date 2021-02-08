// Revised by qlou on 09//11/15: set isAndOr = true for pseudoTree in wmbe::init()
// Revised by qlou on 09/15/15: make getOrder, getPseudotree be const function
// Revised by qlou on 10/31/15: add a function "duplicateInBucket" to find out whether the node lives in more than one mini-bucket
// Revised by qlou on 1/4/16: add a function "eliminateFromNow" to allow exact elimination from a given status.
// Revised by qlou on 1/7/16: add a function "getCurrentModel" to get current model (factors) to enable exact elimination from that

#ifndef __MEX_MINIBUCKET_H
#define __MEX_MINIBUCKET_H

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <set>

#include <list>

#include "graphmodel.h"
#include "alg.h"


namespace mex {

// Graphical Model Algorithm: Weighted Mini-bucket Elimination 
// 


  namespace detail {
		class jgNode {
			public:
			typedef jgNode* nodeID;

			// two parameters added by Qi on 1/6/18 for the augmented model
			// they help identify whether the forward message is sent from a SUM node to a MAX node
			// vid: variable with which this node (mini-bucket) associate.
			// isMsgFwdRep: whether its forward message have to be raised to the K-th power to compute the total bound
			// 	true for two cases: 1) a root SUM variable 2) a SUM variable whose parent is a MAX variable
			graphModel::vindex vid = -1; // default -1
			bool isMsgFwdRep = false; // default: false

			//
			VarSet clique;
			Factor theta;
			nodeID parent;
			vector<nodeID> children;
			Factor msgFwd;
			Factor msgBwd;
			double weight;
			vector<graphModel::findex> origFactors;
			Factor belief;
			jgNode() { clear(); };
			void clear() {
				clique.clear();
				theta = msgFwd = msgBwd = belief = Factor(0.0);
  			parent = NULL; children.clear();
				weight = 0.0;
			}
		};
		struct nodePair : public std::pair<jgNode::nodeID,jgNode::nodeID> {
    	nodePair(jgNode::nodeID ni, jgNode::nodeID nj) { if (ni<nj) {first=nj; second=ni;} else {first=ni;second=nj;} }
  	};
	}

// TODO: convert nodeID from jgNode* to index into vector<jgNode>

class wmbe : public graphModel, public gmAlg, virtual public mxObject {

// TODO: resolve "containing" gmo, containing theta, and  "being" a gmodel

public:
  typedef graphModel::findex    	findex;        // factor index
  typedef graphModel::vindex    	vindex;        // variable index
  typedef graphModel::flist     	flist;      // collection of factor indices
	typedef detail::jgNode 			 		jgNode;
	typedef detail::jgNode::nodeID	nodeID;
	typedef size_t            			bucketID;
	typedef vector<nodeID> 					matchList;
	typedef vector<nodeID> 					nodeList;

public:
  wmbe() : graphModel() { setProperties(); }
  wmbe(vector<Factor> fs) : graphModel(fs) { setProperties(); }
  wmbe(const graphModel& gm) : graphModel(gm) { setProperties(); }
  //wmbe(const wmbe& gm) : { setProperties(); *this = gm; }
  virtual wmbe* clone() const { wmbe* gm = new wmbe(*this); return gm; }

  //wmbe& operator= (wmbe const& gm);          // assignment (deep copy) 

  void            setOrder(const VarOrder& order) { _order = order; }
  // TODO: how to express constrained elimination orders???
  //void            setOrder(const VarOrder& order, const vector<double>& weights=vector<double>());
  //void            setOrder(graphModel::OrderMethod method, const vector<double>& weights=vector<double>());
  //
  //void            setWeight(size_t vIndex, double weight);
  const VarOrder& getOrder()  const              { return _order; }

  // Qi: 6/27/16
  const vector<size_t>& getPriority() const { return _priority; };

  void   setIBound(size_t i) { _iBound=i ? i : std::numeric_limits<size_t>::max(); }
  size_t getIBound() const   { return _iBound; }
  void   setSBound(size_t s) { _sBound=s ? s : std::numeric_limits<size_t>::max(); }
  size_t getSBound() const   { return _sBound; }

  const vector<vindex>& getPseudotree()   const                { return _parents; }
  void                  setPseudotree(const vector<vindex>& p) { _parents = p;    }

  /////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // set max variables and set their weights as well
  void setMaxVars(const VarSet& maxVars ) {
	  _maxVars = maxVars;
  }
  VarSet getMaxVars() { return _maxVars;}
  size_t nMaxVar() const { return _maxVars.size(); }

  // Build mini-bucket data structure
  void init(const VarSet& vs) { init(); }            // (meaning here?)
  void init();
  void init(bool isAndOr); // Qi: almost identical to the previous init() except allowing using AND/OR structure. 6/16/2016
  virtual void run() { }  // TODO: add

  void build();
  void setMsgFwdRep(const int K); // Qi: 1/6/18
  void setTheta();
  void dump(std::ostream&);
  nodeID addNodeBasic( VarSet vs );
  double score( Var, nodeID, nodeID );

  template <class MapType>
  double maxSequential(MapType& config);

  template <class MapType>
  std::pair<double,double> sampleSequential(MapType& config);

  template <class MapType>
  std::pair<double,double> sampleMixture(MapType& config);

  /////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
	// add by Qi Lou on 6/7/2016
	template <class MapType>
	std::vector<double> sampleMixtureMore(MapType& config, std::vector<double>& condBD, std::vector<double>& thetaVec);


	// add by Qi Lou on 6/26/2016
	template <class MapType>
	double sampleMixtureConditioned(int ID, MapType& config, std::vector<double>& logQxVec);

	// add by Qi Lou on 10/18/2016
	// AND/OR version of sampleMixtureConditioned
	template <class MapType>
	void conditionalMixtureSampling(int ID, MapType& config, std::vector<double>& logQxVec, const std::list<Var>& Undone);
	// QLOU on 11/15/2017
	// customized version for "randheur"
	template <class MapType>
	double conditionalMixtureSampling(const int ID, MapType& config, const std::list<vindex>& rts,
			const std::vector<std::list<vindex> >& ascList, const std::vector<int>& varDepthVec, const int maxDepth);
	// Qi: 1/13/18, customized version for mmapIS
	template <class MapType>
	void conditionalMixtureSampling(const int ID, MapType& config, const int K,
			const std::list<Var>& unDoneMax, const std::list<Var>& unDoneSum, double& logFx, double& logQx, double&, double&);
	// Qi: 2/12/18, fast version for mmapIS
	template <class MapType>
	void conditionalMixtureSampling(const int ID, MapType& config, const int K, const int start_loc, const int end_loc,
			const std::vector<mex::Var>& pstree, const std::vector<bool>& varTypeVec, const bool MAX_VAR,
			double& logFx, double& logQx, double& aveEst, double& prodEst);


  /////////////////////////////////////////////////////////////////
  // Setting properties (directly or through property string)
  /////////////////////////////////////////////////////////////////
  MEX_ENUM( Property , iBound,sBound,Order,Distance,DoMatch,DoMplp,DoFill,DoJG,DoHeur );

  virtual void setProperties(std::string opt=std::string()) {
    if (opt.length()==0) {
      setProperties("iBound=4,sBound=inf,Order=MinWidth,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=1");
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
        //case Property::Distance: distMethod = Factor::Distance(asgn[1].c_str()); _byScope=false; break;
        //case Property::DoMplp:  _doMplp  = atol(asgn[1].c_str()); break;
        //case Property::DoMatch: _doMatch = atol(asgn[1].c_str()); break;
        //case Property::DoFill:  _doFill = atol(asgn[1].c_str()); break;
        //case Property::DoJG:    _doJG = atol(asgn[1].c_str()); break;
        //case Property::DoHeur:  _doHeur = atol(asgn[1].c_str()); break;
        default: break;
      }
    }
  }



  // Weight init
  MEX_ENUM( ElimType , MaxUpper,SumUpper,SumLower,SumBethe );

  void setElimType(Var v, ElimType e);
  void setElimType(ElimType e) { for (size_t v=0;v<nvar();++v) setElimType(var(v),e); }

  // TODO: deprecated methods
  void setWeightPositiveFirst(bucketID b, double tot);
  void setWeightNegativeFirst(bucketID b, double tot);
  void setWeightPositiveUniform(bucketID b, double tot);
  void setWeightNegativeUniform(bucketID b, double tot, double pos=-1.0);


  // Message passing updates 
  // TODO: "run" method with simple controls
  void msgForward(double dampTheta, double stepWeights);
  void msgBackward();

  // Fine-grain algorithmic control
  Factor computeNodeBelief(nodeID n) const;
  void   updateTheta(bucketID b, double damp, const matchList& match, bool updateBeliefs);
  void   updateWeights(bucketID b, double step);
  double msgForward(bucketID b, double dampTheta, double stepWeights);
  double msgForwardAug(bucketID b, double dampTheta, double stepWeights, const int K=1); // Qi: 12/2/2017
  void   msgBackward(bucketID b, double dampTheta, double stepWeights);

  void   reparameterize();

  //template <class MapType>
  //double logHeurToGo(Var v, MapType vals) const;  // evaluate remaining cost given context

  template <class MapType>
  double heuristicTheta(Var v, MapType vals) const;   // evaluate heuristic costs for a bucket
  template <class MapType>
  double heuristicIn(Var v, MapType vals) const;  //  "" ""

  /////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Qi: 1/7/18: the following is for the augmented model case
  template <class MapType>
  double heuristicInAug(Var v, MapType vals, int K, bool isMaxVar) const;
  // debug
  template <class MapType>
  unsigned duplicateInBucket(Var v, MapType vals) const;  // nodes in more than one mini-buckets
  unsigned getNumberOfMiniBuckets(Var v) const;  // get the number of mini-buckets for a variable
  double eliminateFromNow(std::set<Var>& frontier, std::set<Var>& subtree) const;
  vector<Factor> getCurrentModel(std::set<Var>& frontier, std::set<Var>& subtree) const;

  double memory();      // return amt of memory (mb) used by the algorithm  :  TODO: iterative? sparse? ???


private:
  // Contains:
  VarOrder       _order;           // v = order[i] is the ith variable eliminated
  vector<size_t> _priority;        // i = priority[v] is the step in which variable v is eliminated
  vector<double> _varWeight;       // TODO: need to keep?
  vector<vector<nodeID> > _nodes;  // nodes[i][k] is the kth node in bucket i
  std::list<jgNode>   _nodeList;   // nodeList[n] is the node associated with ID n
	vector<vector<vector<nodeID> > > _match;  // match[b] = list of sets of nodes to match on in bucket b
	vector<vindex> _parents;  			 // pseudotree data structure  (TODO use?)
  vector<nodeID> _factorIn;
  vector<ElimType> _elim;          // elimination method / weight
  vector<vector<nodeID> > _heuristic;  // nodes whose messages are part of the heuristic function

  const jgNode& node(nodeID n) const { return *n; }
        jgNode& node(nodeID n)       { return *n; }

  size_t _iBound;
  double _sBound;
  double _obj; 
  double _dampTheta;   // damping factor (1.0 = no damping; 0.0 = no updates)
  double _dampWeight;  // step size for weight updates

  Factor::Distance distMethod;
  graphModel::OrderMethod ordMethod;

  /////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // max variables
  VarSet _maxVars;
  // no. of replicates of sum variables, for mmapIS only
  int _K = 1; // initialized to 1, will be assigned in "setMsgFwdRep"

public:

#ifdef MEX  
  // MEX Class Wrapper Functions //////////////////////////////////////////////////////////
  //bool        mxCheckValid(const mxArray*);   // check if matlab object is compatible with this object
  //void        mxSet(mxArray*);     // associate with A by reference to data
  //mxArray*    mxGet();             // get a pointer to the matlab object wrapper (creating if required)
  //void        mxRelease();         // disassociate with a matlab object wrapper, if we have one
  //void        mxDestroy();         // disassociate and delete matlab object
  //void        mxSwap(mbe& gm);     // exchange with another object and matlab identity
  /////////////////////////////////////////////////////////////////////////////////////////
  bool mxCheckValid(const mxArray* GM) { throw std::runtime_error("NOT IMPLEMENTED"); }
  void mxSet(mxArray* GM) { throw std::runtime_error("NOT IMPLEMENTED"); }
  mxArray* mxGet() { if (!mxAvail()) { graphModel::mxGet(); } return M_; }
  void mxRelease() { throw std::runtime_error("NOT IMPLEMENTED"); }
  void mxDestroy() { throw std::runtime_error("NOT IMPLEMENTED"); }
  void mxSwap(wmbe& gm) { graphModel::mxSwap( (graphModel&)gm ); }
#endif

  // Can be an optimization algorithm or a summation algorithm....
  double ub() const { throw std::runtime_error("NOT IMPLEMENTED"); }
  double lb() const { throw std::runtime_error("NOT IMPLEMENTED"); }
  vector<index> best() const { throw std::runtime_error("Not implemented"); }

  double logZ()   const { throw std::runtime_error("NOT IMPLEMENTED"); }
  double logZub() const { throw std::runtime_error("NOT IMPLEMENTED"); }
  double logZlb() const { throw std::runtime_error("NOT IMPLEMENTED"); }

  // No beliefs defined currently
  //Factor belRef;
  const Factor& belief(size_t f)  const { throw std::runtime_error("Not implemented"); }
  const Factor& belief(Var v)     const { throw std::runtime_error("Not implemented"); }
  //const Factor& belief(size_t f)  const { return belRef = computeNodeBelief(_factorIn[f]).marginal(_factors[f].vars()); }
  //const Factor& belief(Var v)     const { return belRef = computeNodeBelief( *_nodes[_priority[v]].begin() ).marginal(v); }
  const Factor& belief(VarSet vs) const { throw std::runtime_error("Not implemented"); }
  const vector<Factor>& beliefs() const { throw std::runtime_error("Not implemented"); }


  size_t _iter;
  double _dObj,_dObjCurr;

  double iter()          { return ((double)_iter)/edges().size(); }
  bool   iter_boundary() { return (_iter % edges().size()) == 0; }
  double dObj()          { return std::abs(_dObj); }



};


/////////////////// IMPLEMENTATION ////////////////////////////////////////

/*
template <class MapType>
double wmbe::heuristicPre(Var v, MapType config) const {
  size_t b = _priority[v];  // find bucket to examine
  double heur = 0.0;
  for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
    heur += node(*ni).msgFwd[ sub2ind(node(*ni).msgFwd.vars(),config) ]; // compute contributions of forward messages
  }
  return heur;
}
template <class MapType>
double wmbe::heuristicPost(Var v, MapType config) const {
  size_t b = _priority[v];  // find bucket to examine
  double heur = 0.0;
  for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
    heur += node(*ni).theta[ sub2ind(node(*ni).theta.vars(),config) ];
    //heur -= node(*ni).msgFwd[ sub2ind(node(*ni).msgFwd.vars(),config) ]; // subtract forward message as already included?
    for (nodeList::const_iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci)
      heur += node(*ci).msgFwd[ sub2ind(node(*ci).msgFwd.vars(),config) ];
  }
  return heur;
}

*/

template <class MapType>
double wmbe::heuristicTheta(Var v, MapType config) const {
  size_t b = _priority[v];  // find bucket to examine
  double heur = 0.0;
  for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
    heur += node(*ni).theta[ sub2ind(node(*ni).theta.vars(),config) ]; // compute contributions of theta
  }
  return heur;
}

template <class MapType>
double wmbe::heuristicIn(Var v, MapType config) const {
  size_t b = _priority[v];  // find bucket to examine
  double heur = 0.0;
  for (nodeList::const_iterator ci=_heuristic[b].begin();ci!=_heuristic[b].end();++ci) 
    heur += node(*ci).msgFwd[ sub2ind(node(*ci).msgFwd.vars(),config) ];
  //for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) 
  //  for (nodeList::const_iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci)
  //    heur += node(*ci).msgFwd[ sub2ind(node(*ci).msgFwd.vars(),config) ];
  return heur;
}

/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
// Qi: 1/7/18: the following is for the augmented model case
template <class MapType>
double wmbe::heuristicInAug(Var v, MapType config, int K, bool isMaxVar) const {
	// TODO "MapType& config"? to save time?
  size_t b = _priority[v];	// find bucket to examine
	double heur = 0.0;
	for (nodeList::const_iterator ci=_heuristic[b].begin();ci!=_heuristic[b].end();++ci) {
		// MAX variable, only when the message is from a SUM descendant
		if (isMaxVar && node(*ci).isMsgFwdRep) {
			heur += K * node(*ci).msgFwd[ sub2ind(node(*ci).msgFwd.vars(),config) ];
		} else {
			heur += node(*ci).msgFwd[ sub2ind(node(*ci).msgFwd.vars(),config) ];
		}
		// SUM variable, no change, although it's possible  node(*ci).isMsgFwdRep = true
	}
	//for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni)
	//	for (nodeList::const_iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci)
	//		heur += node(*ci).msgFwd[ sub2ind(node(*ci).msgFwd.vars(),config) ];
	return heur;
}

// debug, find out whether the node lives in more than one mini-bucket
// search results should be improved if so
// return no. of mini-buckets containing that variable
template <class MapType>
unsigned wmbe::duplicateInBucket(Var v, MapType config) const {
//	unsigned nMiniBuckets = 0;
//	size_t b = _priority[v];	// find bucket to examine
//	// each node here seems to be a mini-bucket
//	// actually, _nodes[b].size() would serve the purpose, see wmbe.cpp
//	for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
//		if ((node(*ni).clique).contains(v)) ++nMiniBuckets;
//	}
//	unsigned nNodes = _nodes[b].size();
//
//	// debug
//	assert( nMiniBuckets == nNodes );
//
//	return nMiniBuckets;
	return _nodes[_priority[v]].size();
}
template <class MapType>
double wmbe::maxSequential(MapType& config) { 
  VarSet done;
  double logP=0.0;
  for (int b=_nodes.size()-1;b>=0;--b) {
    Factor bel(0.0);
    for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
      VarSet condv = node(*ni).theta.vars() & done;
      bel += node(*ni).theta.condition( condv , sub2ind(condv,config) );
      if (node(*ni).parent != NULL)
        bel -= node(*ni).msgFwd[ sub2ind(node(*ni).msgFwd.vars(),config) ];
      for (nodeList::iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
        bel += node(*ci).msgFwd.condition( condv, sub2ind(condv,config) );
      }
    }
    assert( bel.nvar() == 1 );
    config[ bel.vars()[0] ] = bel.argmax(); 
    logP += bel[ config[ bel.vars()[0] ] ];
    done |= bel.vars();
  }
  return logP;
}


template <class MapType>
std::pair<double,double> wmbe::sampleSequential(MapType& config) { 
  VarSet done;
  double logP=0.0, logQ=0.0;
  for (int b=_nodes.size()-1;b>=0;--b) {
    Factor bel(0.0);
    for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
      VarSet condv = node(*ni).theta.vars() & done;
      bel += node(*ni).theta.condition( condv , sub2ind(condv,config) );
      if (node(*ni).parent != NULL)
        bel -= node(*ni).msgFwd[ sub2ind(node(*ni).msgFwd.vars(),config) ];
      for (nodeList::iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
        bel += node(*ci).msgFwd.condition( condv, sub2ind(condv,config) );
      }
    }
    assert( bel.nvar() == 1 );
    Factor proposal = bel; proposal -= proposal.max(); proposal.exp(); proposal /= proposal.sum();
    config[ bel.vars()[0] ] = proposal.sample();   
    logQ += std::log(proposal[ config[ bel.vars()[0] ] ]);
    logP += bel[ config[ bel.vars()[0] ] ];
    done |= bel.vars();
  }
  return std::pair<double,double>(logP,logQ);
}

template <class MapType>
std::pair<double,double> wmbe::sampleMixture(MapType& config) {
  VarSet done;
  double logPx=0.0, logQx=0.0;
  for (int b=_nodes.size()-1;b>=0;--b) {
    Factor qMix(0.0);
    for (nodeList::const_iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
      Factor qComponent(0.0);
      VarSet condv = node(*ni).theta.vars() & done;
      qComponent += node(*ni).theta.condition( condv , sub2ind(condv,config) );
      qComponent -= node(*ni).msgFwd[ sub2ind(node(*ni).msgFwd.vars(),config) ];
      for (nodeList::const_iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
        qComponent += node(*ci).msgFwd.condition( condv, sub2ind(condv,config) );
      }
      //qComponent -= qComponent.max();  // TODO: subtract & re-add max for numerical stability?
      qComponent *= (1.0 / node(*ni).weight);
      qComponent.exp();
      if (qComponent.sum() == 0.0) {  // Sampling failure! p(x)=0
        // return std::pair<double,double>(-infty(),-1e100);   // Simpler: return p(x) and arbitrary q(x)
        for (int bb=b;bb>=0;--bb) {   // More involved: finish sampling x uniformly and return p(x) and q(x)
          config[_order[bb]] = size_t( var(_order[bb]).states() * mex::randu() );
          logQx -= std::log( var(_order[bb]).states() );
        }
        return std::pair<double,double>(-infty(),logQx); 
      }
      //if (qComponent.sum() == 0.0) qComponent += 1.0/qComponent.nrStates(); // Alternate form of uniform above
      qMix += node(*ni).weight * qComponent;
    }
    qMix /= qMix.sum();
    size_t x_value = qMix.sample();
    config[ qMix.vars()[0] ] = x_value; 
    logQx += std::log( qMix[ x_value ] ); 
    done |= qMix.vars();
  }
  logPx = logP( config );
  return std::pair<double,double>(logPx,logQx);

/* TODO: slightly faster?  mostly lookup (non-selected components), vs unary factor constructions for other impl.
  VarSet done;
  double logPx=0.0, logQx=0.0;
  for (int b=_nodes.size()-1;b>=0;--b) {
		double wTot = 0.0; for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) wTot += node(*ni).weight;
		double wDraw = wTot * mex::randu();
		nodeList::iterator nDraw=_nodes[b].begin();
		for (;nDraw!=_nodes[b].end();++nDraw) if ( (wDraw -= node(*nDraw).weight) < 0 ) break;
		// Draw sample from selected node belief *nDraw
		{ Factor bel(0.0);	// TODO: *ni not set
      VarSet condv = node(*nDraw).theta.vars() & done;
      bel += node(*nDraw).theta.condition( condv , sub2ind(condv,config) );
			bel -= node(*nDraw).msgFwd[ sub2ind(node(*nDraw).msgFwd.vars(),config) ];	
      for (nodeList::iterator ci=node(*nDraw).children.begin();ci!=node(*nDraw).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
        bel += node(*ci).msgFwd.condition( condv, sub2ind(condv,config) );
    	}
			bel -= bel.max();
			bel *= (1.0 / node(*nDraw).weight);
			bel.exp();
			if (bel.sum() == 0.0) bel += 1.0/bel.nrStates();
			config[ bel.vars()[0] ] = bel.sample();			
    	done |= bel.vars();
		}
		// Compute the overall mixture probability q(x)
		double qMix = 0.0;
    for (nodeList::iterator ni=_nodes[b].begin();ni!=_nodes[b].end();++ni) {
    	Factor bel(0.0);
      VarSet condv = node(*ni).theta.vars() & done;
      bel += node(*ni).theta.condition( condv , sub2ind(condv,config) );
			bel -= node(*ni).msgFwd[ sub2ind(node(*ni).msgFwd.vars(),config) ];	
			if (node(*ni).parent == NULL) logPx += node(*ni).msgFwd[0]; // sub2ind(node(*ni).msgFwd.vars(),config) ];
      for (nodeList::iterator ci=node(*ni).children.begin();ci!=node(*ni).children.end();++ci) {
        VarSet condv = node(*ci).msgFwd.vars() & done;
        bel += node(*ci).msgFwd.condition( condv, sub2ind(condv,config) );
      }
			logPx += bel[0]; // config[ bel.vars()[0] ] ];   // contribution of this node to log P(x)
			qMix += std::exp( bel[0] * (1.0 / node(*ni).weight) + std::log(node(*ni).weight/wTot) ); 
// config[ bel.vars()[0] ] ] * (1.0 / node(*ni).weight);	// mixture component probability
    }
    logQx += std::log(qMix);
  }
  return std::pair<double,double>(logPx,logQx);
*/
}

/////// TODO QI LOU CODE  /////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * add by Qi Lou on 6/7/2016
 * return conditional samples at each variable
*/
template <class MapType>
std::vector<double> wmbe::sampleMixtureMore(MapType& config, std::vector<double>& condBD, std::vector<double>& thetaVec) {
	/*
	 * A more informative version of the function "sampleMixture"
	 * condBD is a vector (with size nvar+1) of upper bounds on logZ_{X_A | X_B} w.r.t. the current configuration, where X = [X_A, X_B]
	 * when B is the empty set, logZ_{X_A | X_B} = logZ.
	 * condBD[i] is (upper) bound when the i-th variable in the elimination order and the rest variables behind it are instantiated.
	 * condBD[0] is the exact value for that configuration, condBD[nvar] is the wmb (upper) bound on logZ
	 * logQxVec[i] is the value (w.r.t. to current configuration) of the conditional distribution
	 * given the i-th variable in the elimination order and the rest variables behind it, i.e., q(X_{0,...,i-1} | X_{i,...,nvar-1}).
	 * logQxVec[0] is not properly defined. logQxVec[nvar] is log{q(X)}
	 * thetaVec[i] is the sum of costs when x_i... x_{nvar-1} are instantiated by config
	 * thetaVec[0] is the value of that config
	 * thetaVec[nvar] is not properly defined, set to 0.0 for implementation convenience
	 */
	// for convenience, the length of those two vectors are nvar+1
	assert(condBD.size() == nvar()+1);
	assert(_order.size() == nvar() );
	std::fill(condBD.begin(), condBD.end(), 0.0); // for safety, set all to be zeros
	std::vector<double> logQxVec(condBD); // initialized to zeros
	std::fill(thetaVec.begin(), thetaVec.end(), 0.0); // for safety, set all to be zeros

	VarSet done;
	double logPx = 0.0, logQx = 0.0;
	for (int b = _nodes.size() - 1; b >= 0; --b) {
		Factor qMix(0.0);
		for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
			Factor qComponent(0.0);
			VarSet condv = node(*ni).theta.vars() & done;
			qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
			qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
			for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
				VarSet condv = node(*ci).msgFwd.vars() & done;
				qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
			}
//			qComponent -= qComponent.max();  // this causes a bug
			qComponent *= (1.0 / node(*ni).weight);
			qComponent.exp();
			if (qComponent.sum() == 0.0)
				qComponent += 1.0 / qComponent.nrStates();
			qMix += node(*ni).weight * qComponent;
		}
		qMix /= qMix.sum();
		size_t tmp = qMix.sample();
		config[qMix.vars()[0]] = tmp; //qMix.sample();
		double tmp_val = std::log(qMix[tmp]); //config[ qMix.vars()[0] ] ] );
		logQx += tmp_val;
		done |= qMix.vars();
		// q(x_i | par(x_i))
		logQxVec[b+1] = tmp_val;
	}
	// q(x_0,...,x_{i-1} | x_i, ... x_{nvar-1} ) = q(x_0 | par(x_0)) *...* q(x_{i-1} | par(x_{i-1}))
	for (int i=1; i<logQxVec.size(); ++i) {
		logQxVec[i] += logQxVec[i-1]; // logQxVec[0] = 0.0
	}

	// Compute condBD
	for (int i = _order.size() - 1; i >= 0; --i) {
		auto v = _order[i]; // v = _order[i] is the i-th variable eliminated
		// priority[v] is the step in which variable v is eliminated
		size_t b = _priority[v];	// find bucket to examine
		assert(i==int(b)); // debug

		// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
		double heurTheta = heuristicTheta(var(v), config); // pass argument by value, won't change config;
		// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
		double heurIn = heuristicIn(var(v), config);  // pass argument by value, won't change config;

		thetaVec[i] = heurTheta + thetaVec[i+1]; // thetaVec[_order.size()] = 0.0;
		condBD[i] = heurIn + thetaVec[i];
	}

	// compute (upper) bound on logZ
	auto v = _order.back();
	Var X = var(v);
	mex::vector<uint32_t> tuple(nvar());
	Factor db(X,0.0);
	for (size_t i=0; i<X.states(); ++i) {
		tuple[ X ] = i;
		db[i] = heuristicTheta(X,tuple) +  heuristicIn(X,tuple);
	}
	condBD.back() = db.logsumexp(); // wmb bound on logZ

	//
	logPx = logP(config);
//	condBD[0] = logPx;
	// check, these  quantities should be the same
	assert( condBD[0]<logPx+1e-5 && condBD[0]>logPx-1e-5);
	assert( condBD[0]<thetaVec[0]+1e-5 && condBD[0]>thetaVec[0]-1e-5);

	return logQxVec;
}

/*
 * add by Qi Lou on 6/26/2016
 * condition on the last several variables (based on on the order),
 * and sample from the mixture distribution of the rest variables.
*/
template <class MapType>
double wmbe::sampleMixtureConditioned(int ID, MapType& config, std::vector<double>& logQxVec) {
	/*
	 * only support chain-type pseudo tree
	 * ID is the unique variable id for the foremost variable conditioned in the order.
	 * ID < 0 means no conditioned variables.
	 * input config is already instantiated by conditioned variables
	 * the input logQxVec is initialized with conditional distribution of those conditioned variables
	 * i.e., if conditioned on the i-th variable and its ancestors, for any j>=i, logQxVec[j] = q(x_{j-1} | par(x_{j-1}) ) is already set;
	 * for any j<i, logQxVec[j] has not been set;
	 *
	 * logQxVec[i] is the value (w.r.t. to current configuration) of the conditional probability
	 * given the i-th variable in the elimination order and the rest variables behind it, i.e., q(X_{0,...,i-1} | X_{i,...,nvar-1}).
	 * logQxVec[0] is not properly defined. logQxVec[nvar] is log{q(X)}
	 */
//	std::vector<double> logQxVec(nvar()+1, 0.0); // initialized to zeros
	assert( logQxVec.size() == nvar()+1  );
	logQxVec[0] = 0.0;
	// position of the variable in the elimination order
	int loc = nvar(); // default value means no conditioning, ID = -1
	if ( ID >= 0 ) loc = _priority[ID]; // Y = _order[loc]
//	std::fill(logQxVec.begin(), logQxVec.begin()+loc, 0.0); // for safety

	VarSet done;
	// add all conditioned variables
	for (int b = loc; b < nvar(); ++b) {
		done |= var(_order[b]);
	}
// debug
	assert( _nodes.size() == nvar() );
//	double logPx = 0.0, logQx = 0.0;
//	for (int b = _nodes.size() - 1; b >= 0; --b) {
	for (int b = loc - 1; b >= 0; --b) {
		Factor qMix(0.0);
		for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
			Factor qComponent(0.0);
			VarSet condv = node(*ni).theta.vars() & done;
			qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
			qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
			for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
				VarSet condv = node(*ci).msgFwd.vars() & done;
				qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
			}
//			qComponent -= qComponent.max(); // this causes a bug
			qComponent *= (1.0 / node(*ni).weight);
			qComponent.exp();
			if (qComponent.sum() == 0.0)
				qComponent += 1.0 / qComponent.nrStates();
			qMix += node(*ni).weight * qComponent;
		}
		qMix /= qMix.sum();
		size_t tmp = qMix.sample();
		config[qMix.vars()[0]] = tmp; //qMix.sample();
		double tmp_val = std::log(qMix[tmp]); //config[ qMix.vars()[0] ] ] );
//		logQx += tmp_val;
		done |= qMix.vars();
		// q(x_i | par(x_i))
		logQxVec[b+1] = tmp_val;
	}
	// q(x_0,...,x_{i-1} | x_i, ... x_{nvar-1} ) = q(x_0 | par(x_0)) *...* q(x_{i-1} | par(x_{i-1}))
	for (int i=1; i<nvar()+1; ++i) {
//	for (int i=1; i<loc+1; ++i) {
		logQxVec[i] += logQxVec[i-1]; // logQxVec[0] = 0.0
	}

/*
	// Compute condBD
	for (int i = _order.size() - 1; i >= 0; --i) {
		auto v = _order[i]; // v = _order[i] is the i-th variable eliminated
		// priority[v] is the step in which variable v is eliminated
		size_t b = _priority[v];	// find bucket to examine
		assert(i==int(b)); // debug

		// heuristicTheta are original (re-parameterized) theta in the bucket of the current node
		double heurTheta = heuristicTheta(var(v), config); // pass argument by value, won't change config;
		// heuristicIn are messages from descendants (pass into and pass by) the bucket of the node
		double heurIn = heuristicIn(var(v), config);  // pass argument by value, won't change config;

		thetaVec[i] = heurTheta + thetaVec[i+1]; // thetaVec[_order.size()] = 0.0;
		condBD[i] = heurIn + thetaVec[i];
	}

	// compute (upper) bound on logZ
	auto v = _order.back();
	Var X = var(v);
	mex::vector<uint32_t> tuple(nvar());
	Factor db(X,0.0);
	for (size_t i=0; i<X.states(); ++i) {
		tuple[ X ] = i;
		db[i] = heuristicTheta(X,tuple) +  heuristicIn(X,tuple);
	}
	condBD.back() = db.logsumexp(); // wmb bound on logZ
*/

	// logPx
	return logP(config);
//	condBD[0] = logPx;
	// check, these  quantities should be the same
/*	assert( condBD[0]<logPx+1e-5 && condBD[0]>logPx-1e-5);
	assert( condBD[0]<thetaVec[0]+1e-5 && condBD[0]>thetaVec[0]-1e-5);*/

//	return logQxVec;
}

/*
 * add by Qi Lou on 10/18/2016
 * AND/OR version of sampleMixtureConditioned (replace it eventually?)
 * sample from the conditional mixture distribution.
*/
template <class MapType>
void wmbe::conditionalMixtureSampling(int ID, MapType& config, std::vector<double>& logQxCond, const std::list<Var>& UnDone) {
	/*
	 * support any pseudo tree structure
	 * ID is the unique variable id for the foremost variable conditioned in the order.
	 * sample ID's descendants on the pseudo tree
	 * ID < 0 means no conditioned variables.
	 * input config is instantiated by conditioned variables
	 * logQxCond[i] = log( q(x_i | par(x_i)) ), par(x_i) is the set of all ancestors on the pseudo tree
	 * UnDone is a list of variables that will be sampled, ancestors prior to descendants in the list
	 */
	// position of the variable in the elimination order
//	int loc = nvar(); // default value means no conditioning, ID = -1
	VarSet done; // the set of all ancestors of UnDone on the pseudo tree
	if ( ID >= 0 ) {
//		loc = _priority[ID]; // Y = _order[loc]
	// add all conditioned variables
//	for (int b = loc; b < nvar(); ++b) {
//		done |= var(_order[b]);
//	}
		done |= var(ID);
		auto par = _parents[ID];
		while( par >= 0 && par < nvar() ) {
			done |= var(par);
			par = _parents[par];
		}
	}
//	double logPx = 0.0, logQx = 0.0;
//	for (int b = loc - 1; b >= 0; --b) {

	// UnDone is an ordered list
	for (const auto& X : UnDone) {
	    int b = _priority[X];
		Factor qMix(0.0);
		for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
			Factor qComponent(0.0);
			VarSet condv = node(*ni).theta.vars() & done;
			qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
			qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
			for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
				VarSet condv = node(*ci).msgFwd.vars() & done;
				qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
			}
//			qComponent -= qComponent.max(); // this causes a bug
			qComponent *= (1.0 / node(*ni).weight);
			qComponent.exp();
			if (qComponent.sum() == 0.0)
				qComponent += 1.0 / qComponent.nrStates();
			qMix += node(*ni).weight * qComponent;
		}
		qMix /= qMix.sum();
		size_t tmp = qMix.sample();
		config[qMix.vars()[0]] = tmp; //qMix.sample();
		double tmp_val = std::log(qMix[tmp]); //config[ qMix.vars()[0] ] ] );
//		logQx += tmp_val;
		done |= qMix.vars();
		// q(x_i | par(x_i))
//		logQxVec[b+1] = tmp_val;
		logQxCond[X.label()] = tmp_val;
	}
	// q(x_0,...,x_{i-1} | x_i, ... x_{nvar-1} ) = q(x_0 | par(x_0)) *...* q(x_{i-1} | par(x_{i-1}))
//	for (int i=1; i<nvar()+1; ++i) {
//		logQxVec[i] += logQxVec[i-1]; // logQxVec[0] = 0.0
//	}
	// logPx
//	return logP(config);
}



template <class MapType>
double wmbe::conditionalMixtureSampling(const int ID, MapType& config, const std::list<vindex>& rts,
		const std::vector<std::list<vindex> >& ascList, const std::vector<int>& varDepthVec, const int maxDepth) {
	// 11/14/2017: this is a customized version mainly for "stochastic lookahead"
	// ascList: associated list
	// varDepthVec: depth of each variable in the pseudo tree
	// maxDepth: max depth to reach when sampling
	/*
	 * support any pseudo tree structure
	 * ID is the unique variable id for the foremost variable conditioned in the order.
	 * sample ID's descendants on the pseudo tree
	 * ID < 0 means no conditioned variables.
	 * input config is instantiated by conditioned variables
	 * logQxCond[i] = log( q(x_i | par(x_i)) ), par(x_i) is the set of all ancestors on the pseudo tree
	 * UnDone is a list of variables that will be sampled, ancestors prior to descendants in the list
	 */
	VarSet done; // the set of all ancestors of UnDone on the pseudo tree
	if ( ID >= 0 ) {
		done |= var(ID);
		auto par = _parents[ID];
		while( par >= 0 && par < nvar() ) {
			done |= var(par);
			par = _parents[par];
		}
	}
	//
	double logFx = 0.0, logQx = 0.0;
	// depth-first search to traverse variables that require to be sampled
	std::stack<vindex> stack;
//	auto childrenList = rts; // if conditioned on nothing
//	if (ID >= 0) childrenList = ascList[ID];
	const auto& childrenList = (ID >= 0)? ascList[ID] : rts;
	for (auto it = childrenList.begin(); it != childrenList.end(); ++it) {
		stack.push(*it);
	}
//	for (int b = loc - 1; b >= 0; --b) {
	while(!stack.empty()) {
		auto id = stack.top(); stack.pop();
	// UnDone is an ordered list
//	for (const auto& X : UnDone) {
	    int b = _priority[id];
		Factor qMix(0.0);
		for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
			Factor qComponent(0.0);
			VarSet condv = node(*ni).theta.vars() & done;
			qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
			qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
			for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
				VarSet condv = node(*ci).msgFwd.vars() & done;
				qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
			}
//			qComponent -= qComponent.max(); // this causes a bug
			qComponent *= (1.0 / node(*ni).weight);
			qComponent.exp();
			if (qComponent.sum() == 0.0)
				qComponent += 1.0 / qComponent.nrStates();
			qMix += node(*ni).weight * qComponent;
		}
		qMix /= qMix.sum(); // TODO: can be removed?
		size_t tmp = qMix.sample();
		config[qMix.vars()[0]] = tmp; //qMix.sample();
		double tmp_val = std::log(qMix[tmp]); //config[ qMix.vars()[0] ] ] );
		logQx += tmp_val;
		done |= qMix.vars();
		// q(x_i | par(x_i))
//		logQxCond[X.label()] = tmp_val;

		auto X = var(id);
		int depth = varDepthVec[id]; // depth of current variable in the pseudo tree
		// compute factors being fully instantiated
		logFx += heuristicTheta(X, config);
		if (depth == maxDepth) {
			// compute heuristics for tip variables, should be 0.0 for leaf variables
			logFx += heuristicIn(X, config);
		} else {
			assert(depth < maxDepth);
		// sample within a depth limit
//		if (depth < maxDepth) {
			for (auto it = ascList[id].begin(); it != ascList[id].end(); ++it) stack.push(*it);
		}
	}
	// logPx
//	return logP(config);
	return (logFx - logQx);
}

template <class MapType>
void wmbe::conditionalMixtureSampling(const int ID, MapType& config, const int K,
		const std::list<Var>& unDoneMax, const std::list<Var>& unDoneSum,
		double& logFx, double& logQx, double& aveEst, double& prodEst) {
	// 1/13/18: this is a customized version mainly for "mmapIS"
	/*
	 * support any pseudo tree structure
	 * ID is the unique variable id for the foremost variable conditioned in the order.
	 * sample ID's descendants on the pseudo tree
	 * ID < 0 means no conditioned variables.
	 * input config is instantiated by conditioned variables
	 * logQxCond[i] = log( q(x_i | par(x_i)) ), par(x_i) is the set of all ancestors on the pseudo tree
	 * UnDone is a list of variables that will be sampled, ancestors prior to descendants in the list
	 *
	 * K: no. of samples to draw for SUM variables, can be different from "_K"!!!
	 *
	 * logFx: thetas instantiated by MAX variables
	 * logQx: probability of this partial MAP configuration to be sampled
	 * aveEst: arithmetic mean of K importance weights for the SUM subproblem
	 * prodEst: product of K importance weights for the SUM subproblem
	 *
	 * TODO: when reach a zero-probability configuration, generate dummy descendants to save time? e.g., fixed or uniformly random?
	 */

	// TODO: be more efficient?
	VarSet done; // the set of all ancestors of UnDone on the pseudo tree
	if ( ID >= 0 ) {
		done |= var(ID);
		auto par = _parents[ID];
		while( par >= 0 && par < nvar() ) {
			done |= var(par);
			par = _parents[par];
		}
	}

	const double eps = 1e-5;
	// re-set to 0.0;
	logFx = 0.0, logQx = 0.0, aveEst = 0.0, prodEst = 0.0;
	/**********************************************************************/
	// unDoneMax, unDoneSum are ordered
	for (const auto& X : unDoneMax) {
		// debug
//		assert(_maxVars.contains(X));
		auto id = X.label();
	    int b = _priority[id];
		Factor qMix(0.0);
		for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
			Factor qComponent(0.0);
			VarSet condv = node(*ni).theta.vars() & done;
			qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
			qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
			for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
				VarSet condv = node(*ci).msgFwd.vars() & done;
				// replicate messages from SUM to MAX
				if (node(*ci).isMsgFwdRep){
					qComponent += _K * node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
					// debug
//					std::cout << "message from child " << node(*ci).vid << " replicated " << _K << " times\n";
				} else {
					qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
				}
			}
			qComponent *= (1.0 / node(*ni).weight);
			qComponent.exp();
			if (qComponent.sum() == 0.0)
				qComponent += 1.0 / qComponent.nrStates();
			qMix += node(*ni).weight * qComponent;
			// debug
//			std::cout << "id: " << id << ", weight: " << node(*ni).weight << "\n";
		}
		// double check
		assert (!qMix.isnan());
		double qsum = qMix.sum(); // qsum should be 1
/*		std::cout.precision(10);
		std::cout << "qsum: " << qsum << std::endl;*/
		if (!( qsum < 1+eps && qsum > 1-eps)) {
			std::cout.precision(10);
			std::cout << "qsum: " << qsum << std::endl;
			std::cout << "qMix vals:\n";
			for (unsigned i=0; i<qMix.numel(); ++i) {
				std::cout << qMix[i] << ", ";
			}
			std::cout << std::endl;
			assert(false);
		}
		qMix /= qsum; // no need to do this if everything goes correctly
//		qMix /= qMix.sum();
		size_t tmp = qMix.sample();
		config[qMix.vars()[0]] = tmp; //qMix.sample();
		logQx += std::log(qMix[tmp]); // q(x_i | par(x_i))
		done |= qMix.vars();
		// compute the factors being instantiated
		logFx += heuristicTheta(X, config);
	}

	// we can set K=0 to sample MAX nodes only
	if (K > 0) {
		auto preDone = done;
		Factor logDiff(Var(0, K), 0.0);
		// sample K times for SUM variables
		for (int j = 0; j < K; ++j) {
			done = preDone; // re-initialize
			double lq = 0.0, lf = 0.0;
			for (const auto& X : unDoneSum) {
				// debug
//				assert(!_maxVars.contains(X));
				auto id = X.label();
				int b = _priority[id];
				Factor qMix(0.0);
				for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
					Factor qComponent(0.0);
					VarSet condv = node(*ni).theta.vars() & done;
					qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
					qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
					for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
						VarSet condv = node(*ci).msgFwd.vars() & done;
						qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
					}
					qComponent *= (1.0 / node(*ni).weight);
					qComponent.exp();
					if (qComponent.sum() == 0.0)
						qComponent += 1.0 / qComponent.nrStates();
					qMix += node(*ni).weight * qComponent;
				}
				// double check
				assert (!qMix.isnan());
				double qsum = qMix.sum();
				assert(qsum < 1+eps && qsum > 1-eps); // should already be normalized
				qMix /= qsum; // no need to do this if everything goes correctly
//				qMix /= qMix.sum();
				size_t tmp = qMix.sample();
				config[qMix.vars()[0]] = tmp; //qMix.sample();
				lq += std::log(qMix[tmp]); // q(x_i | par(x_i))
				done |= qMix.vars();
				// compute the factors being instantiated
				lf += heuristicTheta(X, config);
			}
			logDiff[j] = lf - lq; // importance weight
			prodEst += logDiff[j];
		}
		aveEst = logDiff.logsumexp() - std::log(K);
	}
}
//
template <class MapType>
void wmbe::conditionalMixtureSampling(const int ID, MapType& config, const int K, const int start_loc, const int end_loc,
		const std::vector<mex::Var>& pstree, const std::vector<bool>& varTypeVec, const bool MAX_VAR,
		double& logFx, double& logQx, double& aveEst, double& prodEst) {
	// 2/21/18: this is a faster version mainly for "mmapIS"
	/*
	 * this function samples all ID's descendants on the pseudo tree
	 * support any pseudo tree structure
	 * ID: the unique variable id for the foremost variable conditioned in the order.
	 * 		ID < 0 means no conditioned variables.
	 * config: an instantiation of conditioned variables
	 * K: no. of samples to draw for SUM variables, can be different from "_K"!!!
	 * pstree: see "_pseudotree" in mmapIS.h
	 * start_loc: location of ID's corresponding variable in "pstree"
	 * end_loc: "end_loc-1" is the last descendant of "ID" in "pstree"
	 * MAX_VAR: a constant for max variables
	 * logFx: thetas instantiated by MAX variables
	 * logQx: probability of this partial MAP configuration to be sampled
	 * aveEst: arithmetic mean of K importance weights for the SUM subproblem
	 * prodEst: product of K importance weights for the SUM subproblem
	 *
	 * TODO: when reach a zero-probability configuration, generate dummy descendants to save time? e.g., fixed or uniformly random?
	 */

	// TODO: be more efficient?
	VarSet done; // the set of all ancestors of UnDone on the pseudo tree
	bool varType = MAX_VAR; // initialized to max
	if ( ID >= 0 ) {
		done |= var(ID);
		auto par = _parents[ID];
		while( par >= 0 && par < nvar() ) {
			done |= var(par);
			par = _parents[par];
		}
		varType = varTypeVec[ID];
	}

	const double eps = 1e-5;
	// re-set to 0.0
	logFx = 0.0, logQx = 0.0, aveEst = 0.0, prodEst = 0.0;
	Factor logDiff(Var(0, K), 0.0);
	/**********************************************************************/
	int loc = start_loc + 1; // sample from its first descendant!!!
	while (loc < end_loc) {
		auto X = pstree[loc];
		auto id = X.label();
		varType = varTypeVec[id];
		if (varType == MAX_VAR) {
			// debug
//			assert(_maxVars.contains(X));
			// sample this MAX variable once
			int b = _priority[id];
			Factor qMix(0.0);
			for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
				Factor qComponent(0.0);
				VarSet condv = node(*ni).theta.vars() & done;
				qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
				qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
				for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
					VarSet condv = node(*ci).msgFwd.vars() & done;
					// replicate messages from SUM to MAX
					if (node(*ci).isMsgFwdRep) {
						qComponent += _K * node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
					} else {
						qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
					}
				}
				qComponent *= (1.0 / node(*ni).weight);
				qComponent.exp();
				if (qComponent.sum() == 0.0)
					qComponent += 1.0 / qComponent.nrStates();
				qMix += node(*ni).weight * qComponent;
			}
			// double check
			assert(!qMix.isnan());
			double qsum = qMix.sum(); // qsum should be 1
			if ( qsum > 1 + eps || qsum < 1 - eps ) {
				std::cout.precision(10);
				std::cout << "qsum: " << qsum << std::endl;
				std::cout << "qMix vals:\n";
				for (unsigned i = 0; i < qMix.numel(); ++i) {
					std::cout << qMix[i] << ", ";
				}
				std::cout << std::endl;
				assert(false);
			}
			qMix /= qsum; // no need to do this if everything goes correctly
			//		qMix /= qMix.sum();
			size_t tmp = qMix.sample();
			config[qMix.vars()[0]] = tmp; //qMix.sample();
			logQx += std::log(qMix[tmp]); // q(x_i | par(x_i))
			done |= qMix.vars();
			// compute the factors being instantiated
			logFx += heuristicTheta(X, config);
			//
			++ loc;
		} else {
			// debug
//			assert(!_maxVars.contains(X));
			// if this is a sum variable, sample K times, each time touch all its sum descendants
			// we can set K=0 to sample MAX nodes only
			int cur_loc = loc;
			bool cur_varType = varType;
			for (; loc < end_loc; ++ loc) {
				// find all its SUM descendants
				cur_varType = varTypeVec[pstree[loc].label()];
				// var type changed to MAX
				if (cur_varType != varType) break;
			}
			if (K > 0) {
//				auto preDone = done;
//				Factor logDiff(Var(0, K), 0.0);
				for (int j = 0; j < K; ++j) {
					auto sumDone = done; // re-initialize
					double lq = 0.0, lf = 0.0;
//					for (const auto& X : unDoneSum) {
					for (int hh = cur_loc; hh < loc; ++hh) {
						auto Y = pstree[hh];
						// debug
//						assert(!_maxVars.contains(Y));
//						auto id = X.label();
						int b = _priority[Y.label()];
						Factor qMix(0.0);
						for (nodeList::const_iterator ni = _nodes[b].begin(); ni != _nodes[b].end(); ++ni) {
							Factor qComponent(0.0);
							VarSet condv = node(*ni).theta.vars() & sumDone;
							qComponent += node(*ni).theta.condition(condv, sub2ind(condv, config));
							qComponent -= node(*ni).msgFwd[sub2ind(node(*ni).msgFwd.vars(), config)];
							for (nodeList::iterator ci = node(*ni).children.begin(); ci != node(*ni).children.end(); ++ci) {
								VarSet condv = node(*ci).msgFwd.vars() & sumDone;
								qComponent += node(*ci).msgFwd.condition(condv, sub2ind(condv, config));
							}
							qComponent *= (1.0 / node(*ni).weight);
							qComponent.exp();
							if (qComponent.sum() == 0.0)
								qComponent += 1.0 / qComponent.nrStates();
							qMix += node(*ni).weight * qComponent;
						}
						// double check
						assert (!qMix.isnan());
						double qsum = qMix.sum();
						assert(qsum < 1+eps && qsum > 1-eps); // should already be normalized
						qMix /= qsum; // no need to do this if everything goes correctly
		//				qMix /= qMix.sum();
						size_t tmp = qMix.sample();
						config[qMix.vars()[0]] = tmp; //qMix.sample();
						lq += std::log(qMix[tmp]); // q(x_i | par(x_i))
						sumDone |= qMix.vars();
						// compute the factors being instantiated
						lf += heuristicTheta(Y, config);
					}
//					logDiff[j] = lf - lq; // importance weight
					logDiff[j] += lf - lq;
					prodEst += lf - lq;
				}
//				aveEst = logDiff.logsumexp() - std::log(K);
			}
		}
	}
	// do this at the end
//	for (auto& v : logDiff) prodEst += v;
	if (K>0) aveEst = logDiff.logsumexp() - std::log(K);
}

//////////////////////////////////////////////////////////////////////////////////////////////
}       // namespace mex
#endif  // re-include



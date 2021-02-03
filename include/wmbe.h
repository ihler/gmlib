#ifndef __MEX_MINIBUCKET_H
#define __MEX_MINIBUCKET_H

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>

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

  void   setIBound(size_t i) { _iBound=i ? i : std::numeric_limits<size_t>::max(); }
  size_t getIBound() const   { return _iBound; }
  void   setSBound(size_t s) { _sBound=s ? s : std::numeric_limits<size_t>::max(); }
  size_t getSBound() const   { return _sBound; }

  const vector<vindex>& getPseudotree()   const                { return _parents; }
  void                  setPseudotree(const vector<vindex>& p) { _parents = p;    }

  // Build mini-bucket data structure
  void init(const VarSet& vs) { init(); }            // (meaning here?)
  void init();
  virtual void run() { }  // TODO: add

  void build();
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
  void   msgBackward(bucketID b, double dampTheta, double stepWeights);

  void   reparameterize();

  //template <class MapType>
  //double logHeurToGo(Var v, MapType vals) const;  // evaluate remaining cost given context

  template <class MapType>
  double heuristicTheta(Var v, MapType vals) const;   // evaluate heuristic costs for a bucket
  template <class MapType>
  double heuristicIn(Var v, MapType vals) const;  //  "" ""

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




//////////////////////////////////////////////////////////////////////////////////////////////
}       // namespace mex
#endif  // re-include



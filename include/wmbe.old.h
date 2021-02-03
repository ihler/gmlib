#ifndef __MEX_MINIBUCKET_H
#define __MEX_MINIBUCKET_H

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>

#include "factorgraph.h"
#include "alg.h"
#include "mplp.h"


namespace mex {

// Graphical Model Algorithm: Weighted Mini-bucket Elimination 
// 

class wmbe : public graphModel, public gmAlg, virtual public mxObject {

// TODO: resolve "containing" gmo, containing theta, and  "being" a gmodel

public:
  typedef graphModel::findex        findex;        // factor index
  typedef graphModel::vindex        vindex;        // variable index
  typedef graphModel::flist            flist;      // collection of factor indices

	typedef size_t nodeID;
	typedef size_t bucketID;

private:
// Contains:
VarOrder       _order;
vector<size_t> _priority;
vector<double> _varWeight;			// TODO: need to keep?

// Nodes (as collections of vectors :( TODO )
vector<VarSet> _clique;
vector<Factor> _theta;
//vector<factor> belief;
vector<findex> _parent;
vector<vector<findex> > _children;
vector<double> _weight;
vector<Factor> _msgFwd;
vector<Factor> _msgBwd;
vector<nodeID> _nodesVacant;

vector<vector<findex> > _nodes;   // 
//vector<vector<VarSet> > _match;   // TODO ???
vector<vector<vector<findex> > > _match;

vector<vindex> _parents;  // pseudotree data structure  (TODO use?)

double _lnZ;
size_t _ibound;

double _obj; 
double _dampTheta;   // damping factor (1.0 = no damping; 0.0 = no updates)
double _dampWeight;  // step size for weight updates


public:
  wmbe() : graphModel() { setProperties(); }
  wmbe(const graphModel& gm) : graphModel(gm), _gmo(gm)  { clearFactors(); setProperties(); }
  virtual wmbe* clone() const { wmbe* gm = new wmbe(*this); return gm; }

  //void            setOrder(const VarOrder& order);
  void            setOrder(const VarOrder& order, const vector<double>& weights=vector<double>());
  void            setOrder(graphModel::OrderMethod method, const vector<double>& weights=vector<double>());
  // TODO: how to express constrained elimination orders???
  void            setWeight(size_t vIndex, double weight);
  const VarOrder& getOrder()                    { return _order; }

  void   setIBound(size_t i) { _iBound=i ? i : std::numeric_limits<size_t>::max(); }
  size_t getIBound() const   { return _iBound; }
  void   setSBound(size_t s) { _sBound=s ? s : std::numeric_limits<size_t>::max(); }
  size_t getSBound() const   { return _sBound; }

  const vector<vindex>& getPseudotree() { return _parents; }
  void                  setPseudotree(const vector<vindex>& p) { _parents=p; }

  // TODO: what to do about internal graphmodel?
  void setModel( const graphModel& gm ) { _gmo = gm; }
  void setModel( const vector<Factor>& fs ) { _gmo = graphModel(fs); }

  // Build mini-bucket data structure
  void init(const VarSet& vs) { init(); }            // !!! inefficient
  void init();

  // Message passing updates 
  void msgForward(double dampTheta, double stepWeights);
  void msgBackward(); 

  // Simulate build process & evaluate memory requirements
  size_t simulateMemory(vector<VarSet>* cliques = NULL , const VarSet* cond=NULL, size_t MemCutoff = std::numeric_limits<size_t>::max(), bool* isExact=NULL);

  template <class MapType>
  double logHeurToGo(Var v, MapType vals) const;  // evaluate remaining cost given context

  

protected:
  void _setNodeWeight(findex node, double weight);
  graphModel _gmo;


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
  bool mxCheckValid(const mxArray* GM) { 
    return graphModel::mxCheckValid(GM);              // we must be a graphmodel
    // TODO !!! check algorithm elements
  }
  void mxSet(mxArray* GM) {
    if (!mxCheckValid(GM)) throw std::runtime_error("incompatible Matlab object type in mbe");
    graphModel::mxSet(GM);                            // initialize base components  
    // Check for algorithmic specialization
    const mxArray* m_alg = mxGetField(GM,0,"Alg");
    if (m_alg == NULL) throw std::runtime_error("Object does not appear to a wmbe algorithm");  // !!!
    _order.mxSet( (mxArray*) mxGetField(m_alg,0,"order") );		
    _clique.mxSet( (mxArray*) mxGetField(m_alg,0,"clique") );
    _theta.mxSet( (mxArray*) mxGetField(m_alg,0,"theta") );
    _parent.mxSet( (mxArray*) mxGetField(m_alg,0,"parent") );	
    _children.mxSet( (mxArray*) mxGetField(m_alg,0,"children") );
    _weight.mxSet( (mxArray*) mxGetField(m_alg,0,"wt") );
    _msgFwd.mxSet( (mxArray*) mxGetField(m_alg,0,"msgFwd") );
    _msgBwd.mxSet( (mxArray*) mxGetField(m_alg,0,"msgBwd") );
		// TODO: vacancy stack

    _nodes.mxSet( (mxArray*) mxGetField(m_alg,0,"nodes") );	
    _match.mxSet( (mxArray*) mxGetField(m_alg,0,"match") );	
  }
  mxArray* mxGet() { if (!mxAvail()) { graphModel::mxGet(); } return M_; }
  void mxRelease() { throw std::runtime_error("NOT IMPLEMENTED"); }
  void mxDestroy() { throw std::runtime_error("NOT IMPLEMENTED"); }
  void mxSwap(wmbe& gm) { graphModel::mxSwap( (graphModel&)gm ); }
#endif

  // Can be an optimization algorithm or a summation algorithm....
  double ub() const { assert(elimOp==ElimOp::MaxUpper); return _logZ; }
  double lb() const { return _lb; }
  vector<index> best() const { throw std::runtime_error("Not implemented"); }

  double logZ()   const { return _logZ; }
  double logZub() const { assert(elimOp==ElimOp::SumUpper); return _logZ; }
  double logZlb() const { assert(elimOp==ElimOp::SumLower); return _logZ; }

  // No beliefs defined currently
  const Factor& belief(size_t f)  const { throw std::runtime_error("Not implemented"); }
  const Factor& belief(Var v)     const { throw std::runtime_error("Not implemented"); }
  const Factor& belief(VarSet vs) const { throw std::runtime_error("Not implemented"); }
  const vector<Factor>& beliefs() const { throw std::runtime_error("Not implemented"); }

  const graphModel& gmOrig() const { return _gmo; }

  void build(const graphModel& gmo, size_t iBound, const VarOrder& elimOrder);
  virtual void run() { }  // !!! init? or run?


  //MEX_ENUM( ElimOp , MaxUpper,SumUpper,SumLower );

  MEX_ENUM( Property , iBound,sBound,Order,Distance,DoMatch,DoMplp,DoFill,DoJG,DoHeur );

  bool    _byScope;
  bool    _doMatch;
  bool    _doFill;
  bool    _doJG;
  bool    _doHeur;
  Factor::Distance distMethod;
  graphModel::OrderMethod ordMethod;
  size_t   _iBound;
  size_t   _sBound;
  double   _logZ;
  double   _lb;
  vector<flist> atElim;
  vector<double> atElimNorm;

  /////////////////////////////////////////////////////////////////
  // Setting properties (directly or through property string)
  /////////////////////////////////////////////////////////////////


  virtual void setProperties(std::string opt=std::string()) {
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
        //case Property::DoMplp:  _doMplp  = atol(asgn[1].c_str()); break;
        case Property::DoMatch: _doMatch = atol(asgn[1].c_str()); break;
        case Property::DoFill:  _doFill = atol(asgn[1].c_str()); break;
        case Property::DoJG:    _doJG = atol(asgn[1].c_str()); break;
        case Property::DoHeur:  _doHeur = atol(asgn[1].c_str()); break;
        default: break;
      }
    }
  }

  // Scoring function for bucket aggregation
  //   Unable to combine => -3; Scope-only => 1.0; otherwise a positive double score
  double score(const vector<Factor>& fin, const Var& VX, size_t i, size_t j, const vector<Factor>& tmp);
  double scoreSim(const vector<VarSet>& fin, const Var& VX, size_t i, size_t j);

  size_t SizeOf( const vector<Factor>& fs ) {
    size_t MemUsed=0; for (size_t f=0;f<fs.size();++f) MemUsed+=fs[f].nrStates(); return MemUsed;
  }
  size_t SizeOf( const vector<VarSet>& fs ) {
    size_t MemUsed=0; for (size_t f=0;f<fs.size();++f) MemUsed+=fs[f].nrStates(); return MemUsed;
  }

  // helper class for pairs of sorted indices
  struct sPair : public std::pair<size_t,size_t> {
    sPair(size_t ii, size_t jj) { if (ii<jj) {first=jj; second=ii; } else {first=ii; second=jj;} }
  };


  size_t _iter;
  double _dObj,_dObjCurr;
  double _bound;

  double iter()          { return ((double)_iter)/edges().size(); }
  bool   iter_boundary() { return (_iter % edges().size()) == 0; }
  double dObj()          { return std::abs(_dObj); }

  Factor nodeBelief(findex n) {

  }

  Factor computeNodeBelief(nodeID n) const;
  void   updateTheta(bucketID i, double damp, const vector<findex>& match, vector<Factor>& bel);
  void   updateWeights(bucketID i, vector<Factor>& bel);
  void   msgForward(bucketID i, double dampTheta, double stepWeights);
  void   msgBackward(bucketID i, double dampTheta, double stepWeights);

  void   reparameterize();

  // Simulate Memory()

private:
  nodeID helperAddCliqueBasic(VarSet c);		// Add clique c to buckets, but doesn't add parents
  void   helperAddCliqueParent(Var x, nodeID a); 						// add / connect to parent of node a
  nodeID helperAddClique(VarSet c);					// Add clique c and parents on down (?)
  nodeID helperMergeCliques(nodeID a, nodeID b); 	// merge two nodes in the joingraph
  


};


//////////////////////////////////////////////////////////////////////////////////////////////
}       // namespace mex
#endif  // re-include



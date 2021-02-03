#ifndef __MEX_pbp_H
#define __MEX_pbp_H

#define USE_LOG

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdarg>
#include <cstring>
#include <cmath>

#include <fstream>
#include <iostream>

#include "factorgraph.h"
#include "alg.h"
#include "indexedHeap.h"

/*
*/

namespace mex {

// // Particle type represents a point in a multidimensional, continuous domain.
// // Since dimensionality can vary, allocate storage for up to 4 values. This is
// // inefficient for particles of lower dimensionality, but doing this rather than
// // using a vector avoids a level of heap allocation / pointer dereferencing.
// typedef struct Particle {
  // Particle() : dims(0) {}
  // Particle(double a) : x0(a), dims(1) {}
  // Particle(double x0, double x1) : x0(x0), x1(x1), dims(1) {}
  // Particle(double x0, double x1, double x2) : x0(x0), x1(x1), x2(x2), dims(3) {}
  // Particle(double x0, double x1, double x2, double x3) : x0(x0), x1(x1), x2(x2), x3(x3), dims(4) {}
  // Particle& operator= (Particle const& rhs) {
    // if (this != &rhs) {
      // x0=rhs.x0; x1=rhs.x1; x2=rhs.x2; x3=rhs.x3; dims=rhs.dims;
    // }
    // return *this;
  // }
  // double x0, x1, x2, x3;
  // size_t dims;
  // double& operator[](size_t i) {
    // return (i == 0 ? x0 :
           // (i == 1 ? x1 :
           // (i == 2 ? x2 : x3)));
  // }
// } Particle;
typedef vector<double> Particle;
typedef vector<Particle> ParticleDomain;

// Factor graph algorithm specialization for loopy belief propagation
// 

class pbp : public gmAlg, public factorGraph, virtual public mxObject {
public:
  typedef factorGraph::findex        findex;        // factor index
  typedef factorGraph::vindex        vindex;        // variable index
  typedef factorGraph::flist         flist;         // collection of factor indices

public:
  /// Constructors : from nothing, copy, list of factors, or input iterators
  pbp()                                 : factorGraph() { setProperties(); }
  // pbp(const factorGraph& fg)            : factorGraph(fg) { setProperties(); }
  pbp(vector<Factor> fs, vector<void *> fn, vector<vector<double> > fn_params, vector<ParticleDomain> const& p, vector<Factor> const& w) : factorGraph(fs), _p(p), _w(w), _fn(fn), _fn_params(fn_params) { setProperties(); }
  // template <class InputIterator>
  // pbp(InputIterator f, InputIterator l) : factorGraph(f,l) { setProperties(); }

  virtual pbp* clone() const            { pbp* fg = new pbp(*this); return fg; }

#ifdef MEX  
  // MEX Class Wrapper Functions //////////////////////////////////////////////////////////
  //void        mxInit();            // initialize any mex-related data structures, etc (private)
  //inline bool mxAvail() const;     // check for available matlab object
  bool        mxCheckValid(const mxArray*);   // check if matlab object is compatible with this object
  void        mxSet(mxArray*);     // associate with A by reference to data
  mxArray*    mxGet();             // get a pointer to the matlab object wrapper (creating if required)
  void        mxRelease();         // disassociate with a matlab object wrapper, if we have one
  void        mxDestroy();         // disassociate and delete matlab object
  void        mxSwap(pbp& gm);     // exchange with another object and matlab identity
  /////////////////////////////////////////////////////////////////////////////////////////
#endif

  Factor& belief(size_t f) { return _beliefs[f]; }  //!!! const
  Factor belief(size_t f, ParticleDomain const& particles) {
    // Recompute incoming messages for factor f at the specified particle locations.
    // NOTE: this only works when f is a "local factor."
    const set<EdgeID>& nbrs = neighbors(f);
    vector<Factor> incoming(nbrs.size());
    size_t j = 0;
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i,++j) {
      computeMsg(i->reverse(), incoming[j], particles);
    }
    // Combine messages with local factor to get the belief.
    Var v = belief(f).vars()[0];
    VarSet vs = (belief(f).vars() - v) + Var(v, particles.size());
    Factor bel(vs);
    computeFactor(f, bel, v, particles);
#ifdef USE_LOG
    bel.log();
#endif
    for (vector<Factor>::const_iterator i=incoming.begin(); i!=incoming.end(); ++i) {
#ifdef USE_LOG
      bel += *i;
#else
      bel *= *i;
      bel /= bel.sum(); // Normalize as we go for stability.
#endif
    }
#ifdef USE_LOG 
    bel -= bel.logsumexp();
    bel.exp();
#endif
    return bel;
  }

  const Factor& belief(size_t f)  const { return _beliefs[f]; }
  const Factor& belief(Var v)     const { return belief(localFactor(v)); }
  const Factor& belief(VarSet vs) const { throw std::runtime_error("Not implemented"); }
  const vector<Factor>& beliefs() const { return _beliefs; }
  const ParticleDomain& particles(vindex vidx) const { return _p[vidx]; }
  const Factor& weights(vindex vidx) const { return _w[vidx]; }
  bool isDiscrete(vindex vidx) const { return _p[vidx].empty(); }

  // Not a bound-producing algorithm but can try to produce a good config
   double lb() const { throw std::runtime_error("Not available"); }
   double ub() const { throw std::runtime_error("Not available"); }
  vector<index> best() const { throw std::runtime_error("Not available"); }

  // Gives an estimate of the partition function, but not a bound
  double logZ()   const { return _lnZ; }
  double logZub() const { throw std::runtime_error("Not available"); }
  double logZlb() const { throw std::runtime_error("Not available"); }
        


  MEX_ENUM( Schedule , Fixed,Random,Flood,Priority );

  MEX_ENUM( Property , Schedule,Distance,MHiter,StopIter,StopObj,StopMsg,StopTime );
  virtual void setProperties(std::string opt=std::string()) {
    if (opt.length()==0) {
      setProperties("Schedule=Priority,Distance=HPM,StopIter=10,StopObj=-1,StopMsg=-1,StopTime=-1");
      return;
    }
    std::vector<std::string> strs = mex::split(opt,',');
    for (size_t i=0;i<strs.size();++i) {
      std::vector<std::string> asgn = mex::split(strs[i],'=');
      switch( Property(asgn[0].c_str()) ) {
        case Property::Schedule:  _sched  = Schedule(asgn[1].c_str()); break;
        case Property::Distance:  _dist   = Factor::Distance(asgn[1].c_str()); break;
        // case Property::MHiter:    setMHIter( strtol(asgn[1].c_str(),NULL) ); break;
        case Property::StopIter:  setStopIter( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopObj:   setStopObj( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopMsg:   setStopMsg( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopTime:  setStopTime( strtod(asgn[1].c_str(),NULL) ); break;
        default: break;
      }
    }
  }

  void setStopIter(double d) { _stopIter = d * nFactors(); }      // stop when d*(# factors) updates have been done
  // void setMHIter(long d) { _mhIter = d; }      // Perform d iterations of M-H sampling when updating particles.
  void setStopObj(double d)  { _stopObj = d;  }                   // stop when objective change is less than d
  void setStopMsg(double d)  { _stopMsg = d;  }                   // stop when message updates sare less than d
  void setStopTime(double d) { _stopTime = d;  }                  // stop after d seconds


  /// Initialize the data structures
  virtual void init(const VarSet& vs) { init(); } // !!! inefficient
  virtual void init() { 
    _iter = 0;
    // Create factors evaluated at the intial particle sets.
    for (size_t f=0;f<nFactors();++f) {
        if (_fn[f]) {
            computeFactor(f, _factors[f]);
        }
    }
    _iter = 0;
    _beliefs=vector<Factor>(_factors);                            // copy initial beliefs from factors
    _msg=vector<Factor>(2*nEdges());                              // initialize messages to the identity
    for (size_t e=0;e<2*nEdges();++e) if (edge(e)!=EdgeID::NO_EDGE) {  // f'n of the right variables
#ifdef USE_LOG
      _msg[e]=Factor( factor(edge(e).first).vars() & factor(edge(e).second).vars(), 0.0 );
#else
      _msg[e]=Factor( factor(edge(e).first).vars() & factor(edge(e).second).vars(), 1.0 );
#endif
    }
    _msgNew=vector<Factor>(_msg);                                 // copy that as "updated" message list

    _lnZ = 0.0;                                                   // compute initial partition f'n estimate
    for (size_t f=0;f<nFactors();++f) {
      belief(f) /= belief(f).sum();                               // normalize the beliefs
      _lnZ += (belief(f)*log0(factor(f))).sum() + objEntropy(f);   // and compute the free energy estimate
    }

    if (_sched==Schedule::Priority) {                             // for priority scheduling
      for (size_t e=0;e<2*nEdges();++e)                           //  initialize all edges to infinity
       if (edge(e)!=EdgeID::NO_EDGE) priority.insert( std::numeric_limits<double>::infinity() , e);
    } else {
      for (size_t f=0;f<nFactors();++f) forder.push_back(f);      // for fixed scheduling, get a default order
    }
  }


// Useful to make iteration / step # a member
// Then later calls can "pick up" in the middle for more time (?)

  /// Run the algorithm until one of the stopping criteria is reached
  virtual void run() {
    double startTime = timeSystem();
    double print = startTime + 1;
    double dObj=infty(), dMsg=infty();                            // initialize termination values
    size_t n=0;  if (forder.size()) n = _iter % forder.size();

    for (; dMsg>=_stopMsg && _iter<_stopIter && std::abs(dObj)>=_stopObj; ) {
      if (_stopTime > 0 && _stopTime <= (timeSystem()-startTime)) break;       // time-out check

      size_t f;                                                   // factor index
      if (_sched==Schedule::Priority) {                           // priority schedule =>
        f=edge(priority.top().second).second;                     //   get next factor for update from queue
        priority.pop();  
      } else {                                                    // fixed schedule =>
        f=forder[n];                                              //   get next factor from list
        if (++n == forder.size()) n=0; 
      }
      std::cout << f << std::endl;

      if (_sched!=Schedule::Flood) {                              // For non-"flood" schedules,
        Factor logF = log0(factor(f));                            // compute new belief and update objective:
        if (_iter % nFactors()==0) dObj=0.0;                      //   dObj measures per round of all factors
        double delta=0.0;
        delta -= (belief(f)*logF).sum() + objEntropy(f);          //   remove old contribution
        acceptIncoming(f);                                        //   accept all messages into factor f
        delta += (belief(f)*logF).sum() + objEntropy(f);          //   re-add new contribution
        _lnZ += delta; dObj += delta;                             //   and update total
      }
      updateOutgoing(f);                                          //   update outgoing messages from factor f

      if (_sched==Schedule::Priority) dMsg=priority.top().first;  // priority schedule => easy to check msg stop
      else if (_stopMsg>0 && n==0) {                              // else check once each time through all factors
        dMsg=0.0;
        for (size_t e=0;e<2*nEdges();++e) dMsg=std::max(dMsg, _msgNew[e].distance(_msg[e], _dist));
      }

      if (_sched==Schedule::Flood && n==0) {                      // for flooding schedules, recalculate all
        dObj = _lnZ; _lnZ = 0.0;                                  //   the beliefs and objective now
        for (size_t f=0;f<nFactors();++f) {
          acceptIncoming(f);
          _lnZ += (belief(f)*log0(factor(f))).sum() + objEntropy(f);
        }
        dObj -= _lnZ;
      }

      if (timeSystem()>print) { print=timeSystem()+1; std::cout<<"iter "<<_iter/nFactors()<<"; lnZ: "<<_lnZ<<"\n"; }
      _iter++;
    }
    std::cout<<"pbp: "<<_iter/nFactors()<<"it, "<<timeSystem()-startTime<<"sec\n";
  }

  void reparameterize() {

    for (size_t f=0;f<nFactors();++f) _factors[f].log();

    for (size_t e=0;e<2*nEdges();++e) {
      if (edge(e)!=EdgeID::NO_EDGE && isVarNode(edge(e).second)) {
#ifdef USE_LOG
        // Don't take message zeros as gospel (otherwise leads to numerical roundoff errors)
        for (size_t i=0;i<_msg[e].nrStates();++i) if (_msg[e][i]<-10) _msg[e][i]=-10;
        _factors[edge(e).first]  -= _msg[e];
        _factors[edge(e).second] += _msg[e];
#else
        Factor lnMsg = log(_msg[e]);
        for (size_t i=0;i<_msg[e].nrStates();++i) if (_msg[e][i]<-10) _msg[e][i]=-10;
        _factors[edge(e).first]  -= lnMsg;
        _factors[edge(e).second] += lnMsg;
#endif
      }
    }

    double Ztot=0.0;
    for (size_t f=0;f<nFactors();++f) { double Z=_factors[f].max(); _factors[f]-=Z; Ztot+=Z; }
    Ztot /= nFactors();
    for (size_t f=0;f<nFactors();++f) { _factors[f]+=Ztot;  _factors[f].exp(); }

    _lnZ = 0.0;
    // reset all messages to 1
  }

protected:  // Contained objects
  vector<ParticleDomain> _p;                       // particle locations.
  vector<Factor> _w;                               // particle weights.
  vector<void *> _fn;                              // function pointers for cont. factors.
  vector<vector<double> > _fn_params;              // parameters of the factor fns.

  vector<Factor>   _beliefs;                       // store calculated messages and beliefs
  vector<Factor>   _msg;
  vector<Factor>   _msgNew;

  indexedHeap      priority;                       // store priority schedule of edges
  vector<findex>   forder;                         // or fixed order of factors

  double           _lnZ;                           // current objective function value

  Schedule         _sched;                         // schedule type
  Factor::Distance _dist;                           // message distance measure for priority
  double           _stopIter, _stopObj, _stopMsg, _stopTime;   // and stopping criteria
  size_t           _iter;

  /*
  void updateMsg(Edge e) {
    _msgNew[e.idx] = (belief(e.src)/_msg[e.rev]).marginal( belief(e.dst).vars() );
  }
  void acceptMsg(Edge e) {
    belief(e.dst) *= _msgNew[e.idx]/_msg[e.idx];    // update belief to reflect new message
    _msg[e.idx]=_msgNew[e.idx];                      // move message into accepted set
  }
  */

  /// Calculate the entropy contribution to the free energy from node n
  double objEntropy(size_t n) {
    double obj = belief(n).entropy();
    if (!isVarNode(n)) {
      VarSet vs=adjacentVars(n);
      for (VarSet::const_iterator i=vs.begin();i!=vs.end();++i)
        obj -= belief(n).marginal(*i).entropy();
    }
    return obj;
  }

  /// Re-calculate the belief at node n from the current incoming messages
  void calcBelief(size_t n) {
    const set<EdgeID>& nbrs = neighbors(n);        // get all incoming edges
#ifdef USE_LOG
    belief(n)=log(factor(n));                          // calculate local factor times messages
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) belief(n) += _msg[i->ridx];
    belief(n) -= belief(n).logsumexp(); belief(n).exp();
#else
    belief(n)=factor(n);                          // calculate local factor times messages
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) belief(n) *= _msg[i->ridx];
    belief(n) /= belief(n).sum();                 // and normalize
#endif
  }

  /// Accept all the incoming messages into node n, and recompute its belief
  void acceptIncoming(size_t  n) {                //
    const set<EdgeID>& nbrs = neighbors(n);        // get the list of neighbors
    double lnZn=0.0;
#ifdef USE_LOG
    belief(n)=log(factor(n));                          //   and start with just the local factor
#else
    belief(n)=factor(n);                          //   and start with just the local factor
#endif
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) {
      _msg[i->ridx] = _msgNew[i->ridx];           // accept each new incoming message
#ifdef USE_LOG
      belief(n) += _msg[i->ridx];                 //   and include it in the belief
#else
      belief(n) *= _msg[i->ridx];                 //   and include it in the belief
      double Zn=belief(n).sum(); belief(n)/=Zn;   // normalize belief as we go for stability
      lnZn+=std::log(Zn); 
#endif
      if (_sched==Schedule::Priority) 
        priority.erase(i->ridx);                  // accepted => remove from priority queue
    } 
#ifdef USE_LOG 
    double Zn=belief(n).logsumexp(); belief(n)-=Zn;   // normalize belief as we go for stability
    belief(n).exp();
#endif
  }

  /// Recompute new messages from node n to its neighbors
  void updateOutgoing(size_t n) {                  //
    const set<EdgeID>& nbrs = neighbors(n);        // get the list of neighbors
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) {
      computeMsg(*i, _msgNew[i->idx]);
      // Update priority in schedule.
#ifdef USE_LOG
      if (_sched==Schedule::Priority)
        priority.insert( exp(_msgNew[i->idx]).distance(exp(_msg[i->idx]),_dist) , i->idx );
#else
      if (_sched==Schedule::Priority)
        priority.insert( _msgNew[i->idx].distance(_msg[i->idx],_dist) , i->idx );
#endif
    }
  }

  // /// Update the particle set for variable vidx by resampling from its belief.
  // void updateParticles(size_t vidx) {
    // // TODO: make mhiter a property for the pbp class.
    // new_p = _p[vidx]; // Copy current ParticleDomain to a candidate set.
    // for (unsigned int iter=0; iter < _mhIter; iter++) {
      // Particle proposal;
      // for (unsigned int pidx=0; pidx < new_p.size(); pidx++) {
        // proposal = mh_sample(new_p[pidx]);
        // density = mh_density(proposal);
        // if () {
          // // Accept the proposal.
          // new_p[pidx] = proposal;
        // }
      // }
      // // particles <- mhstep(particles)
    // }
    // // store the new particles.
    // // update the importance weights (== the belief, assuming we mixed)
    // // update all neighboring particle factors.
    // // store the new incoming messages to the variable.
  // }

  /// Compute new message along edge e.
  void computeMsg(EdgeID const& e, Factor &out) {                  //
    size_t a = e.first;
    size_t b = e.second;
    VarSet sumVars = belief(a).vars() - belief(b).vars();
    Factor weights;
    for (VarSet::const_iterator j=sumVars.begin();j!=sumVars.end();++j) {
      if (!isDiscrete(vindex(*j))) {
        weights *= _w[*j];
      }
    }
#ifdef USE_LOG
    out = (log(belief(a)*weights)-_msg[e.ridx]).logsumexp(sumVars);
    out -= out.max();   // normalize message
#else
    out = (belief(a)*weights/_msg[e.ridx]).sumOut(sumVars);
    out /= out.sum();   // normalize message
#endif
  }


  /// Compute new message along edge e defined over a specified particle set.
  void computeMsg(EdgeID const& e, Factor &out, ParticleDomain const& particles) {                  //
    size_t a = e.first;
    size_t b = e.second;
    size_t dstVar = belief(b).vars()[0];
    VarSet sumVars = belief(a).vars() - belief(b).vars();
    Factor weights;
    for (VarSet::const_iterator j=sumVars.begin();j!=sumVars.end();++j) {
      weights *= _w[*j];
    }
    Factor newFactor((belief(a).vars() - Var(dstVar, 0)) + Var(dstVar, particles.size())); // Correct variable dims.
    computeFactor(a, newFactor, dstVar, particles);
#ifdef USE_LOG
    // After backing out the reverse message and the old factor, marginalize the
    // destination variable -- otherwise, the Factor will not adjust to the
    // correct dimensionality when we add in the newFactor term.
    out = ((log(belief(a)*weights)-_msg[e.ridx]-log(factor(a))).logsumexp(Var(dstVar,0))+log(newFactor)).logsumexp(sumVars);
    out -= out.max();   // normalize message
#else
    out = ((belief(a)*weights/_msg[e.ridx]/factor(a)).sumOut(Var(dstVar,0))*newFactor).sumOut(sumVars);
    out /= out.sum();   // normalize message
#endif
  }
  
  /// Compute new message from node a to b over a specified particle set.
  void computeMsg(size_t a, size_t b, Factor &out, ParticleDomain const& particles) {                  //
    EdgeID e = edge(a, b);
    computeMsg(e, out, particles);
  }

  /// Compute new message from node a to b.
  void computeMsg(size_t a, size_t b, Factor &out) {                  //
    EdgeID e = edge(a, b);
    computeMsg(e, out);
  }

  /// Compute a discrete factor based on the current particle set.
  void computeFactor(size_t fidx, Factor &out) {
    computeFactor(fidx, out, -1, ParticleDomain());
  }

  /// Compute a discrete factor using the provided particle set for variable vidx.
  void computeFactor(size_t fidx, Factor &out, size_t vidx, ParticleDomain const& particles) {
    VarSet const& vars = out.vars();
    const VarSet::vsize* dims = vars.dims();
    ParticleDomain const* pp[out.nvar()];
    unsigned int d;
    for (d=0; d<out.nvar(); ++d) {
      if (vars[d] == vidx) {
        pp[d] = &particles; // !!!
      } else if (!isDiscrete(vindex(vars[d]))) {
        pp[d] = &_p[vars[d]];
      } else {
        // This variable is discrete. Create fake particles that count its values.
        // Don't forget to delete these at the end of the method.
        ParticleDomain *fake = new ParticleDomain();
        for (unsigned int v=0; v<dims[d]; v++) {
            fake->push_back(Particle(1, v));
        }
        pp[d] = fake;
      }
    }
    // ind[i] is the index of the "current" particle for varible vars[i]
    vector<unsigned int> ind(vars.nvar(), 0);
    size_t linearInd = 0;
    vector<Particle> ps(vars.nvar());
    // Initialize ps to hold the 0th particle for each variable.
    for (d=0; d<out.nvar(); ++d) {
      ps[d] = (*pp[d])[0];
    }
    bool done = false;
    while (1) {
      // Evaluate the function at the current set of particles.
      for (d=0; d<out.nvar(); ++d) {
        ps[d] = (*pp[d])[ind[d]]; // Set the particle for the dth variable.
      }
      if (fidx < _fn.size()) {
          out[linearInd++] = ((double (*)(vector<Particle> const&, vector<double> const& params))_fn[fidx])(ps, _fn_params[fidx]);
      } else {
          // This is an internally-created local factor.
          out[linearInd++] = 1;
      }
      // Advance the index to select the next set of particles.
      for (d=0; d<vars.nvar(); ++d) {
        if (ind[d] == dims[d]-1) {
          if (d == vars.nvar()-1) {
            done = true;
            break;
          }
          ind[d] = 0;
        } else {
          ind[d]++;
          break;
        }
      }
      if (done) {
        break;
      }
    }

    // clean up the fake particle domains.
    for (d=0; d<out.nvar(); ++d) {
      if (_p[vars[d]].empty()) {
        delete pp[d];
      }
    }
  }
};


#ifdef MEX
//////////////////////////////////////////////////////////////////////////////////////////////
// MEX specific functions, and non-mex stubs for compatibility
//////////////////////////////////////////////////////////////////////////////////////////////
bool pbp::mxCheckValid(const mxArray* GM) { 
  //if (!strcasecmp(mxGetClassName(GM),"graphModel")) return false;
  // hard to check if we are derived from a graphmodel without just checking elements:
  return factorGraph::mxCheckValid(GM);                      // we must be a factorGraph
  // !!! check that we have beliefs, or can make them
}

void pbp::mxSet(mxArray* GM) {
  if (!mxCheckValid(GM)) throw std::runtime_error("incompatible Matlab object type in factorGraph");
  factorGraph::mxSet(GM);                            // initialize base components  

  // Check for algorithmic specialization???
}

mxArray* pbp::mxGet() {
  if (!mxAvail()) {
    factorGraph::mxGet();
  }
  return M_;
}
void pbp::mxRelease() {  throw std::runtime_error("NOT IMPLEMENTED"); }
void pbp::mxDestroy() { throw std::runtime_error("NOT IMPLEMENTED"); }

void pbp::mxSwap(pbp& gm) {
  factorGraph::mxSwap( (factorGraph&)gm );
}
#endif


//////////////////////////////////////////////////////////////////////////////////////////////
}       // namespace mex
#endif  // re-include

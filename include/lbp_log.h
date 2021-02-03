#ifndef __MEX_LBP_H
#define __MEX_LBP_H

#define USE_LOG

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdarg>
#include <cstring>

#include <fstream>

#include "factorgraph.h"
#include "alg.h"
#include "indexedHeap.h"

/*
*/

namespace mex {

// Factor graph algorithm specialization for loopy belief propagation
// 

class lbp : public gmAlg, public factorGraph, virtual public mxObject {
public:
  typedef factorGraph::findex        findex;        // factor index
  typedef factorGraph::vindex        vindex;        // variable index
  typedef factorGraph::flist         flist;         // collection of factor indices

public:
  /// Constructors : from nothing, copy, list of factors, or input iterators
  lbp()                                 : factorGraph() { setProperties(); }
  lbp(const factorGraph& fg)            : factorGraph(fg) { setProperties(); for (size_t f=0;f<_factors.size();++f) _factors[f].log(); }
  lbp(vector<Factor> fs)                : factorGraph(fs) { setProperties(); for (size_t f=0;f<_factors.size();++f) _factors[f].log(); }
  template <class InputIterator>
  lbp(InputIterator f, InputIterator l) : factorGraph(f,l) { setProperties(); for (size_t f=0;f<_factors.size();++f) _factors[f].log(); }

  virtual lbp* clone() const            { lbp* fg = new lbp(*this); return fg; }

#ifdef MEX  
  // MEX Class Wrapper Functions //////////////////////////////////////////////////////////
  //void        mxInit();            // initialize any mex-related data structures, etc (private)
  //inline bool mxAvail() const;     // check for available matlab object
  bool        mxCheckValid(const mxArray*);   // check if matlab object is compatible with this object
  void        mxSet(mxArray*);     // associate with A by reference to data
  mxArray*    mxGet();             // get a pointer to the matlab object wrapper (creating if required)
  void        mxRelease();         // disassociate with a matlab object wrapper, if we have one
  void        mxDestroy();         // disassociate and delete matlab object
  void        mxSwap(lbp& gm);     // exchange with another object and matlab identity
  /////////////////////////////////////////////////////////////////////////////////////////
#endif

  Factor& belief(size_t f) { return _beliefs[f]; }  //!!! const
  const Factor& belief(size_t f)  const { return _beliefs[f]; }
  const Factor& belief(Var v)     const { return belief(localFactor(v)); }
  const Factor& belief(VarSet vs) const { throw std::runtime_error("Not implemented"); }
  const vector<Factor>& beliefs() const { return _beliefs; }

  // Not a bound-producing algorithm but can try to produce a good config
   double lb() const { throw std::runtime_error("Not available"); }
   double ub() const { throw std::runtime_error("Not available"); }
  vector<index> best() const { throw std::runtime_error("Not available"); }

  // Gives an estimate of the partition function, but not a bound
  double logZ()   const { return _lnZ; }
  double logZub() const { throw std::runtime_error("Not available"); }
  double logZlb() const { throw std::runtime_error("Not available"); }
        


  MEX_ENUM( Schedule , Fixed,Random,Flood,Priority );

  MEX_ENUM( Property , Schedule,Distance,StopIter,StopObj,StopMsg,StopTime );
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
        case Property::StopIter:  setStopIter( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopObj:   setStopObj( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopMsg:   setStopMsg( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopTime:  setStopTime( strtod(asgn[1].c_str(),NULL) ); break;
        default: break;
      }
    }
  }

  void setStopIter(double d) { _stopIter = d * nFactors(); }      // stop when d*(# factors) updates have been done
  void setStopObj(double d)  { _stopObj = d;  }                   // stop when objective change is less than d
  void setStopMsg(double d)  { _stopMsg = d;  }                   // stop when message updates sare less than d
  void setStopTime(double d) { _stopTime = d;  }                  // stop after d seconds


  /// Initialize the data structures
  virtual void init(const VarSet& vs) { init(); } // !!! inefficient
  virtual void init() { 
    _iter = 0;
    _beliefs=vector<Factor>(_factors);                            // copy initial beliefs from factors
    _msg=vector<Factor>(2*nEdges());                              // initialize messages to the identity
    for (size_t e=0;e<2*nEdges();++e) if (edge(e)!=EdgeID::NO_EDGE) {  // f'n of the right variables
      _msg[e]=Factor( factor(edge(e).first).vars() & factor(edge(e).second).vars(), 0.0 );
    }
    _msgNew=vector<Factor>(_msg);                                 // copy that as "updated" message list

    _lnZ = 0.0;                                                   // compute initial partition f'n estimate
    for (size_t f=0;f<nFactors();++f) {
      belief(f) -= belief(f).logsumexp();                               // normalize the beliefs
      _lnZ += (exp(belief(f))*log0(factor(f))).sum() + objEntropy(f);   // and compute the free energy estimate
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

      if (_sched!=Schedule::Flood) {                              // For non-"flood" schedules,
        Factor logF = log0(exp(factor(f)));                            // compute new belief and update objective:
        if (_iter % nFactors()==0) dObj=0.0;                      //   dObj measures per round of all factors
        double delta=0.0;
        delta -= (exp(belief(f))*logF).sum() + objEntropy(f);          //   remove old contribution
        acceptIncoming(f);                                        //   accept all messages into factor f
        delta += (exp(belief(f))*logF).sum() + objEntropy(f);          //   re-add new contribution
        _lnZ += delta; dObj += delta;                             //   and update total
      }
      updateOutgoing(f);                                          //   update outgoing messages from factor f

      if (_sched==Schedule::Priority) dMsg=priority.top().first;  // priority schedule => easy to check msg stop
      else if (_stopMsg>0 && n==0) {                              // else check once each time through all factors
        dMsg=0.0;
        for (size_t e=0;e<2*nEdges();++e) dMsg=std::max(dMsg, exp(_msgNew[e]).distance(exp(_msg[e]), _dist));
      }

      if (_sched==Schedule::Flood && n==0) {                      // for flooding schedules, recalculate all
        dObj = _lnZ; _lnZ = 0.0;                                  //   the beliefs and objective now
        for (size_t f=0;f<nFactors();++f) {
          acceptIncoming(f);
          _lnZ += (exp(belief(f))*log0(exp(factor(f)))).sum() + objEntropy(f);
        }
        dObj -= _lnZ;
      }

      if (timeSystem()>print) { print=timeSystem()+1; std::cout<<"iter "<<_iter/nFactors()<<"; lnZ: "<<_lnZ<<"\n"; }
      _iter++;
    }
    std::cout<<"LBP: "<<_iter/nFactors()<<"it, "<<timeSystem()-startTime<<"sec\n";
  }

  void reparameterize() {
    _lnZ = 0.0;

    //!!! for (size_t f=0;f<nFactors();++f) _factors[f].log();
    for (size_t e=0;e<2*nEdges();++e) {
      if (edge(e)!=EdgeID::NO_EDGE && isVarNode(edge(e).second)) {
        _factors[edge(e).first]  -= _msg[e];
        _factors[edge(e).second] += _msg[e];
      }
    }
    //double Ztot=0.0;
    //for (size_t f=0;f<nFactors();++f) { double Z=_factors[f].max(); _factors[f]-=Z; Ztot+=Z; }
    //Ztot /= nFactors();
    //for (size_t f=0;f<nFactors();++f) { _factors[f]+=Ztot; _factors[f].exp(); }

    _lnZ = 0;
    // reset all messages to 1
  }

protected:  // Contained objects
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
    Factor bel = exp(belief(n));
    double obj = bel.entropy();
    if (!isVarNode(n)) {
      VarSet vs=adjacentVars(n);
      for (VarSet::const_iterator i=vs.begin();i!=vs.end();++i)
        obj -= bel.marginal(*i).entropy();
    }
    return obj;
  }

  /// Re-calculate the belief at node n from the current incoming messages
  void calcBelief(size_t n) {
    const set<EdgeID>& nbrs = neighbors(n);        // get all incoming edges
    belief(n)=factor(n); //!!!                          // calculate local factor times messages
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) belief(n) += _msg[i->ridx];
    belief(n) -= belief(n).logsumexp();  // !!!
  }

  /// Accept all the incoming messages into node n, and recompute its belief
  void acceptIncoming(size_t  n) {                //
    const set<EdgeID>& nbrs = neighbors(n);        // get the list of neighbors
    double lnZn=0.0;
    belief(n)=factor(n);  // !!!                        //   and start with just the local factor
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) {
      _msg[i->ridx] = _msgNew[i->ridx];           // accept each new incoming message
      belief(n) += _msg[i->ridx];                 //   and include it in the belief
      if (_sched==Schedule::Priority) 
        priority.erase(i->ridx);                  // accepted => remove from priority queue
    } 
    double Zn=belief(n).logsumexp(); belief(n)-=Zn;   // normalize belief as we go for stability
  }

  /// Recompute new messages from node n to its neighbors
  void updateOutgoing(size_t n) {                  //
    const set<EdgeID>& nbrs = neighbors(n);        // get the list of neighbors
    for (set<EdgeID>::const_iterator i=nbrs.begin();i!=nbrs.end();++i) {
      _msgNew[i->idx] = (log(belief(n))-_msg[i->ridx]).logsumexp( belief(n).vars() - belief(i->second).vars() );
      _msgNew[i->idx] -= _msgNew[i->idx].max();   // normalize message
      if (_sched==Schedule::Priority)             // and update priority in schedule
        priority.insert( exp(_msgNew[i->idx]).distance(exp(_msg[i->idx]),_dist) , i->idx );
    }
  }


};


#ifdef MEX
//////////////////////////////////////////////////////////////////////////////////////////////
// MEX specific functions, and non-mex stubs for compatibility
//////////////////////////////////////////////////////////////////////////////////////////////
bool lbp::mxCheckValid(const mxArray* GM) { 
  //if (!strcasecmp(mxGetClassName(GM),"graphModel")) return false;
  // hard to check if we are derived from a graphmodel without just checking elements:
  return factorGraph::mxCheckValid(GM);                      // we must be a factorGraph
  // !!! check that we have beliefs, or can make them
}

void lbp::mxSet(mxArray* GM) {
  if (!mxCheckValid(GM)) throw std::runtime_error("incompatible Matlab object type in factorGraph");
  factorGraph::mxSet(GM);                            // initialize base components  

  // Check for algorithmic specialization???
}

mxArray* lbp::mxGet() {
  if (!mxAvail()) {
    factorGraph::mxGet();
  }
  return M_;
}
void lbp::mxRelease() {  throw std::runtime_error("NOT IMPLEMENTED"); }
void lbp::mxDestroy() { throw std::runtime_error("NOT IMPLEMENTED"); }

void lbp::mxSwap(lbp& gm) {
  factorGraph::mxSwap( (factorGraph&)gm );
}
#endif


//////////////////////////////////////////////////////////////////////////////////////////////
}       // namespace mex
#endif  // re-include

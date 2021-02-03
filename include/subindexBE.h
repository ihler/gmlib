#ifndef __MEX_SUBINDEX_H
#define __MEX_SUBINDEX_H

#include <iostream>

namespace mex {

// "Big-endian" version (first variable increment is non-consecutive memory)

/*****************************************************************************
Some basics on table-based factor representations:
Each entry in a table corresponds implicitly to some configuration, or tuple,
of the argument variables.  The order in which these tuples occur is critical.
Assuming the factor is F(x1,x2,x3,x4) and arguments are in "sorted" order,
we can imagine two standard orders we could use:
  Little Endian (Matlab)                   Big Endian (UAI)
  0000                                     0000
  1000                                     0001
  0100                                     0010
  1100                                     0011
  ...                                      ...
In other words, in LE format, the span of the first variable (x1) is small,
meaning that increments to x1's value correspond to small changes in the
accessed memory location, while the span of x4 is large (incrementing x4
corresponds to a much larger memory location).

We usually store our tables in Little Endian format, for Matlab compatibility.
(But, this header file contains an alternative, big-endian implementation!)

*****************************************************************************
Basic operations:

The "subindex" class allows one to walk through a function with a large scope,
say F(x1,x2,x3), and use a subindex to simultaneously step through the 
corresponding elements of a function with smaller scope, e.g. G(x1,x3).

 Ex:  sub = [x1,x3]    =>  [0 0] [0 1] [0 0] [0 1] [1 0] [1 1] [1 0] [1 1]
     full = [x1,x2,x3] =>  [000] [001] [010] [011] [100] [101] [110] [111]

The "superindex" class allows one to walk through an implicitly specified
sub-function, or a "slice" of F:  F(x1, x2=X, x3)

The "permuteIndex" class allows walking through a function whose arguments
are specified in a different (unsorted) order, and/or has a different 
Endian-ness.

*****************************************************************************/



class subindex {
  public:
   typedef VarSet::vsize vsize;

   vsize  idx_,end_;
   vsize  Nd;
   vsize  *state;           // vector of variable-indices (values) for the current position (1-based)
   const vsize *dims;       // dimensions of each variable
   bool *skipped;           // variables not included in the sub-index
   vsize *add;              // how much to add when we increment
   vsize *subtract;         // how much to subtract when we wrap each variable
   static const End=vsize(-1);
   

   subindex(const VarSet& full, const VarSet& sub) {
     assert( full >> sub );
     idx_=0; 
     end_=1;
     Nd=full.nvar(); 
     dims     = full.dims(); 
     state    = new vsize[Nd];
     add      = new vsize[Nd];
     subtract = new vsize[Nd];
     skipped  = new bool[Nd];
     // Compute reference index updates
     vsize i,j; 
     for (i=full.nvar()-1,j=sub.nvar()-1;i!=End;--i) {
       state[i]=1;                                         // start with [0] = (1,1,1,1...)
       skipped[i] = ( j==End || sub[j]!=full[i] );         // are we sub-indexing this variable?
       if (i==full.nvar()-1) add[i]=1;                     // how much does adding one to this var
       else add[i]=add[i+1]*(skipped[i+1]?1:dims[i+1]);    //    add to our position?
       subtract[i]=add[i]*((skipped[i]?1:dims[i])-1);      // how much does wrapping back to 1 remove from pos?
       if (!skipped[i]) j--;
       end_ *= dims[i];
     } 
   }

   subindex(const subindex& S) {
     dims = S.dims; idx_=S.idx_; end_=S.end_;
     Nd = S.Nd;
     state    = new vsize[Nd]; std::copy(S.state,S.state+Nd,state);
     add      = new vsize[Nd]; std::copy(S.add,S.add+Nd,add);
     subtract = new vsize[Nd]; std::copy(S.subtract,S.subtract+Nd,subtract);
     skipped  = new bool[Nd];  std::copy(S.skipped,S.skipped+Nd,skipped);
   }

   ~subindex(void) {
     delete[] state; delete[] skipped; delete[] add; delete[] subtract;
   }

   subindex& reset() {
     for (size_t i=0;i<Nd;++i) state[i]=1;
     idx_=0;
     return *this;
   };

   size_t end() { return end_; }

   /// Prefix addition operator
   subindex& operator++ (void) {
   for (size_t i=Nd-1;i!=End;--i) {             // for each variable                   
     if (state[i]==dims[i]) {                   // if we reached the maximum, wrap around to 1
       state[i]=1;                               //   subtract wrap value from position
       if (!skipped[i]) idx_ -= subtract[i];     //   and continue to next variable in sequence
     } else {                                   // otherwise, increment variable index
       ++state[i];                               //   add to our current position
       if (!skipped[i]) idx_ += add[i];          //   and break (leave later vars the same)
       break;
     }
   }
   return *this;
  }
  /// Postfix addition : use prefix and copy
  subindex operator++ (int) { subindex S(*this); ++(*this); return S; }

  /// Conversion to index value
  operator size_t() const { return idx_; };

};


class superindex {
public:
  typedef VarSet::vsize vsize;

  vsize  idx_,end_;
  vsize  Ns;
  vsize  *state;           // vector of variable-indices (values) for the current position (1-based)
  const vsize *dims;       // dimensions of each variable
  vsize *add;              // how much to add when we increment
   

  superindex(const VarSet& full, const VarSet& sub, size_t offset) {
    assert( full >> sub );
    idx_=offset; end_=full.nrStates();
    Ns = sub.nvar();
    const vsize* dimf = full.dims(); dims=sub.dims();
    state=new vsize[Ns]; add=new vsize[Ns]; 
    size_t i,j,d=1;
    for (i=0,j=0;j<Ns;++i) {
      if (full[i]==sub[j]) {
        state[j]=1;
        add[j]=d;
        j++;
      }
      d*=dimf[i];
    }
    end_ = add[Ns-1]*dims[Ns-1]+offset;
  }

  superindex(const superindex& S) {
    dims = S.dims; idx_=S.idx_; end_=S.end_;
    Ns = S.Ns;
    state    = new vsize[Ns]; std::copy(S.state,S.state+Ns,state);
    add      = new vsize[Ns]; std::copy(S.add,S.add+Ns,add);
  }

  ~superindex(void) {
    delete[] state; delete[] add;
  }

   superindex& reset() {
     for (size_t i=0;i<Ns;++i) state[i]=1;
     idx_=0;
     return *this;
   };

   size_t end() { return end_; }

   /// Prefix addition operator
   superindex& operator++ (void) {
   for (size_t i=0;i<Ns;++i) {                  // for each variable                   
     if (state[i]==dims[i] && i<Ns-1) {         // if we reached the maximum, wrap around to 1
       state[i]=1;                               // 
      idx_ -= add[i]*(dims[i]-1);
     } else {                                   // otherwise, increment variable index
       ++state[i];                               //   add to our current position
      idx_+= add[i];
       break;
     }
   }
   return *this;
  }
  /// Postfix addition : use prefix and copy
  superindex operator++ (int) { superindex S(*this); ++(*this); return S; }

  /// Conversion to index value
  operator size_t() const { return idx_; };

};



/*
class permute {
  public:
    permute(const VarSet& vars, const vector<size_t> pi) {
      idx_=0;
      end_=0;
      Nd=vars.nvar();
      dims = vars.dims();
      state = new vsize[Nd]; pi_=new vindex[Nd];
      for (size_t i=0;i<Nd;++i) { state[i]=0; pi_[i]=pi[i]; }
    }
    permute& operator++ (void) {

      for (size_t i=0;i<Nd;i++) {
        state[i]++; if (state[i]==


}
}
}
*/

} // end namespace mex
#endif


struct SNode {
  SNode() { };
  SNode(const SNode& s, vsize v) { tuple=s.tuple; tuple.push_back(v); }
  vector<vsize> tuple;
  double hUp, hLo;
  double priority;
  bool less(SNode s) { return priority < s.priority; }
}

double calcPriority(double up, double lo) { return up + log(1-exp(min(lo,up)-up)); }

friend void searchAStar(const wmbe& mbUp, const wmbe& bmLo) {
  size_t nV = mbUp.nvar();    // total number of variables
  vector<size_t> states(nV);  // store var -> state mapping for sub2ind
  size_t count=0;

  std::priority_queue<SNode> queue;  // initialize queue with empty (root) search node
  queue.push( SNode() );
  while (!queue.empty()) {
    SNode best = queue.top();
    size_t depth = best.tuple.size();
    if (depth==nV) continue;   // leaf node?  (numerical roundoff)
    Var x = mbUp.var( mbUp.order()[nV-depth] );
    for (size_t j=0; j<depth; ++j) states[mbUp.order()[nV-1-j]]=best.tuple[j]; //copy tuple to map

    double sumUp = -infty(), sumLo = -infty();
    flist nodesUp = mbUp.minibucket[x].nodes;
    flist nodesLo = mbLo.minibucket[x].nodes;
    for (vsize v=0; v<x.states(); ++v) {
      double hUpi=hUp, hLoi=hLo;
      SNode add; add.tuple=n.tuple; add.tuple.push_back(v);
      for (size_t j=0;j<nodesUp.size();++j) {
        findex n;
        size_t idx;
        n = nodesUp[j];
        idx = mbUp._nodes[n].clique.sub2ind(states);
        hUpi += mbUp._nodes[n].theta[idx];
        idx = sub2ind(mbUp._nodes[n].msgFwd.vars() , states);
        hUpi -= mbUp._nodes[n].msgFwd[idx];
        for (size_t cc=0;cc<mbUp._nodes[n].children.size();++cc) {
          size_t c = mbUp._nodes[n].children[cc];
          idx = sub2ind( mbUp._nodes[c].msgFwd.vars() , states);
          hUpi += mbUp._nodes[c].msgFwd[idx];
        }
        n = nodesLo[j];
        idx = sub2ind( mbLo._nodes[n].clique , states);
        hLoi += mbLo._nodes[n].theta[idx];
        idx = sub2ind( mbLo._nodes[n].msgFwd.vars() , states);
        hLoi -= mbLo._nodes[n].msgFwd[idx];
        for (size_t cc=0;cc<mbLo._nodes[n].children.size();++cc) {
          size_t c = mbLo._nodes[n].children[cc];
          idx = sub2ind( mbLo._nodes[c].msgFwd.vars() , states);
          hLoi += mbLo._nodes[c].msgFwd[idx];
        }
      }
      add.hUp=hUpi; add.hLo=hLoi; add.priority=calcPriority(hUpi,hLoi);
      queue.push(add);
      // update total over x's states "v":
      if (v==0) sumUp=hUpi; else sumUp+=log(1+exp(hUpi-sumUp)); end;
      if (v==0) sumLo=hLoi; else sumLo+=log(1+exp(hLoi-sumLo)); end;
    } // end for over states of x
    lnZUp += log(1-exp(hUp-lnZUp)+exp(sumUp-lnZUp));
    lnZLo += log(1-exp(hLo-lnZLo)+exp(sumLo-lnZLo));
    if (count==0) std::cout<<"(time) "<<lnZUp<<" / "<<lnZLo<<"\n";
    count = (count+1)%1000;
  }

}

/*  CODE TODO:

(1) type check "vindex", "findex", "flist", etc.
    -- ensure compatibility with or distinction from size_t, etc.

(2) 



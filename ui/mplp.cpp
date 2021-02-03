// Wrapper for UAI competition code
// 

#include <cstdio>
#include <iostream>
#include <fstream>
#include "boost/program_options.hpp"

#include "factorgraph.h"

#include "mplp.h"
#include "mbe.h"

#include "lbp.h"
#include "gbp.h"

using namespace std;
using mex::mxObject;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::vector;
using mex::graphModel;
using mex::factorGraph;
using mex::mplp;

using mex::timeSystem;

namespace po = boost::program_options;

#define c_log10 std::log(10)

double MemLimit;
double mplpTime, mplpIter, mplpObj, mplpGap;
double dt;
int nLocal;
double addEpsilon;


//Usage: mplp  -f <file.uai> -e <file.evid> -i ibound
int main(int argc, char* argv[])
{

  double timeStart = timeSystem();

  const char* probName; // = argv[1];
  const char* taskName;//  = argv[3];
  mex::vector<Factor> bel;

  po::options_description desc("Available options");
  desc.add_options()
    ("help", "print help message")
    ("file,f", po::value<std::string>(), "input problem filename")
    ("evidence,e", po::value<std::string>(), "input evidence filename")
    ("seed,S", po::value<int>(),         "random number initial seed")
    ("ibound,i", po::value<int>(),       "initial i-bound")
    ("orders,o", po::value<int>(),       "number of variable orderings to try")
    ("ordertime,t", po::value<double>(), "max time spend on variable orderings")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("dt", po::value<double>(&dt)->default_value(300),         "file write time interval")
    ("mplps", po::value<double>(&mplpTime)->default_value(300),  "mplp stop (seconds)")
    ("mplpi", po::value<double>(&mplpIter)->default_value(20000), "mplp stop (iterations)")
    ("mplpo", po::value<double>(&mplpObj )->default_value(-1),   "mplp stop (objective)")
    ("mplpg", po::value<double>(&mplpGap )->default_value(-1),   "mplp stop (gap)")
    ("local", po::value<int>(&nLocal )->default_value(0),   "local search passes per rounding (default 0)")
    ("eps",   po::value<double>(&addEpsilon )->default_value(0.0), "use f(x)+eps for factors (default 0.0)")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);

  if (vm.count("help")) { std::cout<<desc<<"\n"; return 1; }
  if (vm.count("file")) { probName=vm["file"].as<std::string>().c_str(); }
  else { std::cout<<"Missing input problem file!\n"; return 1; }
  if (vm.count("seed")) { mex::randSeed( vm["seed"].as<int>() ); }


  /*** READ IN PROBLEM FILE **********************************************************/

  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> flist = Factor::readUai10(is);
  size_t nvar=0;
  for (size_t f=0;f<flist.size();++f)                        // find maximum variable label
    nvar=std::max(nvar,(size_t)(flist[f].vars().rbegin()->label()+1));
  for (size_t f=0;f<flist.size();++f) flist[f] += addEpsilon;

  // Read in (single!) evidence 
  VarSet evVar;
  std::map<uint32_t,size_t> evid;
  ifstream is2;
  if (vm.count("evidence")) { is2.open(vm["evidence"].as<std::string>().c_str()); }
  if (is2.is_open()) {
    std::cout<<"Got evidence file\n";
    int nEvid; is2 >> nEvid;
    if (nEvid > 0) {
      int nEvidVar; is2 >> nEvidVar;
      for (size_t i=0;i<nEvidVar;i++) {
        uint32_t vid; size_t vval; is2>>vid>>vval; 
        evid[vid]=vval; evVar |= Var(vid,0);
      }
      for (size_t f=0;f<flist.size();f++) {
        if (flist[f].vars().intersects(evVar)) {
          VarSet overlap = flist[f].vars() & evVar;
          evVar |= overlap; // correct dimensions in evVar
          flist[f] = flist[f].condition( overlap, sub2ind(overlap,evid) );
        }
      }
      for (size_t i=0;i<nEvidVar;++i) flist.push_back( Factor::delta(evVar[i],evid[evVar[i]]) );
    }
  } else std::cout<<"Evidence file not specified or not found\n";

  std::string outfiles(probName); outfiles += '.'; outfiles += "MPE"; // taskName;
  std::string::size_type start = outfiles.find_last_of('/');
  if (start==std::string::npos) start=0; else ++start;
  outfiles = outfiles.substr(start,std::string::npos);
  const char* outfile = outfiles.c_str();
  std::cout<<"Writing to "<<outfile<<"\n";
  remove( outfile );

  double ln10 = std::log(10);

  /*** PERFORM REQUESTED TASK ********************************************************/
  mex::factorGraph fg(flist);

  size_t ibound = 18, InducedWidth=10000;
  mex::VarOrder order;
  double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
  double mbMemRatio = 1.0*sizeof(double)/1024.0/1024.0;
  mex::mbe mb(fg.factors());


  // Calculate elimination order(s) ////////////////////////////
  double timeOrder = 1; size_t nOrders = 1;
  if (vm.count("ordertime")) { timeOrder=vm["ordertime"].as<double>(); }
  if (vm.count("orders"))    { nOrders=vm["orders"].as<int>(); }
  if (vm.count("ibound"))    { ibound = vm["ibound"].as<int>(); }
    
  double startOrder = timeSystem();
  size_t iOrder = 1;

  if (ibound==0) order = fg.order( mex::graphModel::OrderMethod::Random );
  else {
    order = fg.order( mex::graphModel::OrderMethod::MinWidth );
    InducedWidth = fg.inducedWidth(order);

    // Repeat process until time or count limit reached
    while (iOrder < nOrders && (timeSystem()-startOrder < timeOrder)) {
      mex::VarOrder newOrder = fg.order(mex::graphModel::OrderMethod::MinWidth);
      //mex::VarOrder newOrder = fg.order(mex::graphModel::OrderMethod::MinFill);
      size_t newWidth = fg.inducedWidth(newOrder);
      if (newWidth < InducedWidth) { InducedWidth=newWidth; order=newOrder; }
      ++iOrder;
    }
    std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<"\n";
  }
  mb.setOrder(order);


  // simulate mini-bucket and get resulting maximal cliques
  mb.setProperties("ElimOp=MaxUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=1,DoHeur=0");
  mb.setIBound(ibound); 
  double mbMem = mb.simulateMemory(NULL,NULL, mbCutoff);
  while ( ibound > 0 && mbMem >= mbCutoff ) { 
    std::cout<<"MBE iBound "<<ibound<<" = "<<mbMem*mbMemRatio<<"M\n";
    mb.setIBound(--ibound); 
    mbMem=mb.simulateMemory(NULL,NULL, mbCutoff);
  }
  std::cout<<"MBE iBound "<<ibound<<" = "<<mbMem*mbMemRatio<<"M\n";
  if (1) {
    std::cout<<"Overriding with blank pseudo-tree\n";
    mex::vector<mex::mbe::vindex> pt(order.size()); for (size_t v=0;v<pt.size();++v) pt[v]=-1;
    mb.setPseudotree(pt);
  }
  mb.init();


  //mex::vector<size_t> bestConfig(fg.nvar());
  mex::vector<uint32_t> bestConfig(nvar, 0);
  //ofstream os(outfile); os << "MPE\n1\n";
  //os<<nvar<<" "; for (size_t v=0;v<nvar;++v) os<<"0 "; os<<"\n"; os.flush();

  double bestValue = -1e100;  //-mex::infty();
  double timeLeft, timeStop = timeSystem()+mplpTime;
  bool done = false;
  mb.reparameterize();
  while ( !done && (timeLeft = timeStop - timeSystem()) > 0 ) {
    mb.tighten( mplpIter , std::min(timeLeft,dt) , mplpObj , false);
    if ( timeLeft - (timeStop - timeSystem()) < dt ) done = true;  // if we quit for some other reason...

    double ub = 0.0, lb = 0.0;
    for (size_t f=0;f<mb.nFactors();++f) { ub+=mb.factor(f).max(); }

    // Try rounding in the MBE order
    lb=0.0;
    mex::vector<uint32_t> vals = mb.maxSequentialFast( order );
    for (size_t f=0;f<mb.nFactors();++f) { lb+=mb.factor(f)[ sub2ind(mb.factor(f).vars(),vals) ]; }
    if (lb > bestValue) { 
      bestValue = lb; 
      //bestConfig=vals; 
      for (size_t v=0;v<nvar;++v) if (!evVar.contains(Var(v,0))) bestConfig[v]=vals[v]; else bestConfig[v]=evid[v];
      ofstream os(outfile); os << "MPE\n1\n";
      os<<nvar<<" "; for (size_t v=0;v<nvar;++v) os<<bestConfig[v]<<" "; os<<"\n"; os.close();
      //os<<"-BEGIN-\n1\n"<<nvar<<" "; for (size_t v=0;v<nvar;++v) os<<bestConfig[v]<<" "; os<<"\n"; os.flush();
    }

    // Then try rounding in a random order
    lb=0.0;
    vals = mb.maxSequentialFast( fg.order(mex::graphModel::OrderMethod::Random) );
    for (size_t f=0;f<mb.nFactors();++f) { lb+=mb.factor(f)[ sub2ind(mb.factor(f).vars(),vals) ]; }
    if (lb > bestValue) { 
      bestValue = lb; 
      //bestConfig=vals; 
      for (size_t v=0;v<nvar;++v) if (!evVar.contains(Var(v,0))) bestConfig[v]=vals[v]; else bestConfig[v]=evid[v];
      ofstream os(outfile); os << "MPE\n1\n";
      os<<nvar<<" "; for (size_t v=0;v<nvar;++v) os<<bestConfig[v]<<" "; os<<"\n"; os.close();
      //os<<"-BEGIN-\n1\n"<<nvar<<" "; for (size_t v=0;v<nvar;++v) os<<bestConfig[v]<<" "; os<<"\n"; os.flush();
    }

    // Then try local hill climbing?
    for (size_t rnd=0;rnd<nLocal;++rnd) {
      for (size_t v=0; v<mb.nvar();++v) {
        if (mb.var(v).states()<2) continue;                      //   (make sure they're non-empty)
        if (evVar.contains(Var(v,0))) { bestConfig[v]=evid[v]; continue; }
        const factorGraph::flist& factors = mb.withVariable(mb.var(v)); // get the factors they're involved with
        Factor F(mb.var(v),0.0);
        for (factorGraph::flist::const_iterator f=factors.begin(); f!=factors.end(); ++f) {
          VarSet vs=mb.factor(*f).vars();  vs/=mb.var(v);           // and condition on their neighbors
          F += mb.factor(*f).slice(vs, sub2ind(vs,bestConfig));
        }
      }
      double score = 0.0; 
      for (size_t f=0;f<mb.nFactors();++f) {
        score += mb.factor(f)[sub2ind(mb.factor(f).vars(),bestConfig)];
      }
      if (score > bestValue) bestValue = score;
      //std::cout<<"Rescored; got "<<score<<"\n";
    }                                                        ///// end iterating local search


    std::cout<<"MPLP ("<<timeSystem()-timeStart<<") : ub="<<ub<<" ; lb="<<bestValue<<"\n";
    if (ub-lb < mplpGap) done = true;
  }

return 0;

}



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
double lbpTime, lbpIter, lbpObj, lbpErr;
double gbpTime, gbpIter, gbpObj;
double dt;

MEX_ENUM( Task , MPE,PR,MAR );

void writePR(const char* outfile, double logZ) {
  ofstream os(outfile);
  //os.precision(8); os.setf(ios::fixed,ios::floatfield);
  os<<"PR\n1\n"<<logZ/c_log10<<"\n";
  os.close();
}

void writeMAR(const char* outfile, mex::vector<Factor>& fs) {
  ofstream os(outfile);
  os<<"MAR\n1\n";
  os<<fs.size()<<" ";
  for (size_t f=0;f<fs.size();++f) {
    os<<fs[f].nrStates()<<" ";
    for (size_t i=0;i<fs[f].nrStates();++i) os<<fs[f][i]<<" ";
  }
  os<<"\n";
  os.close();
}


//Usage: solve_uai  -f <file.uai> -e <file.evid> -S <seed> -T <PR|MPE|MAR>
int main(int argc, char* argv[])
{

  double timeStart = timeSystem();

  const char* probName; // = argv[1];
  const char* taskName;//  = argv[3];
  Task task;
  mex::vector<Factor> bel;

  po::options_description desc("Available options");
  desc.add_options()
    ("help", "print help message")
    ("file,f", po::value<std::string>(), "input problem filename")
    ("evidence,e", po::value<std::string>(), "input evidence filename")
    ("seed,S", po::value<int>(),         "random number initial seed")
    ("task,T", po::value<std::string>(), "inference task string")
    ("ibound,i", po::value<int>(),       "initial i-bound")
    ("orders,o", po::value<int>(),       "number of variable orderings to try")
    ("ordertime,t", po::value<double>(), "max time spend on variable orderings")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("dt", po::value<double>(&dt)->default_value(300),         "file write time interval")
    ("lbps", po::value<double>(&lbpTime)->default_value(300),  "loopy belief propagation stop (seconds)")
    ("lbpi", po::value<double>(&lbpIter)->default_value(2000), "loopy belief propagation stop (iterations)")
    ("lbpe", po::value<double>(&lbpErr )->default_value(-1),   "loopy belief propagation stop (msg error)")
    ("lbpo", po::value<double>(&lbpObj )->default_value(-1),   "loopy belief propagation stop (objective)")
    ("gbps", po::value<double>(&gbpTime)->default_value(300),  "gen belief propagation stop (seconds)")
    ("gbpi", po::value<double>(&gbpIter)->default_value(-1),   "gen belief propagation stop (iterations)")
    ("gbpo", po::value<double>(&gbpObj )->default_value(-1),   "gen belief propagation stop (objective)")
    ("ijgp", "use ijgp regions only")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);

  if (vm.count("help")) { std::cout<<desc<<"\n"; return 1; }
  if (vm.count("file")) { probName=vm["file"].as<std::string>().c_str(); }
  else { std::cout<<"Missing input problem file!\n"; return 1; }
  if (vm.count("seed")) { mex::randSeed( vm["seed"].as<int>() ); }
  if (vm.count("task")) { taskName = vm["task"].as<std::string>().c_str(); task=Task(taskName); }
  else { std::cout<<"Missing task!\n"; return 1; }


  /*** READ IN PROBLEM FILE **********************************************************/

  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> flist = Factor::readUai10(is);
  size_t nvar=0;
  for (size_t f=0;f<flist.size();++f)                        // find maximum variable label
    nvar=std::max(nvar,(size_t)(flist[f].vars().rbegin()->label()+1));
  bel.resize(nvar);

  // Read in (single!) evidence 
  VarSet evVar;
  //ifstream is2( (std::string(probName)+".evid").c_str() );
  ifstream is2;
  if (vm.count("evidence")) { is2.open(vm["evidence"].as<std::string>().c_str()); }
  if (is2.is_open()) {
    std::cout<<"Got evidence file\n";
    std::map<uint32_t,size_t> evid;
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
          for (size_t v=0;v<overlap.nvar();++v) 
            bel[overlap[v].label()]=Factor::delta(overlap[v],evid[overlap[v].label()]);
          flist[f] = flist[f].condition( overlap, sub2ind(overlap,evid) );
        }
      }
    }
  } else std::cout<<"Evidence file not specified or not found\n";

  std::string outfiles(probName); outfiles += '.'; outfiles += taskName;
  std::string::size_type start = outfiles.find_last_of('/');
  if (start==std::string::npos) start=0; else ++start;
  outfiles = outfiles.substr(start,std::string::npos);
  const char* outfile = outfiles.c_str();
  std::cout<<"Writing to "<<outfile<<"\n";

  double ln10 = std::log(10);

  /*** PERFORM REQUESTED TASK ********************************************************/


  std::cout<<"Task is "<<task<<"\n";
  //std::cout<<"Task is "<<(const char*)task<<"\n";
  if (task==Task::MPE) { std::cout<<"MPE task not supported\n"; return 1; }

  /*** LOOPY BELIEF PROPAGATION ******************************************************/
  mex::lbp fg(flist); 
  std::cout<<"Model has "<<fg.nvar()<<" variables, "<<fg.nFactors()<<" factors\n";
if (lbpIter != 0 && lbpTime > 0) {
  fg.setProperties("Schedule=Priority,Distance=L1");
  fg.setStopIter(lbpIter); 
  fg.setStopMsg(lbpErr);
  fg.setStopObj(lbpObj);
  fg.setStopTime(lbpTime);
  fg.init();

  fg.run();
  switch (task) {
    case Task::PR: writePR(outfile,fg.logZ()); break;
    case Task::MAR: {
      for (size_t v=0;v<fg.nvar();++v) if (!evVar.contains(Var(v,0))) bel[v]=fg.belief( fg.localFactor(v) );
      writeMAR(outfile,bel);
    } break;
  }
  std::cout<<"LBP "<<fg.logZ()/ln10<<"\n";

  fg.reparameterize();                                         // Convert loopy bp results to model
}
  /*** HIGHER ORDERS ******************************************************/
  if ( gbpIter == 0 ) return 0;

  mex::gbp _gbp(fg.factors());                      // Create a GBP object for later use
  if (vm.count("ijgp")) _gbp.setMinimal(true);      // use "ijgp-like" regions? (true = remove all with c=0)
  else                  _gbp.setMinimal(false);

  /*** BUILD JUNCTION GRAPH & ASSESS MEMORY USE ***************************/
  bool exact = false;
  size_t ibound = 18, InducedWidth=10000;
  mex::VarOrder order;
  if (MemLimit > 0) {
    double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    mex::mbe mb(fg.factors());

    // Calculate elimination order(s) ////////////////////////////
    double timeOrder = 1; size_t nOrders = 1;
    if (vm.count("ordertime")) { timeOrder=vm["ordertime"].as<double>(); }
    if (vm.count("orders"))    { nOrders=vm["orders"].as<int>(); }
    
    double startOrder = timeSystem();
    size_t iOrder = 1;
    order = fg.order( mex::graphModel::OrderMethod::MinWidth );
    InducedWidth = fg.inducedWidth(order);

    // Repeat process until time or count limit reached
    while (iOrder < nOrders && (timeSystem()-startOrder < timeOrder)) {
      mex::VarOrder newOrder = fg.order(mex::graphModel::OrderMethod::MinFill);
      size_t newWidth = fg.inducedWidth(newOrder);
      if (newWidth < InducedWidth) { InducedWidth=newWidth; order=newOrder; }
      ++iOrder;
    }
    std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<"\n";
    mb.setOrder(order);

    // For PR tasks, see if an exact solution is within reach /////////////
    if (task == Task::PR) {
      mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0");
      mb.setIBound(InducedWidth); double mbMem = mb.simulateMemory(NULL,NULL,mbCutoff);
      if (mbMem < mbCutoff) {
        exact = true;
        mb.init();
        writePR(outfile,mb.logZ()); 
        std::cout<<"Exact solution by MBE: "<<mb.logZ()/ln10<<"\n";
        return 0;
      }
    }

    // simulate mini-bucket and get resulting maximal cliques
    mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=1");
    mex::vector<VarSet> cliques;
    double mem = std::numeric_limits<double>::infinity();
    if (vm.count("ibound")) { ibound = vm["ibound"].as<int>(); }
    mb.setIBound(ibound); 
    mbCutoff *= 2;  // leave a little slack for mismatch in criteria
    double mbMem = mb.simulateMemory(&cliques,NULL, mbCutoff);
    if (mbMem < mbCutoff) { _gbp.addRegions(cliques); mem = _gbp.memory(); }
    while ( ibound > 0 && (mbMem >= mbCutoff || mem > MemLimit) ) { 
      std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
      mb.setIBound(--ibound); cliques.clear(); mbMem=mb.simulateMemory(&cliques,NULL, mbCutoff);
      if (mbMem < mbCutoff) { _gbp.clearRegions(); _gbp.addRegions(cliques); mem=_gbp.memory(); }
    }
    std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
    mb = mex::mbe();                                     // clear mini-bucket (no longer needed)
  } else {
    // Didn't build MBE => no way to check ibound
    ibound = 0;
  }

  // Run GBP on the region graph 
  std::cout<<"GBP with "<<_gbp.nRegions()<<" regions; mem "<<_gbp.memory()<<"M\n";
  _gbp.setProperties("Schedule=Fixed,StopIter=200,StopObj=-1,StopMsg=1e-6");
  _gbp.init();
  _gbp.setStopIter(gbpIter); _gbp.setStopObj(gbpObj); _gbp.setStopMsg(-1.0); 

  // Get region indices for single-variable beliefs
  mex::vector<mex::gbp::findex> regions(fg.nvar());
  for (size_t v=0;v<fg.nvar();++v) {
    if (!evVar.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
    //std::cout<<"Region for "<<v<<" is "<<regions[v]<<"\n";
  }

  double gbpLeft, gbpStop = timeSystem()+gbpTime;
  while ( (gbpLeft = gbpStop - timeSystem()) > 0 ) {
    _gbp.setStopTime( std::min( dt , gbpLeft ) );
    _gbp.run();
    //if (_gbp.dObj() < 10)
    switch (task) {
      case Task::PR: writePR(outfile,_gbp.logZ()); break; 
      case Task::MAR: {
        //for (size_t v=0;v<_gbp.nvar();++v) if (!evVar.contains(Var(v,0))) bel[v]=_gbp.factor(regions[v]);
        for (size_t v=0;v<fg.nvar();++v) if (!evVar.contains(Var(v,0))) bel[v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
        writeMAR(outfile,bel);
      } break;
    }
    //else std::cout<<"Skipping output; presumed bad\n";
    std::cout<<"GBP "<<_gbp.logZ()/ln10<<"\n";
    if (_gbp.dObj() < gbpObj) { std::cout<<"Reached objective tolerance\n"; break; }
    if (_gbp.iter() >= gbpIter && gbpIter > 0) { std::cout<<"Reached iteration limit\n"; break; }
    if (_gbp.logZ() == -mex::infty()) { std::cout<<"Model deemed inconsistent\n"; break; }
  }

  if (InducedWidth <= ibound) {
    std::cout<<"Answer should be exact\n";
    return 0;
  }

std::cout<<"Quitting after GBP\n";
return 0;


  /*** ITERATIVE CONDITIONING AND GBP *************************************/
  VarSet cond;
  exact = (ibound >= InducedWidth);  
  
  while (!exact) {
    cond += fg.bestConditioner(order,cond);
    std::cout<<"Conditioning "<<cond<<"\n";

    double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    mex::mbe mb(fg.factors());
    if (order.size()==0) {
      order=fg.order(mex::graphModel::OrderMethod::MinFill);
      InducedWidth = fg.inducedWidth(order);
    }

    mex::vector<Factor> blank;
    _gbp = mex::gbp(blank);
    if (vm.count("ijgp")) _gbp.setMinimal(true); else _gbp.setMinimal(false);  // use "ijgp-like" regions?

    if (vm.count("ibound")) { ibound = vm["ibound"].as<int>(); } mb.setIBound(ibound);
    mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=1"); mb.setOrder(order);
    mex::vector<VarSet> cliques;
    double mem = std::numeric_limits<double>::infinity();
    //std::cout<<"Searching... "; std::cout.flush();
    double mbMem = mb.simulateMemory(&cliques, &cond, mbCutoff);
    if (mbMem < mbCutoff) { _gbp.addRegions(cliques); mem = _gbp.memory(); }
    while ( ibound > 0 && (mbMem >= mbCutoff || mem > MemLimit) ) {
      //std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
      mb.setIBound(--ibound); cliques.clear(); mbMem=mb.simulateMemory(&cliques, &cond, mbCutoff);
      if (mbMem < mbCutoff) { _gbp.clearRegions(); _gbp.addRegions(cliques); mem=_gbp.memory(); }
    }
    // !!! Test for exactness here: (1) set flag
    std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
    mb = mex::mbe();                                     // clear mini-bucket (no longer needed)

    if (task==Task::MAR) {
      for (size_t v=0;v<fg.nvar();++v) {
        if (!evVar.contains(Var(v,0))&& !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
      }
    }

    Factor lnZ(cond);
    bool failed = false;
    mex::vector<mex::vector<Factor> > condMarginals;
    for (size_t i=0;i<lnZ.nrStates();++i) {
      std::map<Var,size_t> val;  ind2sub(cond,i,val);
      mex::vector<Factor> fcond = fg.factors();
      for (size_t f=0;f<fcond.size();++f) {
        VarSet isect = cond & fcond[f].vars();
        if (isect.size() > 0) fcond[f] = fcond[f].condition(isect, sub2ind(isect,val));
      }

/*
  mex::lbp condfg(fcond);
  condfg.setProperties("Schedule=Priority,Distance=L1");
  condfg.setStopIter(lbpIter); condfg.setStopMsg(lbpErr); condfg.setStopObj(lbpObj); condfg.setStopTime(lbpTime);
  condfg.init();
  condfg.run();
  condfg.reparameterize();                                         // First convert loopy bp results to model
  _gbp.setFactors(condfg.factors()); _gbp.setProperties("Schedule=Fixed");
*/

      _gbp.setFactors(fcond); _gbp.setProperties("Schedule=Fixed");
      _gbp.init();
      _gbp.setStopIter(gbpIter); _gbp.setStopObj(gbpObj); _gbp.setStopMsg(-1.0); _gbp.setStopTime( 10000 ); _gbp.setVerbose(0);
      _gbp.run();
      if (_gbp.dObj() >= gbpObj) { failed=true; break; }     // convergence failure on this condition
      lnZ[i] = _gbp.logZ();
      if (task==Task::MAR) {
        condMarginals.push_back(bel);
        for (size_t v=0;v<fg.nvar();++v)
          if (!evVar.contains(Var(v,0)) && !cond.contains(Var(v,0))) 
            condMarginals[i][v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
      }
      for (size_t v=0;v<cond.size();++v) std::cout<<cond[v]<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
    }
    if (failed) {std::cout<<"Failing out\n"; continue; }           // if we failed out, condition on more variables
    double lnZtot = lnZ.logsumexp();
    switch (task) {
      case Task::PR:
        writePR(outfile, lnZtot);
        break;
      case Task::MAR:
        Factor probs = (lnZ - lnZtot).exp();
        for (size_t v=0;v<fg.nvar();++v) {
          if (evVar.contains(Var(v,0))) { } // evidence variables not updated
          else if (cond.contains(Var(v,0)))  { bel[v] = probs.marginal(Var(v,0)); }
          else { 
            bel[v] = condMarginals[0][v] * probs[0];
            for (size_t i=1;i<lnZ.nrStates();++i) bel[v] += condMarginals[i][v] * probs[i];
          }
        }
        writeMAR(outfile, bel);
        break;
    }
    std::cout<<"Conditioning "<<cond<<" => "<<lnZtot<<" ("<<lnZtot/ln10<<")\n";

  }
  // run CCS expansion on larger regions until timeout (?)

  return 0;

}


/*
void condition(it first, it last, varset vs, maptype val) {
  for (;first != last; ++first) { VarSet vf=first->vars()&vs; *first = first->condition(vf, sub2ind(vf,val)); }
}
*/

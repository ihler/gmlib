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
double gbpTime, gbpIter;
double dt;

MEX_ENUM( Task , MPE,PR,MAR );

void writePR(const char* outfile, double logZ) {
  ofstream os(outfile);
  os.precision(8); os.setf(ios::fixed,ios::floatfield);
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


int main(int argc, char* argv[])
{

  double timeStart = timeSystem();

  //if (argc < 4) { cout<<"Usage: "<<argv[0]<<" <file.uai> <seed> <PR|MPE|MAR>\n"; return 0; }
  const char* probName; // = argv[1];
  const char* taskName;//  = argv[3];
  Task task;
  mex::vector<Factor> bel;
  mex::vector<size_t> dims;

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
  bel.resize(nvar); dims.resize(nvar);
  for (size_t f=0;f<flist.size();++f) for (size_t v=0;v<flist[f].nvar();++v) dims[v]=flist[v].vars()[v].states(); 


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

  for (size_t v=0;v<nvar;++v) {
    if (dims[v]==1) { evVar += Var(v,dims[v]); bel[v]=Factor(Var(v,dims[v]),1.0); }
  }

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

    double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    mex::mbe mb(flist);  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMplp=0,DoMatch=0,DoFill=0");
    if (task==Task::PR)       mb.setProperties("DoJG=0");
    else if (task==Task::MAR) mb.setProperties("DoJG=1");

    // Calculate elimination order(s) ////////////////////////////
    double timeOrder = 20; size_t nOrders = 25;
    if (vm.count("ordertime")) { timeOrder=vm["ordertime"].as<double>(); }
    if (vm.count("orders"))    { nOrders=vm["orders"].as<int>(); }
    
    double startOrder = timeSystem();
    size_t iOrder = 1;
    mex::VarOrder order = fg.order( mex::graphModel::OrderMethod::MinWidth );
    size_t InducedWidth = fg.inducedWidth(order);

    // Repeat process until time or count limit reached
    while (iOrder < nOrders && (timeSystem()-startOrder < timeOrder)) {
      mex::VarOrder newOrder = fg.order(mex::graphModel::OrderMethod::MinFill);
      size_t newWidth = fg.inducedWidth(newOrder);
      if (newWidth < InducedWidth) { InducedWidth=newWidth; order=newOrder; }
      ++iOrder;
    }
    std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<"\n";
    mb.setOrder(order);
    mb.setIBound(InducedWidth);

    VarSet cond;
    mex::vector<VarSet> cliques;
    mex::vector<Factor> blank;
    mex::gbp _gbp(blank);
    double mem = -1;
    if (task==Task::MAR) mbCutoff *= 2;

    double  mbMem = mb.simulateMemory(&cliques,NULL,mbCutoff);
    if (task==Task::MAR) { _gbp.addRegions(cliques); mem = _gbp.memory(); }
    std::cout<<"Starting: "<<mbMem*sizeof(double)/1024/1024<<"MB\n";
    while (mbMem > mbCutoff || mem > MemLimit) {
      cond += fg.bestConditioner(order,cond);
      cliques.clear();
      mbMem = mb.simulateMemory(&cliques,&cond,mbCutoff);
      if (task==Task::MAR) { _gbp.clearRegions(); _gbp.addRegions(cliques); mem=_gbp.memory(); }
std::cout<<_gbp.nRegions()<<" => "<<cliques.size()<<" cliques\n";
      std::cout<<"Conditioning "<<cond<<" => "<<mbMem*sizeof(double)/1024/1024<<"MB, "<<mem<<"MB\n";
    }
    std::cout<<"Conditioner has "<<cond.nrStates()<<" values\n";

    Factor lnZ(cond);
    mex::vector<mex::vector<Factor> > condMarginals;
    mex::vector<mex::gbp::findex> regions(fg.nvar());
    if (task==Task::MAR) for (size_t v=0;v<fg.nvar();++v) {
      if (!evVar.contains(Var(v,0))&& !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
    }

    for (size_t i=0;i<lnZ.nrStates();++i) {
      std::map<Var,size_t> val;  ind2sub(cond,i,val);
      mex::vector<Factor> fcond = flist;
      for (size_t f=0;f<fcond.size();++f) {
        VarSet isect = cond & fcond[f].vars();
        if (isect.size() > 0) fcond[f] = fcond[f].condition(isect, sub2ind(isect,val));
      }
      if (task==Task::PR) {
        mex::mbe mbc(fcond); mbc.setOrder(order); mbc.setIBound(InducedWidth); mbc.setProperties("ElimOp=SumUpper");
        mbc.init();
        lnZ[i] = mbc.logZ();
      } else if (task == Task::MAR) {
        _gbp.setFactors(fcond); _gbp.setProperties("Schedule=Fixed");
        _gbp.init();
        _gbp.setStopIter(-1); _gbp.setStopObj(1e-7); _gbp.setStopMsg(-1.0); _gbp.setStopTime( 10000 ); _gbp.setVerbose(0);
        _gbp.run();
        lnZ[i] = _gbp.logZ();
        if (task==Task::MAR) {
          condMarginals.push_back(bel);
          for (size_t v=0;v<fg.nvar();++v)
            if (!evVar.contains(Var(v,0)) && !cond.contains(Var(v,0)))
              condMarginals[i][v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
        }
      }

      for (size_t v=0;v<cond.size();++v) std::cout<<cond[v]<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
    }
    std::cout<<"lnZ scores "<<lnZ<<"\n";
    double lnZtot = lnZ.logsumexp();
    std::cout<<"Final lnZ "<<lnZtot<<"\n";
    if (task==Task::PR) writePR(outfile,lnZtot);
    if (task==Task::MAR) {
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
    }

  return 0;

}




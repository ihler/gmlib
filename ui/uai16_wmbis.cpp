// Wrapper for UAI competition code
// 

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "boost/program_options.hpp"

#include "graphmodel.h"
#include "wmbe.h"


using namespace std;
using mex::mxObject;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::vector;
using mex::graphModel;

using mex::timeSystem;

namespace po = boost::program_options;

#define c_log10 std::log(10)

double MemLimit;
double dt;
int nOrders,nExtra;
int iBound;
double timeOrder;
double stopIter, stopTime;
double nSample;
double sampleTime;
double dampTheta, stepWeights;
double confidence;
int isAdaptive;
//double memUseRandom = mex::infty(); // if memory *way* out of reach, just use random ordering...
double memUseRandom = std::exp(40.0); // if memory *way* out of reach, just use random ordering...

const char* outfile;
int writetofile;

MEX_ENUM( Task , MPE,PR,MAR );

void writePR(const char* outfile, double logZ) {
  ofstream os(outfile);
  os<<"PR\n"<<logZ/c_log10<<"\n";
  os.close();
  std::cout<<"Wrote PR : "<<logZ/c_log10<<"\n";
}

void writeMAR(const char* outfile, mex::vector<Factor>& fs) {
  ofstream os(outfile);
  os<<"MAR\n";
  os<<fs.size()<<" ";
  for (size_t f=0;f<fs.size();++f) {
    os<<fs[f].nrStates()<<" ";
    double Z = fs[f].sum();  // ensure normalized marginals
    for (size_t i=0;i<fs[f].nrStates();++i) os<<fs[f][i]/Z<<" ";
  }
  os<<"\n";
  os.close();
  std::cout<<"Wrote MAR\n";
}




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
    ("orders,o",    po::value<int>(&nOrders)->default_value(1),      "number of variable orderings to try")
    ("order-time,t",po::value<double>(&timeOrder)->default_value(1), "max time spend on variable orderings")
    ("order-rand",  po::value<int>(&nExtra)->default_value(0),   "var order randomness; n=-1 (none), or among best+n")
    ("order-file", po::value<std::string>(), "problem elimination ordering filename")
    ("order-limit", po::value<double>(&memUseRandom)->default_value(mex::infty()), "fail-out memory limit for ordering search")
    ("ibound,i", po::value<int>(&iBound)->default_value(30),       "initial i-bound")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("stopiter,n",  po::value<double>(&stopIter)->default_value(1.0),      "maximum number of iterations")
    ("stoptime",    po::value<double>(&stopTime)->default_value(20.0),     "maximum iteration time (sec)")
    ("damp",      po::value<double>(&dampTheta)->default_value(1.0),   "damping of factor reparameterization, real[0,1]: 0.0=no update, 1.0=no damping")
    ("stepw",     po::value<double>(&stepWeights)->default_value(0.1),   "step size of weight update; 0.0=no update")
    ("adaptive",  po::value<int>(&isAdaptive)->default_value(1),   "Adaptive step sizes; 0=no, 1=yes (default)")
    ("nsample,s",  po::value<double>(&nSample)->default_value(10000),      "max number of samples to draw")
    ("sampletime", po::value<double>(&sampleTime)->default_value(60.0),    "maximum time for sampling (sec)")
    ("confidence,c",  po::value<double>(&confidence)->default_value(0.05),  "confidence (probability of bound exceedence)")
    ("write",  po::value<int>(&writetofile)->default_value(0),  "write results to file (UAI); default 0")
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

  if ((task != Task::PR)&&(task != Task::MAR)) {std::cout<<"Task should be PR or MAR\n"; return 1; }


  /*** READ IN PROBLEM FILE **********************************************************/
  std::cout<<"Reading model file: "<<probName<<"\n";
  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> forig = Factor::readUai10(is);
  mex::wmbe wmb( forig );
  size_t nvar = wmb.nvar();
  bel.resize(nvar); 
  mex::vector<size_t> evid(nvar, 0);
std::cout<<"Graphical model over "<<wmb.nvar()<<" variables with "<<forig.size()<<" factors\n";


  /*** READ IN EVIDENCE FILE *********************************************************/
  VarSet evVar;
  ifstream is2;
  if (vm.count("evidence")) { is2.open(vm["evidence"].as<std::string>().c_str()); }
  if (is2.is_open()) {
    std::cout<<"Got evidence file\n";
    //int nEvid; is2 >> nEvid;
    int nEvid=1;      // 2014 format: single evidence, no # of evidences entry
    std::cout<<"Got "<<nEvid<<" evidences?\n";
    if (nEvid > 0) {
      int nEvidVar; is2 >> nEvidVar;
      for (size_t i=0;i<nEvidVar;i++) {
        uint32_t vid; size_t vval; is2>>vid>>vval; 
        evid[vid]=vval; evVar |= wmb.var(vid);
        bel[vid] = Factor::delta(wmb.var(vid),vval);
      }
      std::cout<<"Evidence on variables "<<evVar<<"\n";
      wmb.condition(evVar,evid);
    }
  } else std::cout<<"Evidence file not specified or not found\n";



  /*** PREPARE REQUESTED TASK ********************************************************/
  std::cout<<"Task is "<<task<<"\n";
  //std::cout<<"Task is "<<(const char*)task<<"\n";
  //if (task==Task::MPE) { std::cout<<"MPE task not supported\n"; return 1; }

  /*** PREPARE OUTPUT FILE ***********************************************************/
  std::string outfiles(probName); outfiles += '.'; outfiles += taskName;
  std::string::size_type start = outfiles.find_last_of('/');
  if (start==std::string::npos) start=0; else ++start;
  outfiles = outfiles.substr(start,std::string::npos);
  outfile = outfiles.c_str();
  //const char* outfile = outfiles.c_str();
  std::cout<<"Writing to "<<outfile<<"\n";
  std::cout<<"Using log10 (UAI)\n";


  /*** ELIMINATION ORDERS ************************************************************/

    double startOrder = timeSystem();
    mex::VarOrder order;
    size_t InducedWidth, iOrder = 1;
    double score = wmb.order( mex::graphModel::OrderMethod::Random, order, 0, memUseRandom );

    const char *orderFile = NULL;       // Check for pre-specified elimination order file
    if (vm.count("order-file")) { orderFile = vm["order-file"].as<std::string>().c_str(); }
    ifstream orderIStream; if (orderFile!=NULL) orderIStream.open(orderFile);
    if (orderIStream.is_open()) {
      // If we were given an input file with an elimination ordering, just use that
      std::cout << "Reading elimination order from "<<orderFile<<"\n";
      size_t ordersize;  orderIStream>>ordersize; assert(ordersize == order.size());
      for (size_t i=0;i<order.size();++i) { size_t tmp;  orderIStream>>tmp; order[i]=tmp;  };
      orderIStream.close();
      
    } else {
      // Otherwise, calculate elimination order(s) ////////////////////////////////////
      double startOrder = timeSystem();
      iOrder = 0;
      // Try to build new orders until time or count limit reached ////////////////////
      while (iOrder < nOrders && (timeSystem()-startOrder < timeOrder)) {
        score = wmb.order(mex::graphModel::OrderMethod::WtMinFill, order, nExtra, score);
        ++iOrder;
      }

      // If we were given an ordering file name but no file, write our order out 
      ofstream orderOStream; if (orderFile!=NULL) orderOStream.open(orderFile);
      if (orderOStream.is_open()) {
        std::cout << "Writing elimination order to "<<orderFile<<"\n";
        orderOStream<<order.size();
        for (size_t i=0;i<order.size();++i) orderOStream<<" "<<order[i];
        orderOStream<<"\n";
        orderOStream.close();
      }
    }

    InducedWidth = wmb.inducedWidth(order);
    std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<", score "<<score<<"\n";
    wmb.setOrder(order);

  /*** SET IBOUND ********************************************************************/
    wmb.setIBound(iBound);

    for (;;--iBound) {
      wmb.init();
      wmb.setIBound(iBound);
      wmb.build();
      double mem = wmb.memory();
      if (mem <= MemLimit || iBound == 1) break;
      std::cout<<"Projected memory too high ("<<iBound<<"=> "<<mem<<"MB)\n";
    }
  
  wmb.setTheta();
  if (task==Task::MPE) {
    for (size_t b=0;b<wmb.getOrder().size();++b) wmb.setWeightPositiveUniform(b, 1e-6);
    if (stepWeights > 0.0) std::cout<<"MPE task => stepWeights = 0.0\n";
    stepWeights = 0.0;
  } else if ((task==Task::PR)||(task==Task::MAR)) { 
    for (size_t b=0;b<wmb.getOrder().size();++b) wmb.setWeightPositiveUniform(b, 1.0);
  } else { std::cout<<"Task not supported\n"; return 1; }
    
  
  //wmb.dump(std::cout);  

  /*** RUN ITERATIONS ****************************************************************/

  double startIter = timeSystem();
  double best = mex::infty(), last = mex::infty();
  mex::vector<uint32_t> xhat( wmb.nvar() );
  double f_xbest = -mex::infty();
  double iter_time = 0.0;

  for (size_t it=0; it < stopIter && (timeSystem()-startIter < stopTime); ++it) {
    double obj = 0.0;
    for (size_t b=0;b<wmb.getOrder().size();++b) obj += wmb.msgForward(b, dampTheta, stepWeights);

    if (isAdaptive && obj > last) { dampTheta /= 2; stepWeights /= 2; }
    last = obj;
    std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<last/c_log10<<" , "<<obj/c_log10<<"\n";
    if (it == 0) iter_time = timeSystem()-startIter;  // Estimate iteration time for stopping
    if (it == stopIter-1) break;
    if (timeSystem() - startIter + 2*iter_time > stopTime) break;

    for (int b=wmb.getOrder().size()-1;b>=0;--b) wmb.msgBackward(b, dampTheta, stepWeights);
  }

  std::cout<<" ==== Beginning WMB importance sampling ====\n";

  double startSample = timeSystem();
  double Ex = 0.0, Ex2 = 0.0;
  double ub, pred, lb;
  double nextOutput = 8;
  std::pair<double,double> logPQ = wmb.sampleMixture( xhat );
  Ex  = std::exp( logPQ.first - logPQ.second - last);
  if (task == Task::MAR) { // If MAR, initialize beliefs
    for (size_t v=0;v<nvar;++v) bel[v] = Ex*Factor::delta(wmb.var(v),xhat[v]) + 1e-300;
  }
  Ex2 = Ex*Ex;
  for (size_t samp=1; samp < nSample && (timeSystem()-startSample < sampleTime); ++samp) {
    std::pair<double,double> logPQ = wmb.sampleMixture( xhat );
    double dEx = std::exp(logPQ.first - logPQ.second - last);
    //std::cout<<"  "<<logPQ.first<<"/"<<logPQ.second<<" => "<<Ex<<"\n";
    Ex  = (Ex * (samp-1))/samp + dEx/samp;
    Ex2 = (Ex2 * (samp-1))/samp + dEx*dEx/samp;
    if (task == Task::MAR) {  // If MAR, update beliefs given new sample
      for (size_t v=0;v<nvar;++v) { bel[v]*=((double)samp-1)/samp; bel[v][sub2ind(bel[v].vars(),xhat)] += dEx/samp; }
    }
    double var = std::max(Ex2 - Ex*Ex, 0.0);
    //double rng = std::sqrt(2*var*std::log(2.0/confidence)/samp) + 7*std::log(2.0/confidence)/3.0/(samp-1);
    // var *= samp / (samp-1) to make unbiased, but then (var/samp) in range calculation => just divide by samp-1
    double rng = std::sqrt(2*var*std::log(2.0/confidence)/(samp-1)) + 7*std::log(2.0/confidence)/3.0/(samp-1);
    //std::cout<<"  "<<Ex2<<"-"<<Ex*Ex<<","<<rng<<","<<last<<"\n";
    ub = std::log(Ex + rng) + last;
    pred = std::log(Ex) + last;
    lb = std::log(std::max(Ex - rng, 0.0)) + last;
    if (samp >= nextOutput) {
      std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<ub/c_log10<<" > "<<std::setw(10)<<pred/c_log10<<" > "<<std::setw(10)<<lb/c_log10<<"\n";
      nextOutput *= 1.25;
      if (writetofile) {  // Output to file if requested
        if (task == Task::PR) if (pred > -mex::infty()) writePR(outfile, pred); else writePR(outfile, ub);
        if (task == Task::MAR) writeMAR(outfile, bel);
      }
    }
  }
  std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<ub/c_log10<<" > "<<std::setw(10)<<pred/c_log10<<" > "<<std::setw(10)<<lb/c_log10<<"\n";

  return 0;

}




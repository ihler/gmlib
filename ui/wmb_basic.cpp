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
double dampTheta, stepWeights;
double addEpsilon;
int isAdaptive;
int startWeights;
double memUseRandom = mex::infty(); // if memory *way* out of reach, just use random ordering...


MEX_ENUM( Task , MPE,PR,MAR,LB );


int main(int argc, char* argv[])
{

  double timeStart = timeSystem();

  //if (argc < 4) { cout<<"Usage: "<<argv[0]<<" <file.uai> <seed> <PR|MPE|MAR>\n"; return 0; }
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
    ("ibound,i", po::value<int>(&iBound)->default_value(30),       "initial i-bound")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("stopiter,n",  po::value<double>(&stopIter)->default_value(1.0),      "maximum number of iterations")
    ("stoptime",    po::value<double>(&stopTime)->default_value(20.0),     "maximum iteration time (sec)")
    ("damp",      po::value<double>(&dampTheta)->default_value(1.0),   "damping of factor reparameterization, real[0,1]: 0.0=no update, 1.0=no damping")
    ("startw",    po::value<int>(&startWeights)->default_value(1),   "Weight init: 0=basic MBE; 1=uniform weights (default)")
    ("stepw",     po::value<double>(&stepWeights)->default_value(0.1),   "step size of weight update; 0.0=no update")
    ("adaptive",  po::value<int>(&isAdaptive)->default_value(1),   "Adaptive step sizes; 0=no, 1=yes (default)")
    ("eps",     po::value<double>(&addEpsilon)->default_value(0.0),   "constant to add to all factors (to avoid hard zeros)")
    ("writefile",  po::value<std::string>(),   "Write out a reparameterized file")
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
	std::cout<<"Reading model file: "<<probName<<"\n";
  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> forig = Factor::readUai10(is);
	if (addEpsilon > 0.0) for (size_t f=0;f<forig.size();++f) forig[f] += addEpsilon;
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

	switch (task) {
		case Task::MPE:
			wmb.setElimType(mex::wmbe::ElimType::MaxUpper);
			if (stepWeights > 0.0) std::cout<<"MPE task => stepWeights = 0.0\n";
			stepWeights = 0.0;
			break;
		case Task::PR:
		case Task::MAR:
			wmb.setElimType(mex::wmbe::ElimType::SumUpper);
			if (startWeights==0) for (size_t b=0;b<nvar;++b) wmb.setWeightPositiveFirst(b,1.0);
			break;
		case Task::LB:
			wmb.setElimType(mex::wmbe::ElimType::SumLower);
			if (startWeights==0) for (size_t b=0;b<nvar;++b) wmb.setWeightNegativeFirst(b,1.0);
			break;
		otherwise:
  		std::cout<<"Task not supported\n"; return 1; 
	}

	
  //wmb.dump(std::cout);	

  /*** RUN ITERATIONS ****************************************************************/

  double startIter = timeSystem();
	double best = mex::infty(), last = mex::infty();
	mex::vector<uint32_t> xbest( wmb.nvar() ), xhat( wmb.nvar() );
	double f_xbest = -mex::infty();

	if (task == Task::LB) { best = -best; last = -last; }			// set defaults to correct extremum for lower bounds

  for (size_t it=0; it < stopIter && (timeSystem()-startIter < stopTime); ++it) {
    double obj = 0.0;
    for (size_t b=0;b<wmb.getOrder().size();++b) obj += wmb.msgForward(b, dampTheta, stepWeights);

		if (task == Task::LB) {
			if (isAdaptive && obj < last) { dampTheta /= 2; stepWeights /= 2; }
			last = obj;
			best = std::max(best,obj);
		} else {
			if (isAdaptive && obj > last) { dampTheta /= 2; stepWeights /= 2; }
			last = obj;
			best = std::min(best,obj);
		}

		if (task == Task::MPE) {
			if (1) { 
				double f_xhat = wmb.maxSequential(xhat);
				if (f_xhat > f_xbest) { xbest = xhat; f_xbest = f_xhat; }	// can also use wmb.logP(xhat) to evaluate
			}
			if (1) { 
				double f_xhat = wmb.sampleSequential(xhat).first;
				if (f_xhat > f_xbest) { xbest = xhat; f_xbest = f_xhat; }
			}
			std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<best<<" , "<<obj<<" , "<<f_xbest<<"\n";
		} else {
			std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<best<<" , "<<obj<<"\n";
		}

		if (it == stopIter) break;
		if (timeSystem()-startIter > stopTime) break;

    for (int b=wmb.getOrder().size()-1;b>=0;--b) wmb.msgBackward(b, dampTheta, stepWeights);
  }

  if (vm.count("writefile")) { 
		wmb.reparameterize();
  	ofstream ofs;
  	ofs.open(vm["writefile"].as<std::string>().c_str()); 
  	if (ofs.is_open()) {
			std::cout<<"Writing reparameterized model as "<<vm["writefile"].as<std::string>()<<"\n";
			Factor::writeUai10( ofs , wmb.factors() );
		}
	}



  return 0;

}




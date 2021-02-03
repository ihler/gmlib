// Wrapper for UAI competition code
// 

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
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
double nSearch;
double searchTime;
double dampTheta, stepWeights;
double addEpsilon;
int isAdaptive;
double memUseRandom = mex::infty(); // if memory *way* out of reach, just use random ordering...


MEX_ENUM( Task , MPE,PR,MAR );


struct searchNode {
	double priority;
	size_t depth;
	mex::vector<uint32_t> tuple;
	double ub, lb, costUB, costLB;
	bool operator< (const searchNode& n) const { return priority < n.priority; }
};
double calcPriority( searchNode n ) {
	//return n.depth;
	//return (n.ub-n.lb)*std::log((double)n.depth);
	return n.ub + std::log( 1 - std::exp( std::min(n.ub,n.lb)-n.ub ) );
}


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
    ("stepw",     po::value<double>(&stepWeights)->default_value(0.1),   "step size of weight update; 0.0=no update")
    ("adaptive",  po::value<int>(&isAdaptive)->default_value(1),   "Adaptive step sizes; 0=no, 1=yes (default)")
    ("eps",     po::value<double>(&addEpsilon)->default_value(0.0),   "constant to add to all factors (to avoid hard zeros)")
    ("nsearch,s",  po::value<double>(&nSearch)->default_value(10000),      "max number of search nodes")
    ("searchtime", po::value<double>(&searchTime)->default_value(60.0),    "maximum time for searching (sec)")
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

	if (task != Task::PR) {std::cout<<"Task should be PR\n"; return 1; }

  /*** READ IN PROBLEM FILE **********************************************************/
	std::cout<<"Reading model file: "<<probName<<"\n";
  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> forig = Factor::readUai10(is);
  if (addEpsilon > 0.0) for (size_t f=0;f<forig.size();++f) forig[f] += addEpsilon;
	mex::wmbe wmbUB( forig ), wmbLB( forig );
  size_t nvar = wmbUB.nvar();
  bel.resize(nvar); 
	mex::vector<size_t> evid(nvar, 0);
std::cout<<"Graphical model over "<<wmbUB.nvar()<<" variables with "<<forig.size()<<" factors\n";


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
        evid[vid]=vval; evVar |= wmbUB.var(vid);
				bel[vid] = Factor::delta(wmbUB.var(vid),vval);
      }
			std::cout<<"Evidence on variables "<<evVar<<"\n";
			wmbUB.condition(evVar,evid);
			wmbLB.condition(evVar,evid);
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
    double score = wmbUB.order( mex::graphModel::OrderMethod::Random, order, 0, memUseRandom );

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
        score = wmbUB.order(mex::graphModel::OrderMethod::WtMinFill, order, nExtra, score);
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

    InducedWidth = wmbUB.inducedWidth(order);
    std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<", score "<<score<<"\n";
    wmbUB.setOrder(order);
    wmbLB.setOrder(order);

  /*** SET IBOUND ********************************************************************/
    wmbUB.setIBound(iBound);
    wmbLB.setIBound(iBound);

		for (;;--iBound) {
			wmbUB.init();
			wmbLB.init();
			wmbUB.setIBound(iBound);
			wmbLB.setIBound(iBound);
			wmbUB.build();
			wmbLB.build();
			double mem = wmbUB.memory() + wmbLB.memory();
      if (mem <= MemLimit || iBound == 1) break;
      std::cout<<"Projected memory too high ("<<iBound<<"=> "<<mem<<"MB)\n";
		}
  
  wmbUB.setTheta();
  wmbLB.setTheta();

  switch (task) {
    case Task::PR:
      wmbUB.setElimType(mex::wmbe::ElimType::SumUpper);
      wmbLB.setElimType(mex::wmbe::ElimType::SumLower);
      break;
    otherwise:
      std::cout<<"Task not supported\n"; return 1;
  }

  //wmb.dump(std::cout);	

  /*** RUN ITERATIONS ****************************************************************/

  double startIter = timeSystem();
  double UB = 0.0, LB = 0.0;

  for (size_t it=0; it < stopIter && (timeSystem()-startIter < stopTime); ++it) {
		UB = 0.0; LB = 0.0;
    for (size_t b=0;b<wmbUB.getOrder().size();++b) UB += wmbUB.msgForward(b, dampTheta, stepWeights);
    for (size_t b=0;b<wmbLB.getOrder().size();++b) LB += wmbLB.msgForward(b, dampTheta, stepWeights);

		//if (isAdaptive && obj > last) { dampTheta /= 2; stepWeights /= 2; }
		//last = obj;
		std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<UB<<" > ln Z > "<<std::setw(10)<<LB<<"\n";
		if (it == stopIter) break;
		if (timeSystem()-startIter > stopTime) break;

		for (int b=wmbUB.getOrder().size()-1;b>=0;--b) wmbUB.msgBackward(b, dampTheta, stepWeights);
		for (int b=wmbLB.getOrder().size()-1;b>=0;--b) wmbLB.msgBackward(b, dampTheta, stepWeights);
  }

	//wmbUB.dump(std::cout);

	std::cout<<" ==== Beginning WMB search process ====\n";

  double startSearch = timeSystem();
	typedef std::priority_queue<searchNode> searchQueue;
	searchQueue queue;
	searchNode root; root.depth = 0; root.tuple.clear(); root.ub = UB; root.lb = LB; root.priority = calcPriority(root); root.costUB = root.costLB = 0.0;
	queue.push( root );
	//queue.insert( std::pair<double,searchNode>( calcPriority(root) , root) );

	mex::vector<uint32_t> tuple(nvar);

	for(size_t nodes=0; nodes < nSearch; ++nodes) {
		if (queue.empty()) break;										// out of nodes => done
		searchNode top = queue.top();	queue.pop();	// get top node & remove from queue

		size_t d = top.depth; 
		Var X = wmbUB.var( order[nvar - d - 1] );		// get next variable in reverse order
		for (size_t v=0;v<d;++v) tuple[ order[nvar-v-1] ] = top.tuple[v]; 	// copy out config for usage
		//std::cout<<"("; for (size_t v=0;v<d;++v) std::cout<<order[nvar-v-1]; std::cout<<") = ["; for (size_t v=0;v<d;++v) std::cout<<top.tuple[v]; std::cout<<"]\n";

		top.depth += 1;
		top.tuple.push_back(0);

		Factor dUB(X,0.0), dLB(X,0.0);
		// TODO: COMPUTATION INCORRECT	
		for (size_t v=0; v<X.states(); ++v) {		// TODO: just compute f(x) directly?
			searchNode nAdd_v = top;
			nAdd_v.tuple.back() = v; 											// Create new node from top + X=v
			tuple[ X ] = v;
			//std::cout<<" Check "<<v<<": "<<top.costUB<<"+"<<wmbUB.heuristicTheta(X,tuple)<<"+"<<wmbUB.heuristicIn(X,tuple)<<"\n";
			nAdd_v.costUB += wmbUB.heuristicTheta(X,tuple);
			nAdd_v.costLB += wmbLB.heuristicTheta(X,tuple);
			nAdd_v.ub = dUB[v] = nAdd_v.costUB + wmbUB.heuristicIn(X,tuple);		// evaluate upper & lower heuristics
			nAdd_v.lb = dLB[v] = nAdd_v.costLB + wmbLB.heuristicIn(X,tuple);
			nAdd_v.priority = calcPriority(nAdd_v);  // nodes;
			if (nAdd_v.depth < nvar) queue.push( nAdd_v );
		}
		// Adjust overall upper & lower bounds
		//std::cout<<"Adjust "<<top.ub<<" : "<<dUB<<"=>"<<dUB.logsumexp()<<"\n   LB: "<<top.lb<<" : "<<dLB<<"=>"<<dLB.logsumexp()<<"\n";
		// ub = log(  exp(ub) - exp(hUB) + exp(dUB) )
		UB += std::log( 1 - std::exp(top.ub - UB) + std::exp( dUB.logsumexp() - UB ) );
		LB += std::log( 1 - std::exp(top.lb - LB) + std::exp( dLB.logsumexp() - LB ) );

		if ((nodes%10000) == 0) {
			std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<UB<<" > ln Z > "<<std::setw(10)<<LB<<"  ("<<queue.size()<<")"<<"\n";\
			if (timeSystem()-startSearch > searchTime) break;
		}

		if (0) break;
	}
	std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<UB<<" > ln Z > "<<std::setw(10)<<LB<<"\n";\

  return 0;

}




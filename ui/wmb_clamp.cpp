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
double nSample;
double sampleTime;
double dampTheta, stepWeights;
double addEpsilon;
int isAdaptive;
double memUseRandom = mex::infty(); // if memory *way* out of reach, just use random ordering...


MEX_ENUM( Task , MPE,PR,MAR );


struct searchNode {
	double priority;
	size_t depth;
	mex::vector<uint32_t> tuple;
	double hUB, hLB;
	bool operator< (const searchNode& n) const { return priority < n.priority; }
};
double calcPriority( searchNode n ) {
	return n.hUB + std::log( 1 - std::exp( std::min(n.hUB,n.hLB)-n.hUB ) );
}

/*
struct searchNode {
  size_t v_parent;    // this node = (parent.X = v)
  Var X;      // variable expanded for children
  double f; // f = g + h ?
  searchNode* parent;
  vector<searchNode*> children;

  //void searchNode( const searchNode& parent, Var x,  
  void expand(Var x) {
    X = x;
    children.resize(X.states());
    for (size_t v=0; v<X.states(); ++v) {
      children[v].parent = this;
      children[v].f = -mex::infty(); // !!!
      children[v].v = v;
    }
  }
  void propagate( cold, cnew ) {
    double fold = f;
    //f = std::log( std::exp(f) - std::exp(cold) + std::exp(cnew) );
    f = std::log( 1 - std::exp(cold - f) + std::exp(cnew - f) );
    if (parent != NULL) parent->propagate( fold, f );
  }
}




*/




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
    ("nsearch,s",  po::value<double>(&nSample)->default_value(10000),      "max number of search nodes")
    ("searchtime", po::value<double>(&sampleTime)->default_value(60.0),    "maximum time for searching (sec)")
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


	typedef std::priority_queue<searchNode*> searchQueue;
	searchQueue queue;
	searchNode root; 
	mex::vector<uint32_t> tuple(nvar);
	queue.push( &root );

	// run wmb on root, get f , beliefs & choose X

	while (1) {
		if (queue.empty()) break;				// out of nodes => done
		searchNode* top = queue.top();	
		queue.pop();										// remove it before adding more nodes...
		top->children.resize(top->X.states());
		VarSet vs(top->X); 	// Build WMB model conditioned on parent's config 
		searchNode* n=top; 
		while (n!=&root) { tuple[n->parent->X]=n->v_parent; vs+=n->parent->X; n=n->parent;}
		for (size_t v=0;v<top->X.states();++v) {
			wmbe wmbv( forig ); tuple[top->X]=v; wmbv.condition(vs,tuple);
			// condition on X=v & run WMB
			top->children[v].f = wmb.obj;
			top->children[v].v_parent = v;
			// Find belief to condition on next as top->children[v].X 
		}		
		// fnew = logsumexp( children.f ); propagate to parents
		wmbe wmbi( forig ); wmbi.condition(evVar,evid);
		searchNode* n = top; while (n != &root) { wmbi.condition(n.parent->X,v); n = n->parent; } // !!! do differently
		// run wmb 5 iter?
		top->f = wmbi.best();
		// check each variable & select one as top->X
		// 
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
  double ub = 0.0, lb = 0.0;

  for (size_t it=0; it < stopIter && (timeSystem()-startIter < stopTime); ++it) {
		ub = 0.0; lb = 0.0;
    for (size_t b=0;b<wmbUB.getOrder().size();++b) ub += wmbUB.msgForward(b, dampTheta, stepWeights);
    for (size_t b=0;b<wmbLB.getOrder().size();++b) lb += wmbLB.msgForward(b, dampTheta, stepWeights);

		//if (isAdaptive && obj > last) { dampTheta /= 2; stepWeights /= 2; }
		//last = obj;
		std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<ub<<" > ln Z > "<<std::setw(10)<<lb<<"\n";
		if (it == stopIter) break;
		if (timeSystem()-startIter > stopTime) break;

		for (int b=wmbUB.getOrder().size()-1;b>=0;--b) wmbUB.msgBackward(b, dampTheta, stepWeights);
		for (int b=wmbLB.getOrder().size()-1;b>=0;--b) wmbLB.msgBackward(b, dampTheta, stepWeights);
  }

	std::cout<<" ==== Beginning WMB search process ====\n";

	//queue.insert( std::pair<double,searchNode>( calcPriority(root) , root) );


	for(size_t nodes=0; nodes < 100; ++nodes) {

		size_t d = top.depth; 
		Var X = wmbUB.var( order[nvar - d - 1] );			// get next variable in reverse order
		for (size_t v=0;v<d;++v) tuple[ order[nvar - v - 1] ] = top.tuple[v]; //top->second.tuple[v];	// copy out config for usage
		std::cout<<"("; for (size_t v=0;v<d;++v) std::cout<<order[nvar-v-1]; std::cout<<") = ["; for (size_t v=0;v<d;++v) std::cout<<top.tuple[v]; std::cout<<"]\n";

		top.depth += 1;
		top.tuple.push_back(0);
		Factor dUB(X,0.0), dLB(X,0.0);
		double hUB = wmbUB.heuristicPre( X, tuple ),  hLB = wmbLB.heuristicPre( X, tuple );
	
		// TODO: COMPUTATION INCORRECT	
		for (size_t v=0; v<X.states(); ++v) {		// TODO: just compute f(x) directly?
			searchNode nAdd_v = top;
			nAdd_v.tuple.back() = v; 											// Create new node from top + X=v
			tuple[ X ] = v;
			dUB[v] = wmbUB.heuristicPost( X, tuple );		// Evaluate upper & lower heuristics
			dLB[v] = wmbLB.heuristicPost( X, tuple );
			nAdd_v.hUB += dUB[v];
			nAdd_v.hLB += dLB[v];													// And add node to queue
			nAdd_v.priority = nodes; //calcPriority(nAdd_v);
			if (nAdd_v.depth < nvar) queue.push( nAdd_v );
		}
		// Adjust overall upper & lower bounds
		std::cout<<"Adjust "<<hUB<<" : "<<dUB<<"=>"<<dUB.logsumexp()<<"\n   LB: "<<hLB<<" : "<<dLB<<"=>"<<dLB.logsumexp()<<"\n";
		// ub = log(  exp(ub) - exp(hUB) + exp(dUB) )
		ub += std::log( 1 - std::exp(hUB - ub) + std::exp( dUB.logsumexp() - ub ) );
		lb += std::log( 1 - std::exp(hLB - lb) + std::exp( dLB.logsumexp() - lb ) );

		if (1) std::cout<<"["<<std::setw(10)<<timeSystem()-startIter<<"] : "<<std::setw(10)<<ub<<" > ln Z > "<<std::setw(10)<<lb<<"\n";

		if (0) break;
	}

  return 0;

}




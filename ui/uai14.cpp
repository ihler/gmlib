// Wrapper for UAI competition code
// 

#include <cstdio>
#include <iostream>
//#include <sstream>
#include <fstream>
#include "boost/program_options.hpp"

#include "factorgraph.h"

#include "mplp.h"
#include "mbe.h"
#include "gibbs.h"

#include "lbp.h"
#include "mixbp.h"
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

VarSet maxVars;

double MemLimit;
double lbpTime, lbpIter, lbpObj, lbpErr;
double gbpTime, gbpIter, gbpObj;
double dt;
int    doCond;
bool   doVerbose;
int    nExtra;
int    iboundInit;

double timeOrder;
int    nOrders;
double memUseRandom = mex::infty();	// if memory *way* out of reach, just use random ordering...

const char* outfile;

MEX_ENUM( Task , MPE,PR,MAR,MMAP );

void writePR(const char* outfile, double logZ) {
  ofstream os(outfile);
  //os.precision(8); os.setf(ios::fixed,ios::floatfield);
  //os<<"PR\n1\n"<<logZ/c_log10<<"\n";		// !!! 2011 PIC version: need "1" evidence
  os<<"PR\n"<<logZ/c_log10<<"\n";
  os.close();
	std::cout<<"Wrote PR : "<<logZ/c_log10<<"\n";
}

void writeMAR(const char* outfile, mex::vector<Factor>& fs) {
  ofstream os(outfile);
  //os<<"MAR\n1\n";		// !!! 2011 PIC version: need "1" evidence
  os<<"MAR\n";
  os<<fs.size()<<" ";
  for (size_t f=0;f<fs.size();++f) {
    os<<fs[f].nrStates()<<" ";
    for (size_t i=0;i<fs[f].nrStates();++i) os<<fs[f][i]<<" ";
  }
  os<<"\n";
  os.close();
	std::cout<<"Wrote MAR\n";
}

void writeMPE(const char* outfile, const mex::vector<mex::index>& xhat) {
  ofstream os(outfile);
  os<<"MPE\n";
  os<<xhat.size()<<" ";
  for (size_t i=0;i<xhat.size();++i) os<<xhat[i]<<" "; 
  os<<"\n";
  os.close();
}

void writeMMAP(const char* outfile, mex::vector<uint32_t> xhat) {
	ofstream os(outfile);
  os<<"MMAP\n";
  os<<maxVars.size();   // output # of MMAP variables, then var/val pairs
  for (size_t i=0;i<maxVars.size();++i) os<<" "<<maxVars[i]<<" "<<xhat[maxVars[i]];
  os<<"\n";
  os.close();
}


bool fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond=NULL);
double solveMBE(const graphModel& gm, const mex::VarOrder& order);
bool tryExactPR(const graphModel& gm, const mex::VarOrder& order);
bool gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond=NULL);

//Usage: solve_uai  -f <file.uai> -e <file.evid> -S <seed> -T <PR|MPE|MAR>
int main(int argc, char* argv[])
{

  double timeStart = timeSystem();

  const char* probName; // = argv[1];
  const char* taskName;//  = argv[3];
  Task task;
  mex::vector<Factor> bel;
	mex::vector<uint32_t> xhat;


  po::options_description desc("Available options");
  desc.add_options()
    ("help", "print help message")
    ("file,f", po::value<std::string>(), "input problem filename")
    ("evidence,e", po::value<std::string>(), "input evidence filename")
    ("query,q", po::value<std::string>(), "input marginal map query filename")
    ("seed,S", po::value<int>(),         "random number initial seed")
    ("task,T", po::value<std::string>(), "inference task string")
    ("ibound,i", po::value<int>(&iboundInit),       "initial i-bound")
    ("orders,o",    po::value<int>(&nOrders)->default_value(1),      "number of variable orderings to try")
    ("order-time,t",po::value<double>(&timeOrder)->default_value(1), "max time spend on variable orderings")
    ("order-rand",  po::value<int>(&nExtra)->default_value(0),   "var order randomness; n=-1 (none), or among best+n")
    ("order-limit", po::value<double>(&memUseRandom)->default_value(mex::infty()), "fail-out memory limit for ordering search")
    ("order-file", po::value<std::string>(), "problem elimination ordering filename")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("dt", po::value<double>(&dt)->default_value(300),         "file write time interval")
    ("lbps", po::value<double>(&lbpTime)->default_value(300),  "loopy belief propagation stop (seconds)")
    ("lbpi", po::value<double>(&lbpIter)->default_value(2000), "loopy belief propagation stop (iterations)")
    ("lbpe", po::value<double>(&lbpErr )->default_value(-1),   "loopy belief propagation stop (msg error)")
    ("lbpo", po::value<double>(&lbpObj )->default_value(-1),   "loopy belief propagation stop (objective)")
    ("gbps", po::value<double>(&gbpTime)->default_value(300),  "gen belief propagation stop (seconds)")
    ("gbpi", po::value<double>(&gbpIter)->default_value(-1),   "gen belief propagation stop (iterations)")
    ("gbpo", po::value<double>(&gbpObj )->default_value(-1),   "gen belief propagation stop (objective)")
//    ("cons", po::value<double>(&conTime)->default_value(300),  "conditional GBP stop (seconds)")
//    ("cono", po::value<double>(&conObj )->default_value(-1),   "conditional GBP stop (objective)")
    ("condition", po::value<int>(&doCond )->default_value(0),"do conditioning search after gbp (0=no,1=incremental,k=at most k states)")
    ("verbose", "verbose output during algorithm execution")
    ("ijgp", "use ijgp regions only")
  ;

	// Add positional defaults for basic UAI14 competition calls
	po::positional_options_description p;
	p.add("file", 1);
	p.add("evidence", 1);
	p.add("query", 1);
	p.add("task", 1);

  po::variables_map vm;
  //po::store(po::parse_command_line(argc,argv,desc),vm);
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm); 
  po::notify(vm);

  /*** ARGUMENT CHECKING *************************************************************/
  if (vm.count("help")) { std::cout<<desc<<"\n"; return 1; }
  if (vm.count("verbose")) doVerbose=true; else doVerbose=false;
  if (vm.count("file")) { probName=vm["file"].as<std::string>().c_str(); }
  else { std::cout<<"Missing input problem file!\n"; return 1; }
  if (vm.count("seed")) { mex::randSeed( vm["seed"].as<int>() ); }
  if (vm.count("task")) { taskName = vm["task"].as<std::string>().c_str(); task=Task(taskName); }
  else { std::cout<<"Missing task!\n"; return 1; }
  if (task==Task::MMAP && !vm.count("query")) { std::cout<<"Missing query for MMAP!\n"; return 1; }

  /*** UAI COMPETITION DEFAULTS BY TASK / RESOURCES **********************************/
	double TotalTime = 60.0;
	const char *inftime = ::getenv("INF_TIME");
	const char *infmem  = ::getenv("INF_MEMORY");
  if (infmem!=NULL && inftime!=NULL) {
		std::cout<<"Using default settings given UAI14 environment variables\n";
		if (vm["memory"].defaulted() && infmem!=NULL) MemLimit = atof(infmem)*1000-100;
  	if (!vm.count("seed")) mex::randSeed( 1234 );
		if (inftime!=NULL) TotalTime = atof(inftime);
		if (TotalTime <= 61.0) {									// Small time limit
			lbpTime = 1.5; lbpErr = 1e-4; 
			nOrders = 100; timeOrder = 1.5; nExtra = 2;
			gbpTime = 300; gbpObj = 1e-2; gbpIter = 100; dt = 0.5;
			iboundInit = 30; MemLimit = std::min(MemLimit, 300.0); memUseRandom=std::exp(28);
			doCond = 1;
		}	else if (TotalTime <= 1801) {					// Moderate time limit
			lbpTime = 10; lbpErr = 1e-4; 
			nOrders = 1000; timeOrder = 60; nExtra = 3;
			gbpTime = 1000; gbpObj = 1e-2; gbpIter = 100; dt = 30;
			iboundInit = 30; memUseRandom=std::exp(30);
			doCond = 1;
		} else {																// Large time limit
			lbpTime = 15; lbpErr = 1e-4; 
			nOrders = 1000; timeOrder = 100; nExtra = 3;
			gbpTime = 1000; gbpObj = 1e-2; gbpIter = 100; dt = 30;
			iboundInit = 30; memUseRandom=std::exp(40);
			doCond = 1;
		}
		//if (task==Task::PR) gbpIter=1;		// 07/12 version (results 07/13): no iterations
		//if (task==Task::MAR) gbpIter=2;		// seemed to work well on marginals (?)
	}

	std::cout<<"Memory limit set to "<<MemLimit<<"mb\n";

  doVerbose=true;

  /*** READ IN PROBLEM FILE **********************************************************/
  std::cout<<"Reading model file: "<<probName<<"\n";
  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> forig = Factor::readUai10(is);
	mex::graphModel gmo( forig );
	size_t nvar = gmo.nvar();  
  bel.resize(nvar);
	xhat.resize(nvar);


  /*** READ IN EVIDENCE FILE *********************************************************/
  VarSet evVar;
  ifstream is2;
  if (vm.count("evidence")) { is2.open(vm["evidence"].as<std::string>().c_str()); }
  if (is2.is_open()) {
    std::cout<<"Got evidence file "<<vm["evidence"].as<std::string>().c_str()<<"\n";
    std::map<uint32_t,size_t> evid;
    //int nEvid; is2 >> nEvid;
    int nEvid=1;			// 2014 format: single evidence, no # of evidences entry
    if (nEvid > 0) {
      int nEvidVar; is2 >> nEvidVar;
      for (size_t i=0;i<nEvidVar;i++) {
        uint32_t vid; size_t vval; is2>>vid>>vval; 
        evid[vid]=vval; evVar |= gmo.var(vid); 
				xhat[vid]=vval;
				bel[vid] = Factor::delta(gmo.var(vid), vval);		// write in belief (in case we forget later)
      }
			std::cout<<"Evidence on variables "<<evVar<<"\n";
			gmo.condition(evVar,xhat);
    }
  } else std::cout<<"Evidence file not specified or not found\n";


  /*** PREPARE REQUESTED TASK ********************************************************/
  std::cout<<"Task is "<<task;

	if (task==Task::MMAP) {
		std::string queryFile = vm["query"].as<std::string>(); 
		std::cout<<": query "<<queryFile;
		ifstream qis(queryFile.c_str()); 
		if (!qis.is_open()) { std::cout<<" does not exist!\n"; return 1; }
		int nMAP; qis >> nMAP;
		for (size_t i=0; i<nMAP; ++i) {
			size_t v; qis >> v;
			maxVars |= gmo.var(v);
		}
		std::cout<<": "<<maxVars.size()<<" query variables";
	}
	std::cout<<"\n";


  /*** PREPARE OUTPUT FILE ***********************************************************/
  std::string outfiles(probName); outfiles += '.'; outfiles += taskName;
  std::string::size_type start = outfiles.find_last_of('/');
  if (start==std::string::npos) start=0; else ++start;
  outfiles = outfiles.substr(start,std::string::npos);
  outfile = outfiles.c_str();
  //const char* outfile = outfiles.c_str();
  std::cout<<"Writing to "<<outfile<<"\n";

  double ln10 = std::log(10);


  /*** LOOPY BELIEF PROPAGATION ******************************************************/
	//mex::graphModel fg(flist);

	if (task==Task::PR || task==Task::MAR || task==Task::MMAP) {
  mex::mixbp _lbp(gmo.factors());  // flist
  std::cout<<"Model has "<<nvar<<" variables, "<<_lbp.nFactors()<<" factors\n";
  if (lbpIter != 0 && lbpTime > 0) {
		std::cout<<"Starting LBP...\n";
    _lbp.setProperties("Schedule=Priority,Distance=L1");
		if (task==Task::MMAP) {
			_lbp.setProperties("Schedule=Random");
			_lbp.setElimType(mex::mixbp::ElimType::Sum); 
			_lbp.setElimType(maxVars,mex::mixbp::ElimType::Max); 
		}
    _lbp.setStopIter(lbpIter); 
    _lbp.setStopMsg(lbpErr);
    _lbp.setStopObj(lbpObj);
    _lbp.setStopTime(lbpTime);
    _lbp.init();
  
    _lbp.run();
    switch (task) {
      case Task::PR: writePR(outfile,_lbp.logZ()); break;
      case Task::MAR: {
        for (size_t v=0;v<nvar;++v) bel[v]=_lbp.belief( _lbp.localFactor(v) );
        //for (size_t v=0;v<nvar;++v) if (!evVar.contains(Var(v,0))) bel[v]=_lbp.belief( _lbp.localFactor(v) );
        writeMAR(outfile,bel);
      } break;
			case Task::MMAP: {
				for (VarSet::const_iterator v=maxVars.begin();v!=maxVars.end();++v) xhat[*v]=_lbp.belief(*v).argmax();
			  writeMMAP(outfile,xhat);
			} break;
    }
    std::cout<<"LBP logZ="<<_lbp.logZ()/ln10<<"\n";

    _lbp.reparameterize();                                         // Convert loopy bp results to model
		gmo=graphModel(_lbp.factors()); 
  }
	} else if (task==Task::MPE) {
	mex::mplp _mplp(gmo.factors()); 
	if (lbpIter != 0 && lbpTime > 0) {
		std::cout<<"Starting MPLP...\n";
		_mplp.setProperties("Schedule=Fixed,Update=Var,StopObj=-1.0,StopMsg=-1.0");
  	_mplp.setStopIter(lbpIter); _mplp.setStopMsg(lbpErr); _mplp.setStopObj(lbpObj); _mplp.setStopTime(lbpTime);
  	_mplp.init();
  	_mplp.run();
		mex::vector<mex::index> best = _mplp.best(); 
		//for (size_t v=0;v<best.size();++v) std::cout<<best[v]<<" "; std::cout<<"\n";
		for (size_t v=0;v<best.size();++v) xhat[v]=best[v];
		std::cout<<"MPLP "<<gmo.logP(best)<<" ("<<_mplp.ub()<<")\n";
    if (gmo.logP(best) > -mex::infty()) writeMPE(outfile,best); 
  	gmo = graphModel(_mplp.beliefs()); 
    if (gmo.logP(xhat) == -mex::infty()) {
			std::cout<<"Trying Gibbs...\n";
			mex::gibbs _gibbs(gmo); 
			_gibbs.setProperties("Best=1,Beliefs=0"); _gibbs.init(); _gibbs.run();
			best = _gibbs.best();
			if (gmo.logP(best) >-mex::infty()) {
				std::cout<<"...got solution with logP="<<gmo.logP(best)<<"\n";
    		for (size_t v=0;v<best.size();++v) xhat[v]=best[v];
    		writeMPE(outfile,best);
			}
		}
	}
	}


  /*** HIGHER ORDERS ******************************************************/
  if ( gbpIter == 0 ) return 0;


  /*** BUILD JUNCTION GRAPH & ASSESS MEMORY USE ***************************/
  bool exact = false;
  size_t ibound = iboundInit, InducedWidth=10000;
  mex::VarOrder order; 				// start with a random order in case the model is very dense
	std::cout<<"Building elimination orders...\n";
	double score = gmo.order( mex::graphModel::OrderMethod::Random, order, 0, memUseRandom );
  if (MemLimit > 0) {

		const char *orderFile = NULL;				// Check for pre-specified elimination order file
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
    	size_t iOrder = 0;			
    	// Try to build new orders until time or count limit reached ////////////////////
    	while (iOrder < nOrders && (timeSystem()-startOrder < timeOrder)) {
    	  score = gmo.order(mex::graphModel::OrderMethod::WtMinFill, order, nExtra, score);
    	  ++iOrder;
   		}
			if (score < memUseRandom) InducedWidth = gmo.inducedWidth(order);
    	std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<", score "<<score<<"\n";

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

    // Try an exact solver and quit if it fits in memory
    if (task==Task::PR && tryExactPR(gmo,order)) return 0;   // TODO: Need to debug without this...

  } // end MemLimit check


  mex::gbp _gbp(gmo.factors());                      // Create a GBP object for later use
  if (vm.count("ijgp")) _gbp.setMinimal(true);      // use "ijgp-like" regions? (true = remove all with c=0)
  else                  _gbp.setMinimal(false);

  bool doneGBP=false, isExact=false;
	if (doVerbose) std::cout<<"\n"<<"Beginning basic GBP...\n";
  while (!doneGBP) { 
   if (MemLimit > 0) { 
     _gbp.clearRegions();
     isExact = gbpPopulateCliques(_gbp,order,ibound,NULL);
   } else {
    // Didn't build MBE => no way to check ibound
    ibound = 0; isExact=false;
   }

  try {
  // Run GBP on the region graph 
  std::cout<<"GBP with "<<_gbp.nRegions()<<" regions; mem "<<_gbp.memory()<<"M\n";
	if (task==Task::MPE) _gbp.setTask(mex::gbp::Task::Max);
  _gbp.setProperties("Schedule=Fixed");
  _gbp.init();
	if (task==Task::MPE) _gbp.setBest( mex::vector<mex::index>(xhat.begin(),xhat.end()) );
  _gbp.setStopIter(gbpIter); _gbp.setStopObj(gbpObj); _gbp.setStopMsg(-1.0); 
  _gbp.setVerbose(doVerbose);
	if (isExact) _gbp.setDamping(-1.0);  		// no damping if we think it's exact
	if (task==Task::MMAP || task==Task::MPE) { 
		_gbp.setDamping(-1.0);  							// no damping for hacky proximal mmap or mpe
		_gbp.setStopObj(-1.0); gbpObj=-1.0;  	// also, no objective stopping (since proximal will update)
	}

  // Get region indices for single-variable beliefs
  mex::vector<mex::gbp::findex> regions(nvar);
  for (size_t v=0;v<nvar;++v) regions[v]=_gbp.regionWith(Var(v,0));

  double gbpLeft, gbpStop = timeSystem()+gbpTime;
  while ( (gbpLeft = gbpStop - timeSystem()) > 0 ) {
    _gbp.setStopTime( std::min( dt , gbpLeft ) );
    _gbp.run();
    switch (task) {
      case Task::PR: writePR(outfile,_gbp.logZStable()); break; 
      case Task::MAR: {
        //for (size_t v=0;v<nvar;++v) bel[v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
        for (size_t v=0;v<nvar;++v) bel[v]=_gbp.beliefs()[v]; //(Var(v,0));
        writeMAR(outfile,bel);
      } break;
			case Task::MPE: {
				if (_gbp.lb() > - mex::infty()) {
					mex::vector<mex::index> best = _gbp.best();
					writeMPE(outfile,best);
				}
			} break;
			case Task::MMAP: {
				for (VarSet::const_iterator v=maxVars.begin();v!=maxVars.end();++v) {
					//bel[*v] = _gbp.belief(*v);
					bel[*v] = _gbp.computeRegionBelief(regions[*v]);
					ind2sub(bel[*v].vars(), bel[*v].argmax(), xhat);
				}
        writeMMAP(outfile,xhat);
				for (VarSet::const_iterator v=maxVars.begin();v!=maxVars.end();++v) {
					_gbp.__factor(regions[*v]) += log(bel[*v]);			// proximal: remove entropy (subtract -log(b)) from factors
					// if GBP is not using GBP_LOG (???) should multiply regions[v] factor by bel[v] 
					// !!!! Note, this is not the right way to do this; it should be done in the gbp object...
				}
      } break;
    }
    std::cout<<"GBP logZ "<<_gbp.logZ()/ln10<<"\n";
    if (_gbp.dObj() < gbpObj) { std::cout<<"Reached objective tolerance\n"; break; }
    if (_gbp.iter() >= gbpIter && gbpIter > 0) { std::cout<<"Reached iteration limit\n"; break; }
    if (_gbp.logZ() == -mex::infty()) { std::cout<<"Model deemed inconsistent\n"; break; }
		
  }
  doneGBP = true;

  if (isExact && _gbp.dObj()<gbpObj) { //InducedWidth <= ibound)
    std::cout<<"Answer should be exact\n";
    return 0;
  }

  } catch (std::exception& e) {
		std::cout<<"Exception: "<<e.what()<<"\n";
		if (_gbp.getDamping() > 0.25) {
			doneGBP=true; std::cout<<"Caught exception (memory overreach?).  Damping on => quitting GBP\n";
		} else {
    	doneGBP=false; MemLimit*=.9; if (ibound) ibound--;
			std::cout<<"Caught exception (memory overreach?).  Trying again with ibound "<<ibound<<" and MemLimit "<<MemLimit<<"\n";
		}
  }
  }

  // If we do not proceed to conditioning search, just exit
  if (!doCond) {
    std::cout<<"Quitting after GBP\n";
    return 0;
  }


  /*** ITERATIVE CONDITIONING AND GBP *************************************/
	// More general way to do this: condition-min (=0), condition-max (=inf)   (# states?)
	// flag to decide whether to discard if non-convergent (?) (see this vs cgbp) (or ok with damping?)
	// decide how much time?  gbpTime, or as a fraction of # states?
	// =1 => previous code; > 1 is subsequent code
	// need to include exact lnZ check here also

	if (task==Task::MMAP || task==Task::MPE) { std::cout<<"Conditioning not supported for MPE/MMAP\n"; return 0; }
	if (doVerbose) std::cout<<"\n"<<"Beginning conditioned GBP...\n";
	if (order.size()==0) order=gmo.order(mex::graphModel::OrderMethod::MinWidth);  // need an order if none yet...
 
	_gbp = mex::gbp( mex::vector<Factor>() );  // !!! blank out GBP object; restore memory
  VarSet cond;

	//!!! TODO: if doCond > 1, we should immediately populate with that many states

  while (!isExact) {
    // add check for best conditioner given cardinality limit !!!  or, only increment if anytime...
    cond += gmo.bestConditioner(order,cond);

		if (doVerbose) std::cout<<"\n";
    std::cout<<"Conditioning "<<cond<<"\n";
		ibound = iboundInit;																// check all iBounds again (in case higher available)

  	bool doneCGBP=false;
  	while (!doneCGBP) {
			bool useMBE=false;
      if (task == Task::PR && fitsMBE(gmo,order,&cond)) {
        std::cout<<"Trying exact via MBE\n";
        useMBE=true; isExact=true;
      }
    	try {
    	Factor lnZ(cond);
   		bool failed = false;
    	mex::vector<mex::vector<Factor> > condMarginals;
     	mex::vector<mex::gbp::findex> regions(nvar);			// !!! dangerous; pull out & match to prev (clear+resize)
    	for (size_t i=0;i<lnZ.nrStates();++i) {
				mex::graphModel gmcond(gmo); 
				std::vector<size_t> val(nvar); ind2sub(cond,i,val);
				gmcond.condition(cond,val);
      	mex::vector<Factor> fcond = gmcond.factors();
/*
      	std::map<Var,size_t> val;  ind2sub(cond,i,val);
      	mex::vector<Factor> fcond = gmo.factors();
      	for (size_t f=0;f<fcond.size();++f) {
        	VarSet isect = cond & fcond[f].vars();
        	if (isect.size() > 0) fcond[f] = fcond[f].condition(isect, sub2ind(isect,val));
      	}
				for (size_t v=0;v<cond.nvar();++v) fcond.push_back(Factor::delta(cond[v],sub2ind(cond[v],val)));
*/

        try { if (useMBE) {         // if a bucket elim pass was good enough, do that:
          mex::graphModel gm(fcond);
          lnZ[i] = solveMBE(gm,order);
          for (size_t v=0;v<cond.size();++v) std::cout<<cond[v]<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
          continue;
        }                     // otherwise we need to do GBP-like updates:
        } catch (std::exception& e) {
					std::cout<<"Exception: "<<e.what()<<"\n";
          std::cout<<"Caught exception; failure in MBE; trying GBP\n";
        }

 	     	mex::lbp fgcond(fcond); fgcond.setProperties("Schedule=Priority,Distance=L1");
 	     	fgcond.setStopIter(lbpIter); fgcond.setStopMsg(lbpErr); fgcond.setStopObj(lbpObj);
 	     	fgcond.setStopTime(0.5);
 	    	fgcond.init();
 	     	fgcond.run();
 	     	fgcond.reparameterize();
	     	fcond = fgcond.factors();

				// TODO:DEBUG: it seems that reinitializing gbp is not equivalent to re-constructing a new one.. FIX
      	mex::gbp _gbp(fcond);
    		if (vm.count("ijgp")) _gbp.setMinimal(true); else _gbp.setMinimal(false);  // use "ijgp-like" regions?
   			if (MemLimit > 0) isExact = gbpPopulateCliques(_gbp,order,ibound,&cond);
				else              isExact = false;
      	_gbp.setProperties("Schedule=Fixed");
      	_gbp.setStopIter(gbpIter); _gbp.setStopObj(gbpObj); _gbp.setStopMsg(-1.0);
      	_gbp.setStopTime(gbpTime); _gbp.setVerbose(false);
      	_gbp.init();

      	if (task==Task::MAR) for (size_t v=0;v<nvar;++v) regions[v]=_gbp.regionWith(Var(v,0));
      	// uncomment duplicate lines below

// !!! make incremental on last output round
//      _gbp.setFactors(fcond); _gbp.setProperties("Schedule=Fixed");
//      _gbp.init();

      	_gbp.run();
      	//if (_gbp.dObj() >= gbpObj) { failed=true; break; }     // convergence failure on this condition
				//TODO: never worry about this?  or, only worry if anytime enabled?
      	lnZ[i] = _gbp.logZStable();
      	if (task==Task::MAR) {
        	condMarginals.push_back(bel);
        	//for (size_t v=0;v<nvar;++v) condMarginals[i][v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
        	for (size_t v=0;v<nvar;++v) condMarginals[i][v]=_gbp.beliefs()[v]; //(Var(v,0));
      	}
      	for (size_t v=0;v<cond.size();++v) std::cout<<cond[v]<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
    	}
    	if (failed) {std::cout<<"Failing out\n"; doneCGBP=true; continue;} // if we failed out, condition on more vars
    	
			doneCGBP=true;
    	double lnZtot = lnZ.logsumexp();
    	switch (task) {
      	case Task::PR:
        	writePR(outfile, lnZtot);
        	break;
      	case Task::MAR:
        	Factor probs = (lnZ - lnZtot).exp();
        	for (size_t v=0;v<nvar;++v) {
          	if (evVar.contains(Var(v,0))) { } // evidence variables not updated
          	else if (cond.contains(Var(v,0)))  { bel[v] = probs.marginal(Var(v,0)); }
          	else {		// TODO: change to incremental update?
            	bel[v] = condMarginals[0][v] * probs[0];
            	for (size_t i=1;i<lnZ.nrStates();++i) bel[v] += condMarginals[i][v] * probs[i];
          	}
        	}
        	writeMAR(outfile, bel);
        	break;
    	}
    	std::cout<<"Conditioning "<<cond<<" => lnZ "<<lnZtot<<" (logZ "<<lnZtot/ln10<<")\n";

    	} catch (std::exception& e) {
				std::cout<<"Exception: "<<e.what()<<"\n";
      	doneCGBP=false; MemLimit*=.9; if (ibound) ibound--;
      	std::cout<<"Caught exception (memory overreach?).  Trying again with ibound "<<ibound<<" and MemLimit "<<MemLimit<<"\n";
      	continue;
				// TODO: this is right if we're only doing this part, but not right if we're running incrementally
    	}

			// !!! TODO: if doCond > 1, quit (non-incremental)?
    	if (isExact) { std::cout<<"Answer should be exact\n"; return 0; }
  	}
  }

  // run CCS expansion on larger regions until timeout (?)

  return 0;

}


bool fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond) {
  double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
  bool isExact = true;
  mex::mbe mb(gm.factors());
  mb.setOrder(order);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
  mb.setIBound(100); double mbMem = mb.simulateMemory(NULL,cond,mbCutoff,&isExact);
  if (mbMem < mbCutoff && isExact) return true;
  return false;
}

double solveMBE(const graphModel& gm, const mex::VarOrder& order) {
  mex::mbe mb(gm.factors());
  mb.setOrder(order);
  mex::vector<mex::mbe::vindex> pt(order.size()); for (size_t v=0;v<pt.size();++v) pt[v]=-1;
  mb.setPseudotree(pt);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
  mb.setIBound(100); //double mbMem = mb.simulateMemory(NULL,NULL,mbCutoff,&isExact);
  std::cout<<"Attempting exact solve\n";
  //std::cout<<"Attempting exact solve; mbMem="<<mbMem<<" vs "<<mbCutoff<<" ("<<MemLimit<<")\n";
  mb.init();
  return mb.logZ();
}


bool tryExactPR(const graphModel& gm, const mex::VarOrder& order) {
  try {
    double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    bool isExact = true;
    mex::mbe mb(gm.factors());
    mb.setOrder(order);
    mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
    mb.setIBound(100); double mbMem = mb.simulateMemory(NULL,NULL,mbCutoff,&isExact);
    if (mbMem < mbCutoff && isExact) {
      std::cout<<"Attempting exact solve; mbMem="<<mbMem<<" vs "<<mbCutoff<<" ("<<MemLimit<<")\n";
      mb.init();
      writePR(outfile,mb.logZ());
      std::cout<<"Exact solution by MBE: "<<mb.logZ()/std::log(10)<<"\n";
      return true;
    }
  } catch (std::exception& e) {
    // Failed (probably for memory reasons) => try GBP
    std::cout<<"Failed (due to memory problem?)  Trying GBP\n";
  }
  return false;
}


bool gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond) {
    bool isExact = false;
    double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    _gbp.clearRegions();
    mex::mbe mb(_gbp.gmOrig().factors());
    mb.setOrder(order);
    mex::vector<mex::mbe::vindex> pt(order.size()); for (size_t v=0;v<pt.size();++v) pt[v]=-1;
    mb.setPseudotree(pt);
    mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=1,DoHeur=0");
    mb.setIBound(ibound);
    mex::vector<VarSet> cliques;
    double mem = std::numeric_limits<double>::infinity();
    mbCutoff *= 2;  // leave a little slack for mismatch in criteria
    double mbMem = mb.simulateMemory(&cliques,cond, mbCutoff, &isExact);
    if (mbMem < mbCutoff) { _gbp.addRegions(cliques); mem = _gbp.memory(); }
    while ( ibound > 0 && (mbMem >= mbCutoff || mem > MemLimit) ) {
      std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
      mb.setIBound(--ibound); cliques.clear(); mbMem=mb.simulateMemory(&cliques,cond, mbCutoff, &isExact);
      if (mbMem < mbCutoff) { _gbp.clearRegions(); _gbp.addRegions(cliques); mem=_gbp.memory(); }
    }
		//ofstream ofs("cliques.mbe.txt");
		//for (size_t c=0;c<cliques.size();++c) ofs<<cliques[c]<<"\n"; std::cout<<"\n";  // output for DEBUG !!!
		//ofs.close();
    std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
    return isExact;
}




/*
void condition(it first, it last, varset vs, maptype val) {
  for (;first != last; ++first) { VarSet vf=first->vars()&vs; *first = first->condition(vf, sub2ind(vf,val)); }
}
*/

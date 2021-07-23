#include <time.h>
#include <Utils.h>
#include <IoData.h>
#include <VarFcnSG.h>
#include <VarFcnMG.h>
#include <VarFcnJWL.h>
#include <ExactRiemannSolverBase.h>
#include <set>
using std::cout;
using std::endl;
/*************************************
 * Main Function
 ************************************/
int verbose = 0;
int main(int argc, char* argv[])
{
  clock_t start_time = clock(); //for timing purpose only

  //! Initialize PETSc and MPI 
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m                 START                    \033[0m\n"); 
  print("\033[0;32m==========================================\033[0m\n");
  print("\n");


  //! Read user's input file
  IoData iod(argc, argv);

  verbose = iod.output.verbose;

  //! Initialize VarFcn (EOS, etc.) 

  std::vector<VarFcnBase *> vf;
  for(int i=0; i<iod.eqs.materials.dataMap.size(); i++)
    vf.push_back(NULL); //allocate memory for the VarFcn pointers

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= vf.size()) {
      print_error("Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    if(it->second->eos == MaterialModelData::STIFFENED_GAS)
      vf[matid] = new VarFcnSG(*it->second);
    else if(it->second->eos == MaterialModelData::MIE_GRUNEISEN)
      vf[matid] = new VarFcnMG(*it->second);
    else if(it->second->eos == MaterialModelData::JWL)
      vf[matid] = new VarFcnJWL(*it->second);
    else {
      print_error("Error: Unable to initialize variable functions (VarFcn) for the specified material model.\n");
      exit_mpi();
    }
  }


  ExactRiemannSolverBase riemann(vf, iod.exact_riemann);
  
  double Vm[5], Vp[5], V[5];
  int idm, idp;
  Vm[0] = iod.bc.inlet.density;
  Vm[1] = iod.bc.inlet.velocity_x;
  Vm[2] = 0.0;
  Vm[3] = 0.0;
  Vm[4] = iod.bc.inlet.pressure;
  idm   = iod.bc.inlet.materialid;
  Vp[0] = iod.bc.outlet.density;
  Vp[1] = iod.bc.outlet.velocity_x;
  Vp[2] = 0.0;
  Vp[3] = 0.0;
  Vp[4] = iod.bc.outlet.pressure;
  idp   = iod.bc.outlet.materialid;

  print("Solving a One-Dimensional Riemann Problem...\n");
  print("Left  State: %e %e %e (MaterialID: %d).\n",Vm[0],Vm[1],Vm[4],idm);
  print("Right State: %e %e %e (MaterialID: %d).\n",Vp[0],Vp[1],Vp[4],idp);


  if(argc==5) {//plot p-u relation
    double pmin = atof(argv[2]);
    double pmax = atof(argv[3]);
    double dp   = atof(argv[4]);
    riemann.PrintStarRelations(Vm[0], Vm[1], Vm[4], idm, Vp[0], Vp[1], Vp[4], idp, pmin, pmax, dp);
    print("Printed the star state relations.\n");
  }

  int id;
  double Vsm[5], Vsp[5];
  riemann.ComputeRiemannSolution(0, Vm, idm, Vp, idp, V, id, Vsm, Vsp);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m           NORMAL TERMINATION             \033[0m\n"); 
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  for(int i=0; i<vf.size(); i++)
    delete vf[i];

  return 0;
}

//--------------------------------------------------------------

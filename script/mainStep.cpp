// =================================================================== //
// this is file mainStep.cpp, created at 10-Jun-2007                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "Simulator3D.h"
#include "TElement.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "Periodic3D.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 double Re = 100;
 double Sc = 200;
 double Fr = 10;
 //double alpha = 1;
 double cfl = 0.1;
 double mu_l = 1.0;
 double rho_l = 1.0;
 //Solver *solverP = new PCGSolver();
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 //Solver *solverP = new PetscSolver(KSPPREONLY,PCLU);
 //Solver *solverP = new PetscSolver(KSPLSQR,PCILU);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 //string meshFile = "rectangle.msh";

 //const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk-step/";
 //const char *datFolder  = "./dat/";
 //string meshDir = (string) getenv("DATA_DIR");
 //meshDir += "/gmsh/3d/singlePhase/" + meshFile;
 //const char *mesh = meshDir.c_str();

 //const char *mesh = "/Users/peixoto/meshes/3d/periodic-cylindrical-3d-Lx6.vtk";
 //const char *mesh = "/Users/peixoto/meshes/3d/periodic-rectangular-3d-channel-Lx2-Ly1-Ly1-23verts.vtk";
 const char *mesh = "/Users/peixoto/meshes/3d/periodic-rectangular-3d-channel-Lx5-Ly1-Ly1-1254verts.vtk";

 Model3D m1;

 //m1.setMeshStep(40,20,2);
 //m1.setAdimenStep();
 //m1.readMSH(mesh);
 m1.readVTK(mesh);
 m1.setInterfaceBC();
 //m1.setTriEdge();
 //m1.mesh2Dto3D();
 m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setVertNeighbour();
 m1.setInOutVert();
 //m1.setStepBC();
 m1.setStepPBC();
 m1.setCStepBC();

 Periodic3D pbc(m1);
 Simulator3D s2(pbc,m1);
 pbc.MountPeriodicVectors(m1);
 //Simulator3D s1(m1);

 s2.setRe(Re);
 s2.setSc(Sc);
 s2.setFr(Fr);
 s2.setCfl(cfl);
 s2.setDtEulerian();
 s2.setMu(mu_l);
 s2.setRho(rho_l);
 s2.setSolverPressure(solverP);
 s2.setSolverVelocity(solverV);
 s2.setSolverConcentration(solverC);

 s2.init();
 s2.assembleSlip();
 s2.matMount();
 s2.matMountC();
 s2.setUnCoupledBC(); 
 //s2.setUnCoupledCBC(); 

 if( (*(argv+1)) == NULL )
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;
 }
 else if( strcmp( *(argv+1),"restart") == 0 )
 {
  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;

  string file = (string) "sim-" + *(argv+2);
  iter = s2.loadSolution("./","sim",atoi(*(argv+2)));
 }

 InOut save(m1,s2); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);

 int nIter = 100;
 int nRe = 5;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << iter << endl;

   s2.stepSL();
   //s2.setCRHS();
   s2.setUnCoupledBC(); 
   s2.setGravity("+X");
   s2.setRHS_PBC();
   s2.setCopyDirectionPBC("RL");
   s2.unCoupledPBC();
   s2.inputVelocityPBC();
   //s2.unCoupledC();

   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTU(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);

   s2.saveOldData();

   s2.timeStep();

   cout << "________________________________________ END of "
	    << iter << endl;

   iter++;
  }
 }

 PetscFinalize();
 return 0;
}

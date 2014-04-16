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

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 double Re = 10000;
 double Sc = 2000;
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

 //const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "/home/gcpoliveira/post-processing/vtk/3d/step/";
 //const char *datFolder  = "./dat/";

 string meshFile = "cuboid-3d-w0.2.msh";

 string meshDir = (string) getenv("MESH3D_DIR");
 meshDir += "/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 //m1.setMeshStep(40,20,2);
 //m1.setAdimenStep();
 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge();
 m1.mesh2Dto3D();
 m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 //m1.setNeighbour(); // error if called together with setVertNeighbour()
 m1.setVertNeighbour();
 m1.setInOutVert();
 m1.setStepBC();
 m1.setCStepBC();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setFr(Fr);
 s1.setCfl(cfl);
 s1.setDtEulerian();
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 s1.init();
 s1.assembleSlip();
 s1.matMount();
 s1.matMountC();
 s1.setUnCoupledBC(); 
 s1.setUnCoupledCBC(); 

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
  iter = s1.loadSolution("./","sim",atoi(*(argv+2)));
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);

 int nIter = 10;
 int nRe = 1;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << iter << endl;

   s1.stepSL();
   s1.setRHS();
   s1.setCRHS();
   s1.unCoupled();
   s1.unCoupledC();

   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTU(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);

   s1.saveOldData();

   s1.timeStep();

   cout << "________________________________________ END of "
	    << iter << endl;

   iter++;
  }
 }

 PetscFinalize();
 return 0;
}

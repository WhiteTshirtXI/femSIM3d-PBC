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
#include "Helmholtz3D.h"
#include "PetscSolver.h"
#include "petscksp.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 real Re = 10000;
 real Sc = 2000;
 real Fr = 10;
 //real alpha = 1;
 //real beta = 0;
 real c1 = 0.3;  // lagrangian
 real c2 = 0.0;  // smooth vel
 real c3 = 0.2;  // smooth coord (fujiwara)
 real d1 = 0.0;  // surface tangent velocity u_n=u-u_t 
 real d2 = 0.00;  // surface smooth cord (fujiwara)
 real cfl = 1;
 real mu_l = 1.0;
 real rho_l = 1.0;
 //Solver *solverP = new PCGSolver();
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 //Solver *solverP = new PetscSolver(KSPPREONLY,PCLU);
 //Solver *solverP = new PetscSolver(KSPLSQR,PCILU);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 string meshFile = "stepSimple.msh";

 //const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 //m1.setMeshStep(40,20,4);
 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge();
 m1.checkTriangleOrientation();
 m1.mesh2Dto3D();
 //m1.setAdimenStep();
 //m1.setSingleElement();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setOFace();
 m1.setVertNeighbour();
 m1.setInOutVert();
 m1.setMapEdge();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();
 m1.tetMeshStats();
 m1.setStepBC();
 m1.setCStepBC();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setFr(Fr);
 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setD1(d1);
 s1.setD2(d2);
 s1.setCfl(cfl);
 s1.setDtEulerian();
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 s1.init();

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
  iter = s1.loadSolution("sim",atoi(*(argv+2)));
 }
 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initSquareChannel();
 h1.assemble();
 h1.setk(0.7);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 //h1.saveVTK("./vtk/","edge");
 h1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);

 int nIter = 1000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   s1.setDtALESinglePhase();

   save.printSimulationReport();

   //s1.stepSL();
   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.movePoints();
   s1.assembleSlip();
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();
   s1.setUnCoupledCBC();
   s1.setRHS();
   s1.setCRHS();
   s1.unCoupled();
   s1.unCoupledC();

   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTU(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Helmholtz3D h2(m1,h1);
  h2.setBC();
  h2.initSquareChannel();
  h2.assemble();
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();
  h2.saveVTK("./vtk/","edge",iter-1);
  h2.setModel3DEdgeSize();

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume(0.005);
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.delete3DPoints();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setVertNeighbour();
  m1.setInOutVert();
  m1.setMapEdge();
  m1.setSurfaceVolume();
  m1.setSurfaceArea();
  m1.setStepBC();
  m1.setCStepBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveVTK(vtkFolder,"sim",iter-1);
  saveEnd.saveSol(binFolder,"sim",iter-1);
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}



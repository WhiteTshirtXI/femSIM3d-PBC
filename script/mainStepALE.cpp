// =================================================================== //
// this is file mainStep.cpp, created at 10-Jun-2007                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "TElement.h"
#include "GMRes.h"
#include "InOut.h"
#include "Helmholtz3D.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 int iter = 1;
 real Re = 1000;
 real Sc = 200;
 real Fr = 10;
 //real alpha = 1;
 real c1 = 0.3;  // lagrangian
 real c2 = 0.0;  // smooth vel
 real c3 = 2.0;  // smooth coord (fujiwara)
 real cfl = 1.0;
 real mu_l = 1.0;
 real rho_l = 1.0;
 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverP = new PetscSolver(KSPBICG,PCJACOBI);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 //string meshFile = "backwardStep.msh";
 string meshFile = "backwardStepHole.msh";
 //string meshFile = "retangle.msh";

 //const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/singlePhase/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 if( (*(argv+1)) == NULL )
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  //m1.setMeshStep(40,20,4);
  m1.readMSH(mesh);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();

  // mesh statistics info 
  m1.setInitSurfaceVolume();
  m1.setSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setSurfaceArea();
  m1.tetMeshStats();

  // boundary conditions
  m1.setGenericBC();

  s1(m1);

  s1.setRe(Re);
  s1.setSc(Sc);
  s1.setFr(Fr);
  s1.setC1(c1);
  s1.setC2(c2);
  s1.setC3(c3);
  s1.setCfl(cfl);
  s1.setDtALESinglePhase();
  s1.setMu(mu_l);
  s1.setRho(rho_l);
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  s1.init();
 }
 else if( strcmp( *(argv+1),"restart") == 0 )
 {
  // load surface mesh
  string aux = *(argv+2);
  string file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();

  s1(m1);

  // load 3D mesh
  file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setGenericBC();

  s1(m1);

  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  iter = s1.loadSolution("./","sim",atoi(*(argv+2)));
 }
 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initBackwardStep();
 h1.assemble();
 h1.setk(2.0);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 //h1.saveVTK(vtkFolder,"edge",0);
 h1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 10000;
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

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   //s1.stepLagrangian();
   s1.stepALE();
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

   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);

   s1.saveOldData();

   s1.timeStep();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Helmholtz3D h2(m1,h1);
  h2.setBC();
  h2.initBackwardStep();
  h2.assemble();
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();
  h2.saveVTK(vtkFolder,"edge",iter-1);
  //h2.saveChordalEdge(datFolder,"edge",iter-1);
  h2.setModel3DEdgeSize();

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  m1.initMeshParameters();

  // 3D operations
  //m1.insert3dMeshPointsByVolume();
  m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight(); // pressure points
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

  // mesh statistics info 
//--------------------------------------------------
//   m1.setInitSurfaceVolume();
//   m1.setSurfaceVolume();
//   m1.setInitSurfaceArea();
//   m1.setSurfaceArea();
//   m1.tetMeshStats();
//-------------------------------------------------- 
  m1.setSurfaceConfig();

  // boundary conditions
  m1.setGenericBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}



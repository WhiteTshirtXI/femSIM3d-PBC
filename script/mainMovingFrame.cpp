// =================================================================== //
// this is file mainRisingBubble.cpp, created at 10-Jun-2009           //
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
#include "Periodic3D.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 // bogdan's thesis 2010 (Bhaga and Weber, JFM 1980)
 int iter = 1;
 double vel = 0;
 double c1 = 0.0;      // lagrangian
 double c2 = 0.0;      // smooth vel
 double c3 = 3.0;      // smooth coord (fujiwara)
 double d1 = 1.0;      // surface tangent velocity u_n=u-u_t 
 double d2 = 0.1;      // surface smooth cord (fujiwara)
 double alpha = 1.0;   // time discrete method

 /* Japonese paper */
 double Re = 571.732; 
 double We = 13.6168;
 double Fr = 31.9275;
 double mu_in=1.2471e-05;
 double mu_out=0.00020235;
 double rho_in=26.016;
 double rho_out=1156.9;

 double cfl = 0.8;

 //string meshFile = "circular.msh";
 
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 //Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 //string meshDir = (string) getenv("DATA_DIR");
 //meshDir += "/gmsh/3d/micro/" + meshFile;
 //const char *mesh = meshDir.c_str();

 //const char *mesh = "/Users/peixoto/meshes/3d/annular-3d.msh";
 const char *mesh = "/Users/peixoto/meshes/3d/circular.msh";
 //const char *mesh = "/Users/peixoto/meshes/gustavo-anjos-meshes/airWater.msh";

 Model3D m1;
 Simulator3D s2;
 Periodic3D pbc;

 if( *(argv+1) == NULL )     
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  const char *mesh1 = mesh;

  m1.readMSH(mesh1);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();
  m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  //m1.setGenericBC();
  m1.setBubbleArrayPeriodicBC(); // <<<
  pbc.MountPeriodicVectors(m1);
  
  s2(m1);

  s2.setRe(Re);
  s2.setWe(We);
  s2.setFr(Fr);
  s2.setC1(c1);
  s2.setC2(c2);
  s2.setC3(c3);
  s2.setD1(d1);
  s2.setD2(d2);
  s2.setAlpha(alpha);
  s2.setMu(mu_in,mu_out);
  s2.setRho(rho_in,rho_out);
  s2.setCfl(cfl);
  s2.init();
  s2.setDtALETwoPhase();
  s2.setSolverPressure(solverP);
  s2.setSolverVelocity(solverV);
  s2.setSolverConcentration(solverC);
 }
 else if( strcmp( *(argv+1),"restart") == 0 ) 
 {
  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;

  // load surface mesh
  string aux = *(argv+2);
  string file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();

  s2(m1);

  // load 3D mesh
  file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.readVTKHeaviside(vtkFile);
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  //m1.setGenericBC();
  m1.setBubbleArrayPeriodicBC(); // <<<
  
  s2(m1);

  s2.setSolverPressure(solverP);
  s2.setSolverVelocity(solverV);
  s2.setSolverConcentration(solverC);

  iter = s2.loadSolution("./","sim",atoi(*(argv+2)));
 }
 else if( strcmp( *(argv+1),"remesh") == 0 ) 
 {
  cout << endl;
  cout << "--------------> RE-MESHING & STARTING..." << endl;
  cout << endl;

  // load old mesh
  Model3D mOld;
  string file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();
  mOld.readVTK(vtkFile);
  mOld.readVTKHeaviside(vtkFile);
  mOld.setMapping();

  // load surface mesh and create new mesh
  file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3DOriginal();
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  //m1.setGenericBC();
  m1.setBubbleArrayPeriodicBC(); // <<<

  s2(m1);

  s2.setSolverPressure(solverP);
  s2.setSolverVelocity(solverV);
  s2.setSolverConcentration(solverC);
  iter = s2.loadSolution("./","sim",atoi(*(argv+2)));
  s2.applyLinearInterpolation(mOld);
 }
 else if( strcmp( *(argv+1),"restop") == 0 )  
 {
  cout << endl;
  cout << "--------------> RE-MESHING (NO ITERATION)..." << endl;
  cout << endl;

  // load old mesh
  Model3D mOld;
  string file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();
  mOld.readVTK(vtkFile);
  mOld.readVTKHeaviside(vtkFile);
  mOld.setMapping();

  // load surface mesh and create new one
  file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3DOriginal();
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();

  s2(m1);
  //file = (string) "sim-" + *(argv+2);
  //const char *sol = file.c_str();
  iter = s2.loadSolution("./","sim",atoi(*(argv+2)));
  s2.applyLinearInterpolation(mOld);

  InOut saveEnd(m1,s2); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveMSH(mshFolder,"newMesh",atoi(*(argv+2)));
  saveEnd.saveSol(binFolder,"sim",atoi(*(argv+2)));
  //saveEnd.saveVTKSurface(vtkFolder,"sim",atoi(*(argv+2)));
  return 0;
 }
 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initRisingBubble();
 h1.assemble();
 h1.setk(0.2);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 h1.setModel3DEdgeSize();

 InOut save(m1,s2); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 double vinst=0;
 double vref=0;
 int nIter = 30;
 int nReMesh = 1;

 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   vinst = s2.getCentroidVelXAverage();
   vref += vinst;
   cout << vref << " " << vinst << endl;
   s2.setUSol(vinst);
   m1.setGenericBC(vref);
   s2.setURef(vref);

   s2.setDtALETwoPhase();

   InOut save(m1,s2); // cria objeto de gravacao
   save.printSimulationReport();

   //s2.stepLagrangian();
   //s2.stepALE();
   s2.stepSLPBCFix(); // <<< 
   s2.movePoints();
   //s2.assemble();
   s2.assemblePBC(); // <<<
   s2.matMount();
   //s2.setUnCoupledBC();
   s2.setUnCoupledPBC(); // <<<
   //s2.setRHS();
   s2.setRHS_PBC(); // <<<
   //s2.setGravity("-Z");
   //s2.setInterface();
   s2.setInterfaceGeo();
   //s2.unCoupled();
   s2.unCoupledPBC(); // <<<

   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s2.saveOldData();

   s2.timeStep();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Helmholtz3D h2(m1,h1);
  h2.setBC();
  h2.initRisingBubble();
  h2.assemble();
  h2.setk(0.2);
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();
  h2.saveChordalEdge(datFolder,"edge",iter-1);
  h2.setModel3DEdgeSize();

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();
  m1.initMeshParameters();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints();

  // surface operations
  m1.smoothPointsByCurvature();

  m1.insertPointsByLength("flat");
  //m1.insertPointsByCurvature("flat");
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance("flat");
  m1.contractEdgesByLength("flat");
  //m1.removePointsByLength();
  m1.flipTriangleEdges();

  m1.removePointsByNeighbourCheck();
  //m1.checkAngleBetweenPlanes();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
  m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInterfaceBC();
  m1.setGenericBC(vref);

  Simulator3D s3(m1,s2);
  s3.applyLinearInterpolation(mOld);
  s2 = s3;
  s2.setSolverPressure(solverP);
  s2.setSolverVelocity(solverV);
  s2.setSolverConcentration(solverC);

  InOut saveEnd(m1,s2); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}



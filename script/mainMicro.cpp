// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "InOut.h"
#include "Laplace3D.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 // set each bubble length
 vector< real > triEdge;
 triEdge.resize(6);
 triEdge[0] = 0.1;   // none
 triEdge[1] = 0.24;   // wall
 triEdge[2] = 0.09;   // bubble 1 
 triEdge[3] = 0.09;   // bubble 2 
 triEdge[4] = 0.022; // bubble 3
 triEdge[5] = 0.022; // bubble 4

 // bogdan's thesis 2010 - case 2
 int iter = 0;
 real Re = 10;
 real Sc = 1;
 real We = 2;
 real Fr = 1.0;
 real c1 = 0.0; // lagrangian
 real c2 = 0.0; // smooth vel
 real c3 = 0.1; // smooth - fujiwara
 real c4 = 0.1; // smooth surface - fujiwara
 real alpha = 1;
 real beta = 1;

 real sigma = 0.1;

 real mu_in = 1;
 real mu_out = 10;

 real rho_in = 1;
 real rho_out = 100;

 real cfl = 0.3;

 //string meshFile = "micro.msh";
 string meshFile = "2micros.msh";
 //string meshFile = "4micros.msh";

 Solver *solverP = new PetscSolver(KSPBICG,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 if( *(argv+1) == NULL )     
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  const char *mesh1 = mesh;
  m1.readMSH(mesh1);
  m1.setInterfaceBC();
  m1.setTriEdge(triEdge);
  m1.checkTriangleOrientation();
  m1.mesh2Dto3D();
  m1.setMiniElement();
  //m1.setQuadElement();
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setMicroWallBC();

  s1(m1);

  s1.setRe(Re);
  s1.setSc(Sc);
  s1.setWe(We);
  s1.setFr(Fr);
  s1.setC1(c1);
  s1.setC2(c2);
  s1.setC3(c3);
  s1.setC4(c4);
  s1.setAlpha(alpha);
  s1.setBeta(beta);
  s1.setSigma(sigma);
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.setCfl(cfl);
  s1.init();
  s1.setDt();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
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
  m1.setTriEdge(triEdge);
  m1.mesh2Dto3D();

  s1(m1);

  // load 3D mesh
  file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
  m1.setMiniElement();
  //m1.setQuadElement();
  m1.readVTKHeaviside(vtkFile);
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setMicroWallBC();

  s1(m1);

  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  iter = s1.loadSolution("sim",atoi(*(argv+2)));
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
  mOld.setOFace();

  // load surface mesh and create new mesh
  file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge(triEdge);
  m1.mesh2Dto3DOriginal();
  m1.setMiniElement();
  //m1.setQuadElement();
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setMicroWallBC();

  s1(m1);

  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
  iter = s1.loadSolution("sim",atoi(*(argv+2)));
  s1.applyLinearInterpolation(mOld);
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
  mOld.setOFace();

  // load surface mesh and create new one
  file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge(triEdge);
  m1.mesh2Dto3DOriginal();
  m1.setMiniElement();
  //m1.setQuadElement();
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();

  s1(m1);
  //file = (string) "sim-" + *(argv+2);
  //const char *sol = file.c_str();
  iter = s1.loadSolution("sim",atoi(*(argv+2)));
  s1.applyLinearInterpolation(mOld);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveMSH(mshFolder,"newMesh",atoi(*(argv+2)));
  saveEnd.saveSol(binFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveVTKTest(vtkFolder,"simCutPlane",atoi(*(argv+2)));
  //saveEnd.saveVTKSurface(vtkFolder,"sim",atoi(*(argv+2)));
  return 0;
 }
 // Point's distribution
 Laplace3D d1(m1);
 d1.init();
 d1.assemble();
 d1.setBC();
 d1.matMountC();
 d1.setUnCoupledCBC(); 
 d1.setCRHS();
 d1.unCoupledC();
 //d1.saveVTK("./vtk/","edge");
 d1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);
 save.printInfo(meshFile.c_str());

 int nIter = 3000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << i*nReMesh+j+iter << endl << endl;
   cout << resetColor();

   s1.setDt();
   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravity("Y");
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveMSH(mshFolder,"newMesh",i*nReMesh+j+iter);
   save.saveVTK(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveVTKTest(vtkFolder,"simCutPlane",i*nReMesh+j+iter);
   save.saveVTKSurface(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveSol(binFolder,"sim",i*nReMesh+j+iter);
   save.saveBubbleInfo(datFolder);
   save.printSimulationReport();
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",i*nReMesh+j+iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << i*nReMesh+j+iter << endl << endl;;
   cout << resetColor();
  }
  Laplace3D d2(m1,d1);
  d2.assemble();
  d2.setBC();
  d2.matMountC();
  d2.setUnCoupledCBC(); 
  d2.setCRHS();
  d2.unCoupledC();
  d2.saveVTK("./vtk/","edge",nReMesh+i*nReMesh+iter-1);
  d2.setModel3DEdgeSize();

  Model3D mOld = m1; 
  m1.setTriEdge(triEdge);

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();

  // 3D operations
  m1.insert3dMeshPointsByDiffusion(2.0);
  m1.remove3dMeshPointsByDiffusion(0.33);
  m1.removePointByVolume(0.005);
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.delete3DPoints();

  // surface operations
  m1.insertPointsByLength();
  //m1.insertPointsByCurvature();
  m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance();
  m1.contractEdgeByLength();
  //m1.removePointsByLength();
  m1.flipTriangleEdge();
  m1.checkNeighbours();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
  m1.setMiniElement();
  //m1.setQuadElement();
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setMicroWallBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveMSH(mshFolder,"newMesh",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveVTK(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveVTKSurface(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveVTKTest(vtkFolder,"simCutPlane",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveSol(binFolder,"sim",nReMesh+i*nReMesh+iter-1);
  //saveEnd.saveVTU(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
  //saveEnd.saveSolTXT(binFolder,"sim",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveMeshInfo(datFolder);
  saveEnd.printMeshReport();
 }

 PetscFinalize();
 return 0;
}



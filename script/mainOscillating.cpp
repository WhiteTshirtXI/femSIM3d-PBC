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
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 
 // set each bubble length
 vector< real > triEdge;
 triEdge.resize(3);
 triEdge[0] = 0.1; // none
 triEdge[1] = 0.77; // wall
 triEdge[2] = 0.09; // bubble 1 

 // Tryggvason (Computations of Multiphase Flows by a FDM/FTM
 int iter = 0;
 real Re = 1000;
 real Sc = 1;
 real We = 1;
 real Fr = 1;
 real c1 = 1.0;  // lagrangian
 real c2 = 0.0; // velocity
 real c3 = 0.1; // coordinates - fujiwara
 real c4 = 0.0; // surface
 real alpha = 1;
 real beta = 1;

 real sigma = 1;

 real mu_in = 1;
 real mu_out = 0.001;

 real rho_in = 1; 
 real rho_out = 0.001;

 real cfl = 0.5;

 //string meshFile = "oscillating.msh";
 //string meshFile = "static.msh";
 string meshFile = "newMesh-last.msh";

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
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
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setWallBC();
  m1.setSphereToEllipsoid(1.01);

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
  //s1.setDt(dt);
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
  m1.readVTKHeaviside(vtkFile);
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setWallBC();

  s1(m1);

  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  iter = s1.loadSolution("sim",atoi(*(argv+2)));
 }
 else
  cout << "The options are: NULL or restart" << endl;


 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);
 save.printInfo(meshFile.c_str());

 int nIter = 13000;
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
   //s1.setGravity("Z");
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
  Model3D mOld = m1; 
  m1.setTriEdge(triEdge);

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  //m1.setNormalAndKappa();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion(2.0);
  //m1.remove3dMeshPointsByDiffusion(0.33);
  //m1.removePointByVolume(0.005);
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  //m1.delete3DPoints();

  // surface operations
  //m1.insertPointsByLength();
  //m1.insertPointsByCurvature();
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance();
  //m1.contractEdgeByLength();
  //m1.removePointsByLength();
  //m1.flipTriangleEdge();
  //m1.checkNeighbours();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
  m1.setMiniElement();
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallBC();

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



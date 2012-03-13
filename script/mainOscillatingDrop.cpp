// =================================================================== //
// this is file mainOscillatingDrop.cpp, created at 10-Jan-2011        //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "TElement.h"
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
 triEdge.resize(2);
 triEdge[0] = 0.9; // wall
 triEdge[1] = 0.2; // bubble 1 

 int iter = 1;
 real Re = 1000;
 real Sc = 1;
 real We = 1;
 real Fr = 1;
 real c1 = 0.5;  // lagrangian
 real c2 = 1.0;  // smooth vel
 real c3 = 0.1;  // smooth coord (fujiwara)
 real d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 real d2 = 0.01;  // surface smooth cord (fujiwara)
 real alpha = 1;
 real beta = 1;

 real sigma = 1;

 real mu_in = 1;
 real mu_out = 0.01;

 real rho_in = 1; 
 real rho_out = 0.001;

 real cfl = 0.5;

 //string meshFile = "sphere.msh";
 string meshFile = "staticDumped5.msh";

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/sphere/" + meshFile;
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
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
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
  s1.setD1(d1);
  s1.setD2(d2);
  s1.setAlpha(alpha);
  s1.setBeta(beta);
  s1.setSigma(sigma);
  //s1.setDtALETwoPhase(dt);
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.setCfl(cfl);
  s1.init();
  s1.setDtALETwoPhase();
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
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.readVTKHeaviside(vtkFile);
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
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
	    << iter << endl << endl;
   cout << resetColor();

   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.setDtALETwoPhase();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravity("Z");
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKQuarter(vtkFolder,"simCutPlane",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   save.printSimulationReport();
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Model3D mOld = m1; 
  m1.setTriEdge(triEdge);

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  //m1.setNormalAndKappa();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  //m1.remove3dMeshPointsByDiffusion();
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
  //m1.removePointByNeighbourCheck();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
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
  saveEnd.saveMSH(mshFolder,"newMesh",iter-1);
  saveEnd.saveVTK(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKSurface(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKQuarter(vtkFolder,"simCutPlane",iter-1);
  saveEnd.saveSol(binFolder,"sim",iter-1);
  //saveEnd.saveVTU(vtkFolder,"sim",iter-1);
  //saveEnd.saveSolTXT(binFolder,"sim",iter-1);
  saveEnd.saveMeshInfo(datFolder);
  saveEnd.printMeshReport();
 }

 PetscFinalize();
 return 0;
}



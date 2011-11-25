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
#include "TElement.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 real Re = 25;
 real Sc = 2;
 real We = 2;
 real Fr = 0.4;
 real c1 = 0.0; // lagrangian
 real c2 = 0.05; // smooth
 real c3 = 0.0;
 real c4 = 0.03; // surface
 real alpha = 1;
 real beta = 1;

 real sigma = 1.0;

 real mu_in = 1.7894E-05;
 real mu_out = 0.001;

 real rho_in = 1.225;
 real rho_out = 1000;

 real cfl = 40;

 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU); 
 //Solver *solverP = new PetscSolver(KSPPREONLY,PCLU); // MUMPS
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 //string meshFile = "squareAnnular.msh";
 string meshFile = "cylinderAnnular.msh";

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1,mOld;
 Simulator3D s1,s2;

 if( *(argv+1) == NULL )     
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  const char *mesh1 = mesh;
  m1.readMSH(mesh1);
  m1.setInterfaceBC();
  m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallAnnularBC();

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
  //s1.setSigma(sigma);
  s1.setMu(mu_l,mu_g);
  s1.setRho(rho_l,rho_g);
  s1.setCfl(cfl);
  s1.init();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
 }
 else if( strcmp( *(argv+1),"restart") == 0 ) 
 {
  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;

  string aux = *(argv+2);
  string file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);

  m1.setInterfaceBC();
  m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallAnnularBC();

  s1(m1);

  s1.setRe(Re);
  s1.setSc(Sc);
  s1.setWe(We);
  s1.setFr(Fr);
  s1.setAlpha(alpha);
  s1.setBeta(beta);
  //s1.setSigma(sigma);
  s1.setCfl(cfl);
  s1.setMu(mu_l,mu_g);
  s1.setRho(rho_l,rho_g);
  s1.init();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

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
  m1.setWallAnnularBC();

  s2(m1,s1);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  file = (string) "sim-" + *(argv+2);
  const char *sol = file.c_str();
  s1.loadSolution(binFolder,sol);
  s1.setCfl(cfl);
  iter = s1.loadIteration(vtkFolder,sol);
 }
 else if( strcmp( *(argv+1),"remesh") == 0 )  
 {
  cout << endl;
  cout << "--------------> RE-MESHING (NO ITERATION)..." << endl;
  cout << endl;

  string aux = *(argv+2);
  string file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallAnnularBC();

  s1(m1);

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
  m1.setWallAnnularBC();

  s2(m1,s1);
  s1 = s2;

  file = (string) "sim-" + *(argv+2);
  const char *sol = file.c_str();
  s1.loadSolution(binFolder,sol);
  s1.setCfl(cfl);
  iter = s1.loadIteration(vtkFolder,sol);

  mOld = m1; 
  m1.mesh2Dto3DOriginal();
  m1.setInterfaceDistance();

  s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveMSH(mshFolder,"newMesh",atoi(*(argv+2)));
  saveEnd.saveSol(binFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveVTKTest(vtkFolder,"simCutPlane",atoi(*(argv+2)));
  //saveEnd.saveVTKSurface(vtkFolder,"sim",atoi(*(argv+2)));
  return 0;
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo("./","info",mesh);
 save.printInfo(meshFile.c_str());

 int nIter = 1;
 int nReMesh = 3;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << resetColor() << iter << endl;

   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.setDtALETwoPhase();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   s1.setGravity("Z");
   //s1.setGravityBoussinesq("Z");
   s1.setInterfaceGeo();
   //s1.setInterfaceGeoTest();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTKTest(vtkFolder,"simCutPlane",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << resetColor() << iter << endl;

   iter++;
  }
  mOld = m1; 

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
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallAnnularBC();

  s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim",iter-1);
  saveEnd.saveVTU(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKSurface(vtkFolder,"sim",iter-1);
  saveEnd.saveMSH(mshFolder,"newMesh",iter-1);
  saveEnd.saveSol(binFolder,"sim",iter-1);
  saveEnd.saveSimTime(iter-1);
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}



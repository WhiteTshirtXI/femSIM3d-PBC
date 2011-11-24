// =================================================================== //
// this is file mainCurvature.cpp, created at 30-Sep-2011              //
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
 triEdge[2] = 0.1; // bubble 1 

 // Tryggvason (Computations of Multiphase Flows by a FDM/FTM
 real Re = 1000;
 real Sc = 1;
 real We = 1;
 real Fr = 1;
 real c1 = 0.0;  // lagrangian
 real c2 = 0.0; // velocity
 real c3 = 0.2; // coordinates - fujiwara
 real c4 = 0.0; // surface
 real alpha = 1;
 real beta = 1;

 real sigma = 1.0;

 real mu_in = 1;
 real mu_out = 0.001;

 real rho_in = 1; 
 real rho_out = 0.01;

 real cfl = 0.8;

 string meshFile = "static.msh";

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

 const char *mesh1 = mesh;
 m1.readMSH(mesh1);
 m1.setInterfaceBC();
 m1.setTriEdge(triEdge);
 m1.checkTriangleOrientation();
 m1.mesh2Dto3D();
 m1.setMiniElement();
 m1.setOFace();
 m1.setSurfaceConfig();
 //m1.setBiggerSphere(1);
 m1.setInitSurfaceVolume();
 m1.setWallBC();

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
 //s1.setDtALETwoPhase(dt);
 s1.setMu(mu_in,mu_out);
 s1.setRho(rho_in,rho_out);
 s1.setCfl(cfl);
 s1.init();
 s1.setDtALETwoPhase();
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);
 save.printInfo(meshFile.c_str());

 int nIter = 1;
 for( int i=1;i<=nIter;i++ )
 {
  cout << color(none,magenta,black);
  cout << "____________________________________ Iteration: " 
       << i << endl << endl;
  cout << resetColor();

  //s1.setDtALETwoPhase();
  //s1.stepLagrangian();
  //s1.stepALE();
  s1.stepALEVel();
  s1.assemble();
  s1.matMount();
  s1.setUnCoupledBC();
  s1.setRHS();
  //s1.setInterface();
  s1.setInterfaceGeo();
  s1.unCoupled();

  InOut save(m1,s1); // cria objeto de gravacao
  save.saveMSH(mshFolder,"newMesh",i);
  save.saveVTK(vtkFolder,"sim",i);
  save.saveVTKTest(vtkFolder,"simCutPlane",i);
  save.saveVTKSurface(vtkFolder,"sim",i);
  save.saveSol(binFolder,"sim",i);
  save.saveBubbleInfo(datFolder);

  s1.saveOldData();

  cout << color(none,magenta,black);
  cout << "________________________________________ END of " 
       << i << endl << endl;;
  cout << resetColor();
 }

 PetscFinalize();
 return 0;
}


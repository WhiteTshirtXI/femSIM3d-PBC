// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "Simulator3D.h"
#include "Interface3D.h"
#include "InOut.h"
#include "Mumps_Petsc.h"
#include "PetscSolver.h"
#include "petscksp.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 const char *dir  = "./";
 const char *mesh = "../../db/mesh/3d/cube-cube2D-1.vtk";
 const char *txt  = "txt/txt";
 const char *bin  = "bin/bin";
 const char *vtk  = "vtk/sim";
 const char *sim  = "dat/sim";
 int iter = 0;

 Model3D m1,mNew,mTemp;
 m1.readVTKSurface(mesh);
 m1.meshAll();
 m1.setMiniElement();
 m1.setCubeCubeBC();
 m1.setCube(1.1,1.9,1E-10);
 m1.setOFace();
 m1.setSurfaceTri();

 Simulator3D s1(m1);

 s1.setRe(400);
 s1.setSc(2);
 s1.setWe(10);
 s1.setAlpha(1);
 s1.setBeta(0.0);
 //s1.setSigma(0.0);
 s1.setCflBubble(20);
 s1.init();

 s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));

 const int restart = 0;

 if( restart == 1 )
 {
  const char *mesh2 = "./vtk/sim-last-0.vtk";

  m1.readVTK(mesh2);
  m1.readVTKCC(mesh2);
  m1.meshRestart();
  m1.setMiniElement2();
  m1.setCubeCubeBC();
  m1.setOFace();
  m1.setSurfaceTri();

  //Simulator3D s2(m1,s1);
  Simulator3D s2(m1);

  s1 = s2; 
  s1.setRe(400);
  s1.setSc(2);
  s1.setWe(10);
  s1.setAlpha(1);
  s1.setBeta(0.0);
  s1.setCflBubble(20);
  s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
  s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
  s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));
  s1.loadSolution(dir,"./vtk/sim-last",0); // set para velocidade no simulador
  iter = s1.loadIteration();
 }

 {
  InOut save(m1,s1); // cria objeto de gravacao
  save.saveVTK(dir,vtk);
  save.saveVTKTri("./vtk/","sim-lastRestart",0);
  save.saveInfo(dir,mesh);
  save.printInfo(dir,mesh);
 }
 
 int nIter = 100;
 int nReMesh = 5;
 for( int i=0;i<nIter;i++ )
 {
  InOut save(*s1.m,s1); // cria objeto de gravacao

  for( int j=0;j<nReMesh;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << i*nReMesh+j+iter << endl;
   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALE2();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravityBoussinesq();
   //s1.setInterface();
   s1.setInterfaceGeo();
   //s1.setInterfaceGeoTest();
   s1.unCoupled();

   save.saveVTK(dir,vtk,i*nReMesh+j+iter);
   save.saveVTKTri(dir,vtk,i*nReMesh+j+iter);
   save.oscillating("oscillating.dat");
   save.oscillatingD("oscillatingD.dat");
   save.oscillatingKappa("oscillatingKappa.dat");
   cout << "________________________time: " << s1.getTime2() << endl;
  }
  mNew = *s1.m; // operador (=) nao esta redefinido em Model3D!!!
  mNew.reMeshAll2();
  mNew.setMiniElement2();
  mNew.setCubeCubeBC();
  mNew.setOFace();
  mNew.setSurfaceTri();

  Simulator3D s2(mNew,s1);
  s2.applyLinearInterpolation(*s1.m);
  s1 = s2; 
  s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
  s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
  s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));

 }

 InOut save(*s1.m,s1); // cria objeto de gravacao
 save.saveVTK("./vtk/","sim-last",0);
 save.saveVTKTri("./vtk/","sim-last",0);
 save.saveSol("./vtk/","sim-last",0);
 save.saveSolTXT("./vtk/","sim-last",0);
 save.saveSimTime(iter+nReMesh*nIter);


 PetscFinalize();
 return 0;
}



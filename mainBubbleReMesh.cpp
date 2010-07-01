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
 //const char *mesh = "../../db/mesh/3d/cube-cube1.vtk";
 const char *mesh = "../../db/mesh/3d/bubble-bubble1.vtk";
 //const char *mesh = "../../db/mesh/3d/bubble4-9-20.vtk";
 //const char *mesh = "../../db/mesh/3d/bubble8-31-2.vtk";
 const char *txt  = "txt/txt";
 const char *bin  = "bin/bin";
 const char *vtk  = "vtk/sim";
 const char *sim  = "dat/sim";

 Model3D m1;
 m1.readVTK(mesh);
 m1.setMiniElement();
 //m1.setCubeCubeBC(1.5);
 m1.setBubbleBubbleBC();
 //m1.setBubbleBC2(); // malha do disco
 //m1.setBubble3DBC();
 m1.setOFace();
 m1.reMeshAll();
 
 Simulator3D s1(m1);

 s1.setRe(100);
 s1.setSc(2);
 s1.setWe(10);
 s1.setAlpha(1);
 s1.setBeta(-2.0);
 //s1.setDt(0.1);
 s1.setCflBubble(10);
 s1.init();

 s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 //s1.setSolverPressure(new PetscSolver(KSPGMRES,PCJACOBI));
 s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));
 
 InOut save(m1,s1); // cria objeto de gravacao

 save.saveVTK(dir,vtk);
 save.saveInfo(dir,mesh);
 save.printInfo(dir,mesh);
 save.saveVTKTri(dir,vtk,0);

 int nReMesh = 2;
 for( int i=0;i<100;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << i*nReMesh+j << endl;
   s1.setCentroid();
   s1.stepLagrangian();
   //s1.stepALE();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravityBoussinesq();
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.unCoupled();
   save.saveVTK(dir,vtk,i*nReMesh+j);
   save.saveVTKTri(dir,vtk,i*nReMesh+j);
   save.oscillating("oscillating.dat");
   save.oscillatingD("oscillatingD.dat");
   save.oscillatingKappa("oscillatingKappa.dat");
   //save.oscillating(0,1,4,"oscillating.dat");
   //save.oscillatingD(0,2,1,3,4,5,"oscillatingD.dat");
  }
 m1.reMeshAll();
 Simulator3D s2(m1,s1);
 s1 = s2; 
 s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));
 }

 PetscFinalize();
 return 0;
}



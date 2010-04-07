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
 const char *mesh = "../../db/mesh/3d/bubble-bubble2.vtk";
 //const char *mesh = "../../db/mesh/3d/bubble-cube2.vtk";
 //const char *mesh = "../../db/mesh/3d/bubble4-9-20.vtk";
 //const char *mesh = "../../db/mesh/3d/bubble8-31-2.vtk";
 const char *txt  = "txt/txt";
 const char *bin  = "bin/bin";
 const char *vtk  = "vtk/sim";
 const char *sim  = "dat/sim";

 Model3D m1;
 m1.readVTK(mesh);
 m1.setCentroid();
 //m1.setCubeCubeBC(1.5);
 m1.setBubbleBubbleBC();
 //m1.setBubbleBC2(); // malha do disco
 //m1.setBubble3DBC();
 m1.setOFace();

 Simulator3D s1(m1);

 s1.setRe(100);
 s1.setSc(2);
 s1.setWe(10);
 s1.setAlpha(1);
 s1.setBeta(0);
 //s1.setDt(0.1);
 s1.setCflBubble(1);
 s1.init();

 /* ASSEMBLE */
 //s1.setSolverPressure(new PetscSolver(KSPCG,PCJACOBI));
 //s1.setSolverPressure(new PetscSolver(KSPCG,PCICC));
 //s1.setSolverPressure(new PetscSolver(KSPGMRES,PCNONE));
 //s1.setSolverPressure(new PetscSolver(KSPGMRES,PCJACOBI));
 //s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 //s1.setSolverPressure(new PetscSolver(KSPBICG,PCJACOBI));
 //s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 //s1.setSolverPressure(new GMRes);
 //s1.setSolverVelocity(new PCGSolver);
 //s1.setSolverConcentration(new PCGSolver);
 //s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));
 
 /* ASSEMBLESLIP */
 s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(dir,vtk);
 save.saveInfo(dir,mesh);
 save.printInfo(dir,mesh);

 for( int i=0;i<1000;i++ )
 {

  for( int j=0;j<10;j++ )
  {
   cout << "____________________________________ Iteration: " << i*10+j << endl;
   //s1.stepLagrangian();
   s1.stepALE();
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();
   s1.setUnCoupledCBC();
   s1.setRHS();
   //s1.setGravity();
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.unCoupled();
   save.saveVTK(dir,vtk,i*10+j);
   save.oscillating("oscillating.dat");
   save.oscillatingD("oscillatingD.dat");
   save.oscillatingKappa("oscillatingKappa.dat");
   //save.oscillating(0,1,4,"oscillating.dat");
   //save.oscillatingD(0,2,1,3,4,5,"oscillatingD.dat");
  }
 }

 PetscFinalize();
 return 0;
}



// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
//#include "GMRes.h"
#include "Simulator3D.h"
#include "InOut.h"
#include "PetscSolver.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

 const char *dir  = "./";
 //const char *mesh = "../../db/mesh/3d/disk8-20-6.vtk";
 const char *mesh = "../../db/mesh/3d/disk6-6-6.vtk";
 const char *txt  = "txt/txt";
 const char *bin  = "bin/bin";
 const char *vtk  = "vtk/sim";

 Model3D m1;
 m1.setMeshDisk(6,6,6);
 m1.setMiniElement();
 m1.setDiskFSBC();
 m1.setPerturbSurf();
 //m1.setPerturbSurf2();
 //m1.setPerturbSurfSquare();
 m1.setOFace();

 Simulator3D s1(m1);

 s1.init();

 s1.setRe(1000);
 s1.setSc(1000);
 s1.setFr(10);
 s1.setBeta(1);
 //s1.setDt(dt);
 s1.setCfl(1);
 s1.setSolverVelocity(new PCGSolver());
 //s1.setSolverPressure(new PCGSolver());
 //s1.setSolverPressure(new GMRes());
 s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 s1.setSolverConcentration(new PCGSolver());
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(dir,vtk);
 save.saveInfo(dir,mesh);
 save.printInfo(dir,mesh);

 for( int i=0;i<1000;i++ )
 {
  cout << "____________________________________ Iteration: " << i << endl;
  //s1.stepLagrangian();
  s1.stepALE();
  s1.matMount();
  s1.matMountC();
  s1.setUnCoupledBC();
  s1.setUnCoupledCBC();
  s1.setRHS();
  s1.setCRHS();
  s1.setGravity();
  s1.unCoupled();
  s1.unCoupledC();
  save.saveVTK(dir,vtk,i);
 };

 PetscFinalize();
 return 0;
}



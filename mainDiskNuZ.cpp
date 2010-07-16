// =================================================================== //
// this is file mainDiskNuCte.cpp, created at 10-Jun-2007              //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "InOut.h"
#include "Simulator3D.h"
#include "PetscSolver.h"
#include "petscpc.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

 const char *mesh = "../../db/mesh/3d/disk6-10-20.vtk";
 const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
 const char *simFolder  = "./sim/";
 const char *pertFolder  = "./pert/";


 Model3D m1;
 m1.setMeshDisk(6,10,20);
 m1.setAdimenDisk();
 m1.setMiniElement();
 m1.setNuZDiskBC();;
 m1.setOFace();
 //m1.readBaseStateNu("NuCte");

 Simulator3D s1(m1);

 s1.setRe(1);
 s1.setCflDisk(10);
 //s1.setSolverVelocity(new PetscSolver(KSPCG,PCJACOBI));
 //s1.setSolverPressure(new PetscSolver(KSPGMRES,PCJACOBI));
 s1.setSolverVelocity(new PCGSolver());
 s1.setSolverPressure(new PCGSolver());

 s1.init();
 s1.assembleNuZ();
 s1.matMount();
 s1.setUnCoupledBC(); 

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./",mesh);
 save.printInfo("./",mesh);
 
 for( int i=0;i<1000;i++ )
 {
  for( int j=0;j<10;j++ )
  {
   cout << "____________________________________ Iteration: " << i*10+j << endl;
   s1.stepSL();
   s1.setRHS();
   s1.unCoupled();
   save.saveVonKarman(simFolder,"vk1",i*10+j,4);
   save.saveVonKarman(simFolder,"vk2",i*10+j,5);
   save.saveVonKarman(simFolder,"vk3",i*10+j,6);
   save.saveVonKarman(simFolder,"vk4",i*10+j,7);
   save.saveVonKarman(simFolder,"vk5",i*10+j,8);
   save.saveVonKarman(simFolder,"vk6",i*10+j,9);
   save.saveVonKarman(simFolder,"vk7",i*10+j,10);
   //save.savePert(pertFolder,"pert",i*10+j,7);
   s1.convergenceCriteria(10E-6);
  }
  save.saveVTK(vtkFolder,"sim",i);
  save.saveSol(binFolder,"UVWPC",i);
 }

 PetscFinalize();
 return 0;
}



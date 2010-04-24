// =================================================================== //
// this is file mainDiskNuCte.cpp, created at 10-Jun-2007              //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "PCGSolver.h"
//#include "petscpc.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
 const char *dir  = "./";
 const char *mesh = "../../db/mesh/3d/disk6-10-20-var-.vtk";
 const char *vtk  = "vtk/sim";
 const char *txt  = "txt/txt";
 const char *bin  = "bin/bin";
 const char *von1  = "sim/vk1";
 const char *von2  = "sim/vk2";
 const char *von3  = "sim/vk3";
 const char *von4  = "sim/vk4";
 const char *von5  = "sim/vk5";
 const char *von6  = "sim/vk6";
 const char *von7  = "sim/vk7";
 const char *vonPert = "pert/sim1";


 Model3D m1;
 m1.readVTK(mesh);
 m1.setNuCDiskBC();
 m1.setAdimenDisk();
 m1.setMiniElement();
 m1.setOFace();

 Simulator3D s1(m1);
 s1.init();
 s1.assembleSlip();

 s1.setRe(1); // Reynolds do disco (~1)
 s1.setSc(2000); // Schmidt da concentracao (~2000)
 s1.setCflDisk(1);
 s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 s1.setSolverPressure(new PetscSolver(KSPBICG,PCJACOBI));
 s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(dir,vtk);
 save.saveInfo(dir,mesh);
 save.printInfo(dir,mesh);

 //save.loadSol(s1,dir,bin,count);
 for( int i=0;i<1000;i++ )
 {
  for( int j=0;j<10;j++ )
  {
   cout << "____________________________________ Iteration: " << i*10+j << endl;
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();
   s1.setUnCoupledCBC();
   s1.stepSL();
   s1.setRHS();
   s1.setCRHS();
   s1.unCoupled();
   s1.unCoupledC();
   s1.assembleK();
   save.saveVonKarman(dir,von1,i*10+j,4);
   save.saveVonKarman(dir,von2,i*10+j,5);
   save.saveVonKarman(dir,von3,i*10+j,6);
   save.saveVonKarman(dir,von4,i*10+j,7);
   save.saveVonKarman(dir,von5,i*10+j,8);
   save.saveVonKarman(dir,von6,i*10+j,9);
   save.saveVonKarman(dir,von7,i*10+j,10);
   s1.convergenceCriteria(10E-06);
   save.saveVTK(dir,vtk,i*10+j);
  }
  //save.saveSol(dir,bin,i);
 }

 PetscFinalize();
 return 0;
}



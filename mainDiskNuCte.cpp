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
#include "Simulator3D.h"
#include "InOut.h"
#include "PCGSolver_Petsc.h"
#include "PetscSolver.h"
#include "GMRes_Petsc.h"
#include "petscpc.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

 int iter = 0;
 real Re = 1;
 real cfl = 10;
 Solver *solverP = new PetscSolver(KSPBICG,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);

 const char *mesh = "../../db/mesh/3d/disk6-10-20.vtk";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *simFolder  = "./sim/";

 Model3D m1;
 m1.setMeshDisk(6,12,40);
 //m1.readVTK(mesh);
 m1.setAdimenDisk();
 m1.setMiniElement();
 m1.setNuCteDiskBC();
 //m1.readAndSetPressureDiskBC("../../db/baseState/nuCte/","p");
 //m1.setCDiskBC();
 m1.setOFace();
 //m1.readBaseStateNu("NuCte");

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setCflDisk(cfl);
 s1.setSolverVelocity(solverV);
 s1.setSolverPressure(solverP);

 s1.init();
 s1.assembleNuCte();
 //s1.assembleC();
 s1.matMount();
 //s1.matMountC();
 s1.setUnCoupledBC(); 
 //s1.setUnCoupledCBC(); 
 
 if( (*(argv+1)) == NULL )
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;
 }
 else if( strcmp( *(argv+1),"restart") == 0 )
 {
  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;
  s1.loadSolution(binFolder,"sim-last");
  iter = s1.loadIteration(vtkFolder,"sim-last");
  //s1.loadSolution(binFolder,"UVWPC",70);
  //iter = s1.loadIteration(vtkFolder,"sim",70);
 }
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);
 save.printInfo(mesh);
 
 int nIter = 1000;
 int nR = 10;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nR;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << i*nR+j+iter << endl;

   s1.stepSL();
   s1.setRHS();
   //s1.setCRHS();
   s1.unCoupled();
   //s1.unCoupledC();
   save.saveVonKarman(simFolder,"vk1",i*nR+j+iter,4);
   save.saveVonKarman(simFolder,"vk2",i*nR+j+iter,5);
   save.saveVonKarman(simFolder,"vk3",i*nR+j+iter,6);
   save.saveVonKarman(simFolder,"vk4",i*nR+j+iter,7);
   save.saveVonKarman(simFolder,"vk5",i*nR+j+iter,8);
   save.saveVonKarman(simFolder,"vk6",i*nR+j+iter,9);
   save.saveVonKarman(simFolder,"vk7",i*nR+j+iter,10);
   save.saveVTK(vtkFolder,"sim",i*nR+j+iter);
   save.saveVTKTest(vtkFolder,"simTest",i*nR+j+iter);
   save.saveVTU(vtkFolder,"sim",i*nR+j+iter);
   save.saveSol(binFolder,"UVWPC",i*nR+j+iter);
   save.saveConvergence("./","convergence");

   cout << "__________________________________________ End: " 
	    << i*nR+j+iter << endl;
  }
 }

 PetscFinalize();
 return 0;
}



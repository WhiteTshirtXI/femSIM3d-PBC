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
#include "PetscSolver.h"
#include "petscpc.h"
#include "colors.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

 int iter = 0;
 real Re = 1;
 real cfl = 0.1;
 real mu_l = 1.0;
 real rho_l = 1.0;
 //Solver *solverP = new PetscSolver(KSPPREONLY,PCLU);
 Solver *solverP = new PetscSolver(KSPBICG,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);

 const char *mesh = "../../db/mesh/3d/disk6-10-20.vtk";
 const char *binFolder  = "./bin/";
 const char *datFolder  = "./dat/";
 const char *vtkFolder  = "./vtk/";
 const char *simFolder  = "./sim/";

 Model3D m1;
 m1.setMeshDisk(6,6,40);
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
 s1.setMu(mu_l);
 s1.setRho(rho_l);
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

  iter = s1.loadSolution("sim",atoi(*(argv+2)));
  s1.setCflDisk(cfl);
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
   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: "
	    << i*nR+j+iter << endl << endl;
   cout << resetColor();

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
   save.saveVonKarman(simFolder,"vk7",i*nR+j+iter,4);
   save.saveVonKarman(simFolder,"vk8",i*nR+j+iter,5);
   save.saveVonKarman(simFolder,"vk9",i*nR+j+iter,6);
   save.saveVonKarman(simFolder,"vk10",i*nR+j+iter,7);
   save.saveVonKarman(simFolder,"vk11",i*nR+j+iter,8);
   save.saveVonKarman(simFolder,"vk12",i*nR+j+iter,10);
   save.saveVTK(vtkFolder,"sim",i*nR+j+iter);
   save.saveSol(binFolder,"sim",i*nR+j+iter);
   save.saveConvergence(datFolder,"convergence");

   cout << color(none,magenta,black);
   cout << "________________________________________ END of "
	    << i*nR+j+iter << endl << endl;;
   cout << resetColor();

  }
 }

 PetscFinalize();
 return 0;
}



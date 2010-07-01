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
#include "InOut.h"

int main(int argc, char **argv)
{
 const char *dir  = "./";
 //const char *mesh = "../../db/mesh/3d/step2-2-2.vtk";
 //const char *mesh = "../../db/mesh/3d/step5-5-2.vtk";
 //const char *mesh = "../../db/mesh/3d/step20-10-2.vtk";
 const char *mesh = "../../db/mesh/3d/step40-20-2.vtk";
 const char *txt  = "txt/txt";
 const char *bin  = "bin/bin";
 const char *vtk  = "vtk/sim";


 Model3D m1;
 //m1.readVTK(mesh);
 m1.setMeshStep(40,20,2);
 //m1.setStepReservInvBC();
 m1.setAdimenStep();
 m1.setMiniElement();
 //m1.setQuadElement();
 m1.setStepBC();
 m1.setOFace();

 Simulator3D s1(m1);

 s1.setRe(10000);
 s1.setSc(2000);
 s1.setFr(10);
 //s1.setDt(dt);
 s1.setCfl(3);
 s1.setSolverVelocity(new PCGSolver());
 s1.setSolverPressure(new PCGSolver());
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(dir,vtk);
 save.saveInfo(dir,mesh);
 save.printInfo(dir,mesh);

 //save.saveVTKFreeFace(m1,vtkDir);
 s1.init();
 s1.assembleSlip();
 s1.matMount();
 s1.matMountC();
 s1.setUnCoupledBC(); 
 s1.setUnCoupledCBC(); 

 //save.loadSol(s1,binDir,count);
 for( int i=0;i<100;i++ )
 {
  for( int j=0;j<10;j++ )
  {
   cout << "____________________________________ Iteration: " << i*10+j << endl;
   s1.stepSL();
   s1.setRHS();
   s1.setCRHS();
   s1.unCoupled();
   s1.unCoupledC();
   save.saveVTK(dir,vtk,i*10+j);
  }
  //save.saveSol(dir,bin,i);
  //save.saveSolTXT(dir,bin,i);
 }

 return 0;
}



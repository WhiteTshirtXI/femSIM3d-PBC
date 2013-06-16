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

#define NUMPHASES 1

int main(int argc, char **argv)
{
 const char *dir  = "./";

 string meshFile = "reserv1.vtk";

 const char *bcFolder = "../../db/mesh/3d/reserv1.bc";
 const char *txtFolder = "txt/txt";
 const char *binFolder = "bin/bin";
 const char *vtkFolder = "vtk/sim";
 const char *datFolder = "dat/sim";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/mesh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 m1.readVTK(mesh);
 m1.readBC(bc);
 m1.setMapping();
 m1.setMiniElement();

 Simulator3D s1(m1);

 s1.setRe(10000);
 //s1.setDt(dt);
 s1.setCfl(1);
 s1.setSolverVelocity(new PCGSolver());
 //s1.setSolverPressure(new PCGSolver());
 s1.setSolverPressure(new GMRes());
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"sim");
 save.saveInfo(datFolder,"info",mesh);

 s1.init();
 s1.assembleSlip();
 s1.matMount();
 //s1.matMountC();
 s1.setUnCoupledBC(); 
 //s1.setUnCoupledBC(); 
 
 //save.loadSol(binFolder,count);
 for( int i=0;i<100;i++ )
 {
  for( int j=0;j<10;j++ )
  {
   cout << "____________________________________ Iteration: " << i*10+j << endl;
   s1.stepSL();
   s1.setRHS();
   //s1.setCRHS();
   s1.unCoupled();
   //s1.unCoupledC();
   save.saveVTK(vtkFolder,"sim",i*10+j);

   s1.saveOldData();

   s1.timeStep();
  }
  //save.saveVTK(vtkFolder,"sim",i*10+j);
  //save.saveSol(binFolder,"bin",i*10+j);
  //save.saveSolTXT(txtFolder,"sol",i*10+j);
 }

 return 0;
}



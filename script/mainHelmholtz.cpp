// =================================================================== //
// this is file mainRisingBubble.cpp, created at 10-Jun-2009           //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Helmholtz3D.h"
#include "PetscSolver.h"
#include "petscksp.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 string meshFile = "2bubbles.msh";
 
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/2bubbles/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;

 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge();
 m1.mesh2Dto3D();
 m1.setMiniElement();
 m1.setOFace();
 m1.setSurfaceConfig();

 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initMicro();
 h1.assemble();
 h1.setk(3.8);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 h1.saveVTK(vtkFolder,"edge1",0);
 h1.saveChordalEdge(datFolder,"edge1",0);

 h1.setk(1.2);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 h1.saveVTK(vtkFolder,"edge2",0);
 h1.saveChordalEdge(datFolder,"edge2",0);

 h1.setk(0.03);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 h1.saveVTK(vtkFolder,"edge3",0);
 h1.saveChordalEdge(datFolder,"edge3",0);

 PetscFinalize();
 return 0;
}


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
 InOut save(m1,s1); // cria objeto de gravacao

 save.saveVTK(dir,vtk);
 save.saveVTKTri(dir,vtk,0);

 PetscFinalize();
 return 0;
}



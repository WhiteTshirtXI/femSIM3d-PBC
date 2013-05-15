// =================================================================== //
// this is file mainCurvature.cpp, created at 30-Sep-2011              //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include "Model3D.h"
#include "InOut.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 const char *mshFolder  = "./msh/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";

 /* mesh */
 const char* mesh = "../../db/gmsh/3d/sphere/kappa/initial.msh";


 Model3D m1;

 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge();

 m1.restoreMappingArrays();
 m1.setNormalAndKappa();
 m1.setSurfaceVolume();
 m1.setSurfaceArea();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();

 InOut save(m1); // cria objeto de gravacao
 save.saveMSH(mshFolder,"newMesh");
 save.saveVTKSurface(vtkFolder,"sim");
 save.saveKappaErrorSphere(datFolder);
 save.saveVolumeCorrection(datFolder);

 return 0;
}



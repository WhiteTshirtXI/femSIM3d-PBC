// =================================================================== //
// this is file mainDiskNuCte.cpp, created at 10-Jun-2007              //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "InOut.h"
#include "colors.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 const char *datFolder  = "./dat/";
 const char *mshFolder  = "./msh/";
 const char *vtkFolder  = "./vtk/";

 Model3D m1;
 //m1.setMeshTorus(1.5,0.5);

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
 save.saveKappaErrorTorus(datFolder);
 save.saveVolumeCorrection(datFolder);

 return 0;
}

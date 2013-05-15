// =================================================================== //
// this is file mainCurvature.cpp, created at 30-Sep-2011              //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "InOut.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 const char *mshFolder  = "./msh/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";

 /* meshes */
 vector<const char*> mesh;
 mesh.resize(43);
 mesh[0]  =  "../../db/gmsh/3d/torus/curvature/0.90.msh";
 mesh[1]  =  "../../db/gmsh/3d/torus/curvature/0.85.msh";
 mesh[2]  =  "../../db/gmsh/3d/torus/curvature/0.80.msh";
 mesh[3]  =  "../../db/gmsh/3d/torus/curvature/0.75.msh";
 mesh[4]  =  "../../db/gmsh/3d/torus/curvature/0.70.msh";
 mesh[5]  =  "../../db/gmsh/3d/torus/curvature/0.65.msh";
 mesh[6]  =  "../../db/gmsh/3d/torus/curvature/0.60.msh";
 mesh[7]  =  "../../db/gmsh/3d/torus/curvature/0.55.msh";
 mesh[8]  =  "../../db/gmsh/3d/torus/curvature/0.50.msh";
 mesh[9]  =  "../../db/gmsh/3d/torus/curvature/0.45.msh";
 mesh[10]  = "../../db/gmsh/3d/torus/curvature/0.40.msh";
 mesh[11]  = "../../db/gmsh/3d/torus/curvature/0.35.msh";
 mesh[12]  = "../../db/gmsh/3d/torus/curvature/0.30.msh";
 mesh[13]  = "../../db/gmsh/3d/torus/curvature/0.25.msh";
 mesh[14]  = "../../db/gmsh/3d/torus/curvature/0.24.msh";
 mesh[15]  = "../../db/gmsh/3d/torus/curvature/0.23.msh";
 mesh[16]  = "../../db/gmsh/3d/torus/curvature/0.22.msh";
 mesh[17]  = "../../db/gmsh/3d/torus/curvature/0.21.msh";
 mesh[18]  = "../../db/gmsh/3d/torus/curvature/0.20.msh";
 mesh[19]  = "../../db/gmsh/3d/torus/curvature/0.19.msh";
 mesh[20]  = "../../db/gmsh/3d/torus/curvature/0.18.msh";
 mesh[21]  = "../../db/gmsh/3d/torus/curvature/0.17.msh";
 mesh[22]  = "../../db/gmsh/3d/torus/curvature/0.16.msh";
 mesh[23] =  "../../db/gmsh/3d/torus/curvature/0.15.msh";
 mesh[24] =  "../../db/gmsh/3d/torus/curvature/0.14.msh";
 mesh[25] =  "../../db/gmsh/3d/torus/curvature/0.13.msh";
 mesh[26] =  "../../db/gmsh/3d/torus/curvature/0.12.msh";
 mesh[27] =  "../../db/gmsh/3d/torus/curvature/0.11.msh";
 mesh[28] =  "../../db/gmsh/3d/torus/curvature/0.10.msh";
 mesh[29] =  "../../db/gmsh/3d/torus/curvature/0.095.msh";
 mesh[30] =  "../../db/gmsh/3d/torus/curvature/0.09.msh";
 mesh[31] =  "../../db/gmsh/3d/torus/curvature/0.085.msh";
 mesh[32] =  "../../db/gmsh/3d/torus/curvature/0.08.msh";
 mesh[33] =  "../../db/gmsh/3d/torus/curvature/0.075.msh";
 mesh[34] =  "../../db/gmsh/3d/torus/curvature/0.07.msh";
 mesh[35] =  "../../db/gmsh/3d/torus/curvature/0.065.msh";
 mesh[36] =  "../../db/gmsh/3d/torus/curvature/0.06.msh";
 mesh[37] =  "../../db/gmsh/3d/torus/curvature/0.055.msh";
 mesh[38] =  "../../db/gmsh/3d/torus/curvature/0.05.msh";
 mesh[39] =  "../../db/gmsh/3d/torus/curvature/0.045.msh";
 mesh[40] =  "../../db/gmsh/3d/torus/curvature/0.04.msh";
 mesh[41] =  "../../db/gmsh/3d/torus/curvature/0.035.msh";
 mesh[42] =  "../../db/gmsh/3d/torus/curvature/0.03.msh";

 for( int i=0;i<(int) mesh.size();i++ )
 {
  cout << color(none,magenta,black);
  cout << "____________________________________ Iteration: " 
       << i << endl << endl;
  cout << resetColor();

  Model3D m1;
  Simulator3D s1(m1);

  m1.readMSH(mesh[i]);
  m1.setInterfaceBC();
  m1.setTriEdge();

  m1.restoreMappingArrays();
  m1.setNormalAndKappa();
  m1.setSurfaceVolume();
  m1.setSurfaceArea();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();

  InOut save(m1,s1); // cria objeto de gravacao
  save.saveMSH(mshFolder,"newMesh",i);
  save.saveVTKSurface(vtkFolder,"sim",i);
  save.saveKappaErrorTorus(datFolder);
  save.saveVolumeCorrection(datFolder);

  cout << color(none,magenta,black);
  cout << "________________________________________ END of " 
       << i << endl << endl;;
  cout << resetColor();
 }

 return 0;
}



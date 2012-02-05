// =================================================================== //
// this is file mainCurvature.cpp, created at 30-Sep-2011              //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "TElement.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";

 string meshFile = (string) getenv("DATA_DIR");
 meshFile += "/gmsh/3d/curvature/0.60.msh";
 const char *mesh = meshFile.c_str();

 // set each bubble length
 vector< real > triEdge;
 triEdge.resize(2);
 triEdge[0] = 0.6; // wall
 triEdge[1] = 0.60; // bubble 1 

 Model3D m1;
 Simulator3D s1;

 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge(triEdge);
 m1.checkTriangleOrientation();

 m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setOFace();
 m1.setSurfaceConfig();
 m1.setBiggerSphere(0.25);

 s1(m1);
 s1.setInterfaceGeo();

 // adding vertex
 Model3D mOld = m1; 
 m1.setTriEdge(triEdge);
 //m1.insertPointByVertex(2);
 //m1.insertPointByVertex(2);
 m1.mesh3DPoints();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setOFace();
 m1.setSurfaceConfig();

 Simulator3D s2(m1,s1);
 s2.applyLinearInterpolation(mOld);
 s1 = s2;
 s1.setInterfaceGeo();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTKSurface(vtkFolder,"sim",0);
 save.saveBubbleInfo(datFolder);

 PetscFinalize();
 return 0;
}



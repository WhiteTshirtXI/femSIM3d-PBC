// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "TElement.h"
#include "GMRes.h"
#include "InOut.h"
#include "Helmholtz3D.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitializeNoArguments();
 
 int iter = 1;
 real c1 = 1.0;  // lagrangian
 real c2 = 0.5;  // smooth vel
 real c3 = 3.0;  // smooth coord (fujiwara)
 real d1 = 0.0;  // surface tangent velocity u_n=u-u_t 
 real d2 = 0.1;  // surface smooth cord (fujiwara)

 real dt = 0.01;

 string meshFile = "sphere.msh";

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/sphere/zalesak/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 const char *mesh1 = mesh;
 m1.readMSH(mesh1);
 m1.setInterfaceBC();
 m1.setTriEdge();
 m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setOFace();
 m1.setSurfaceConfig();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();

 s1(m1);

 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setD1(d1);
 s1.setD2(d2);
 s1.setDt(dt);

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 3000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.stepImposedPeriodicField("rotating",0.0);
   s1.stepALE();
   s1.movePoints();
   s1.setInterfaceGeo();

   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();
  m1.initMeshParameters();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints();

  // surface operations
  m1.smoothPointsByCurvature();

  m1.insertPointsByLength();
  //m1.insertPointsByCurvature();
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance();
  m1.contractEdgeByLength();
  //m1.removePointsByLength();
  //m1.flipTriangleEdge();

  m1.removePointByNeighbourCheck();
  //m1.checkAngleBetweenPlanes();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();

  Simulator3D s2(m1,s1);
  s1 = s2;
  s1.setCentroidVelPos();

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}



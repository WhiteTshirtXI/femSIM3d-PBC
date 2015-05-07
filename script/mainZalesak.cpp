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
 /* This test case applies a prescribed 2D vortex field in a unit cube to
  * test the re-meshing techinique of the surface mesh. 
  *
  * OBS.: - comment stepSL() on Simulator3D::stepALE
  *       - switch to tetrahedralize( (char*) "QYYAp",&in,&out ) on
  *       Model3D::mesh3DPoints
  *
  * Since the field is prescribed, there is no need of calculating the
  * convection in a Euleurian way (stepSL) and the insertion of nodes on
  * the 3D mesh.
  *
  * */

 PetscInitializeNoArguments();
 
 int iter = 1;
 double d1 = 0.0;  // surface tangent velocity u_n=u-u_t 
 double d2 = 0.1;  // surface smooth cord (fujiwara)

 double dt = 0.01;
 double time = 0.0;

 string meshFile = "sphere.msh";

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
 m1.mesh2Dto3D("QYYAp");
 m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setSurfaceConfig();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();

 s1(m1);

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

   time = s1.getTime();

   s1.stepImposedPeriodicField("rotating",0.0,time);
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

  // surface operations
  //m1.smoothPointsByCurvature();

  m1.insertPointsByLength("flat");
  m1.contractEdgesByLength("flat");
  //m1.removePointsByLength();
  //m1.flipTriangleEdges();

  //m1.removePointsByNeighbourCheck();
  //m1.checkAngleBetweenPlanes();
  /* **************************************** */

  m1.setInterfaceBC();
  m1.setMiniElement();
  m1.restoreMappingArrays();
  m1.setSurfaceVolume();
  m1.setSurfaceArea();
  m1.triMeshStats();

  // computing velocity field X^(n+1),time+1 at new nodes too!
  Simulator3D s2(m1,s1);
  s2.stepImposedPeriodicField("rotating",0.0,time);
  s2.saveOldData();
  s1 = s2;
  s1.setCentroidVelPos();

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}



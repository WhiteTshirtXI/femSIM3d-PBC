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
 /* This test case applies a prescribed vortex field in a unit cube to
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
 real d1 = 0.0;   // surface tangent velocity u_n=u-u_t 
 real d2 = 0.0;   // surface smooth cord (fujiwara)

 real dt = 0.02;
 real T = 3.0;
 real time = 0;

 string meshFile = "sphere.msh";
 //string meshFile = "sphere.msh";

 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/sphere/vortex/" + meshFile;
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

 s1.setDt(dt);

 s1.setD1(d1);
 s1.setD2(d2);

 // initial conditions
 s1.stepImposedPeriodicField("3d",T,s1.getTime()); // X,Y and Z --> Sol(n+1)

 int nReMesh = 1;
 while( time < T )
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

   // time step: n+1/4
   Simulator3D s20(m1,s1);
   real stepTime = dt/4.0;
   s20.stepImposedPeriodicField("3d",T,time+stepTime,stepTime); // SolOld(n) --> Sol(n+1/2)
   s20.saveOldData();        // Sol(n+1/2) --> SolOld(n+1/2)

   // time step: n+1/2
   Simulator3D s30(m1,s20);
   stepTime = dt/2.0;
   s30.stepImposedPeriodicField("3d",T,time+stepTime,stepTime); // SolOld(n) --> Sol(n+1/2)
   s30.saveOldData();        // Sol(n+1/2) --> SolOld(n+1/2)
   s30.stepALE();         // SolOld(n+1/2) --> ALE(n+1/2)

   // time step: n using ALE(n+1/2)
   s1.setUALE(s30.getUALE());
   s1.setVALE(s30.getVALE());
   s1.setWALE(s30.getWALE());
   s1.movePoints();

   real field = cos(3.14159265358*time/T);
   cout << endl;
   cout << "                             | T:        " << T << endl;
   cout << "                             | dt:       " << dt << endl;
   cout << "                             | time:     " << time << endl;
   cout << "                             | iter:     " << iter << endl;
   cout << "                             | field:    " << field << endl;
   cout << endl;

   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);

   s1.saveOldData(); // Sol(n+1) --> SolOld(n)

   s1.timeStep();

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

  m1.setMiniElement();
  m1.restoreMappingArrays();

  // computing velocity field X^(n+1),time+1 at new nodes too!
  Simulator3D s2(m1,s1);
  s2.stepImposedPeriodicField("3d",T,time); // X,Y and Z --> Sol(n+1)
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



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
 real c1 = 0.0;   // lagrangian
 real c2 = 1.0;   // smooth vel
 real c3 = 10.0;  // smooth coord (fujiwara)
 real d1 = 0.0;   // surface tangent velocity u_n=u-u_t 
 real d2 = 0.1;   // surface smooth cord (fujiwara)

 real dt = 0.01;
 real T = 3.0;
 real time = 0;

 string meshFile = "sphere.msh";
 //string meshFile = "sphere.msh";

 //const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *binFolder  = "./bin/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/sphere/vortex/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 if( (*(argv+1)) == NULL )
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

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

 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setD1(d1);
 s1.setD2(d2);
 }
 else if( strcmp( *(argv+1),"restart") == 0 )
 {
  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;

  // load surface mesh
  string aux = *(argv+2);
  string file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();

  s1(m1);

  // load 3D mesh
  file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.readVTKHeaviside(vtkFile);
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();

  s1(m1);

  iter = s1.loadSolution("./","sim",atoi(*(argv+2)));
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);


 int nReMesh = 1;
 while( time <= T )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   /* predictor-corrector */
   Simulator3D s20(m1,s1);
   s20.setDt(dt/2.0);
   s20.setTime(s1.getTime()+dt/2.0);
   // in: SolOld^(n),X^(n)
   // out: Sol^(n+1/2)
   s20.stepImposedPeriodicField("3d",T,dt/2.0); 
   s20.saveOldData();
   s20.stepALE();

   // dt variavel
   //s20.setDtALESinglePhase()
   //dt = s20.getDt()
   //s1.setDt(dt)

   // with ALE(n+1/2)
   // compute velocity at time step: n+1 using dt
   s1.movePoints(s20.getUALE(),
                 s20.getVALE(),
				 s20.getWALE());
   s1.setInterfaceGeo();
   s1.stepImposedPeriodicField("3d",T); // X,Y and Z --> Sol(n+1)

   time = s1.getTime();
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
   save.saveSol(binFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   save.savePoint(datFolder,0);
   save.savePoint(datFolder,1);
   save.savePoint(datFolder,2);
   save.savePoint(datFolder,3);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData(); // Sol(n+1) --> SolOld(n)

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
  //m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints();

  // surface operations
  //m1.smoothPointsByCurvature();

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



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
 real d2 = 0.3;   // surface smooth cord (fujiwara)

 real dt = 0.03;
 real T = 3.0;

 //string meshFile = "sphereCenterLow.msh";
 string meshFile = "sphere.msh";

 //const char *binFolder  = "./bin/";
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

 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setD1(d1);
 s1.setD2(d2);

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = T/dt;
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

//--------------------------------------------------
//    /* predictor-corrector */
//    Simulator3D s20(m1,s1);
//    s20.setDt(dt/2.0);
//    s20.setTime(s1.getTime()+dt/2.0);
//    // in: SolOld^(n),X^(n)
//    // out: Sol^(n+1/2)
//    s20.stepImposedPeriodicField("3d",T,dt/2.0); 
//    s20.saveOldData();
//    s20.stepALE();
//-------------------------------------------------- 

   /* predictor-multicorrector */
   // points at position n
   // compute velocity of time step: n+1/2 
   Model3D m10(m1);
   Simulator3D s10(m10,s1);
   s10.setTime(s1.getTime()+dt/4.0); // t^(n+1/4)
   // in: SolOld(n)
   // out: Sol(n+1/2)
   s10.stepImposedPeriodicField("3d",T,dt/4.0); // SolOld(n) --> Sol(n+1/4)
   s10.saveOldData();        // Sol(n+1/4) --> SolOld(n+1/4),t=t^(n+1/4)+1/4
   s10.setDt(dt/2.0);        // compute X^(n+1/2) from X^(n+1/4)
   s10.stepALE();            // 1st) SolOld(n+1/4) --> ALE(n+1/4)
   // time step: n+1/2 using ALE(n+1/4)
   // 2nd) result: X^(n+1/2) using dt/2
   s10.movePoints(s10.getUALE(),
                  s10.getVALE(),
				  s10.getWALE());

   // points at X^(n+1/2)
   // compute velocity at time step: n+1/4 
   Model3D m20(m10);
   Simulator3D s20(m20,s10);
   s20.setTime(s10.getTime()+dt/4.0); // t^(n+3/4)
   // in: SolOld(n+1/2),X^(n+1/2)
   // out: Sol(n+3/4)
   s20.stepImposedPeriodicField("3d",T,dt/4.0); 
   s20.saveOldData();   // Sol(n+3/4) --> SolOld(n+3/4)
   s20.setDt(dt/2.0);   // compute X^(n) from X^(n+1/2)
   s20.stepALE();       // SolOld(n+3/4) --> ALE(n+3/4)


   // with ALE(n+3/4)
   // compute velocity at time step: n+1 using dt
   s1.movePoints(s20.getUALE(),
                 s20.getVALE(),
				 s20.getWALE());
   s1.setInterfaceGeo();
   s1.stepImposedPeriodicField("3d",T); // X,Y and Z --> Sol(n+1)

   real time = s1.getTime();
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

  //m1.removePointByNeighbourCheck();
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



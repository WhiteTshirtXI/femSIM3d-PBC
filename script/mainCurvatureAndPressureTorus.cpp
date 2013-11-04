// =================================================================== //
// this is file maintorus/curvature.cpp, created at 30-Sep-2011              //
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

 double Re = 100;
 double Sc = 1;
 double We = 1;
 double Fr = 1;
 double c1 = 0.0;  // lagrangian
 double c2 = 0.0;  // smooth vel
 double c3 = 0.0;  // smooth coord (fujiwara)
 double d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 double d2 = 0.0;  // surface smooth cord (fujiwara)
 double alpha = 1;

 double mu_in = 1;
 double mu_out = 0.001;

 double rho_in = 1; 
 double rho_out = 0.01;

 double cfl = 0.8;

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";

 //string meshFile = (string) getenv("MESHLIB");
 //meshFile += "/gmsh/3d/torus/curvature/";
 //const char *mesh = meshFile.c_str();

 /* meshes */
 vector<const char*> mesh;
 mesh.resize(7);
 mesh[0]  = "../../db/gmsh/3d/torus/curvature/0.10.msh";
 mesh[1]  = "../../db/gmsh/3d/torus/curvature/0.09.msh";
 mesh[2]  = "../../db/gmsh/3d/torus/curvature/0.08.msh";
 mesh[3]  = "../../db/gmsh/3d/torus/curvature/0.07.msh";
 mesh[4]  = "../../db/gmsh/3d/torus/curvature/0.06.msh";
 mesh[5]  = "../../db/gmsh/3d/torus/curvature/0.05.msh";
 mesh[6]  = "../../db/gmsh/3d/torus/curvature/0.04.msh";

 for( int i=0;i<(int) mesh.size();i++ )
 {
  cout << color(none,magenta,black);
  cout << "____________________________________ Iteration: " 
       << i << endl << endl;
  cout << resetColor();

  Model3D m1;
  Simulator3D s1;

  m1.readMSH(mesh[i]);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();
  m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setGenericBC();

  s1(m1);

  s1.setRe(Re);
  s1.setSc(Sc);
  s1.setWe(We);
  s1.setFr(Fr);
  s1.setC1(c1);
  s1.setC2(c2);
  s1.setC3(c3);
  s1.setD1(d1);
  s1.setD2(d2);
  s1.setAlpha(alpha);
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.setCfl(cfl);
  s1.init();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  //s1.stepLagrangian();
  s1.stepALE();
  s1.setDtALETwoPhase();
  s1.movePoints();
  s1.assemble();
  s1.matMount();
  s1.setUnCoupledBC();
  s1.setRHS();
  //s1.setInterface();
  s1.setInterfaceGeo();
  s1.unCoupled();

  InOut save(m1,s1); // cria objeto de gravacao
  save.saveMSH(mshFolder,"newMesh",i);
  save.saveVTK(vtkFolder,"sim",i);
  save.saveVTKSurface(vtkFolder,"sim",i);
  save.saveKappaErrorTorus(datFolder);
  save.saveBubbleInfo(datFolder);
  save.chordalPressure(datFolder,"chordalPressure",i);
  save.crossSectionalPlane(datFolder,"XZ",i);

  cout << color(none,magenta,black);
  cout << "________________________________________ END of " 
       << i << endl << endl;;
  cout << resetColor();
 }

 PetscFinalize();
 return 0;
}



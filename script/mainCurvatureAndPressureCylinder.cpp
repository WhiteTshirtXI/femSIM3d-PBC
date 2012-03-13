// =================================================================== //
// this is file maincylinder/curvature.cpp, created at 30-Sep-2011              //
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

 real Re = 100;
 real Sc = 1;
 real We = 1;
 real Fr = 1;
 real c1 = 0.0;  // lagrangian
 real c2 = 0.0;  // smooth vel
 real c3 = 0.0;  // smooth coord (fujiwara)
 real d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 real d2 = 0.0;  // surface smooth cord (fujiwara)
 real alpha = 1;
 real beta = 1;

 real sigma = 1.0;

 real mu_in = 1;
 real mu_out = 0.001;

 real rho_in = 1; 
 real rho_out = 0.01;

 real cfl = 0.8;

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";

 //string meshFile = (string) getenv("MESHLIB");
 //meshFile += "/gmsh/3d/cylinder/curvature/";
 //const char *mesh = meshFile.c_str();

 /* meshes */
 vector<const char*> mesh;
 mesh.resize(43);
 mesh[0]  = "../../db/gmsh/3d/cylinder/curvature/0.90.msh";
 mesh[1]  = "../../db/gmsh/3d/cylinder/curvature/0.85.msh";
 mesh[2]  = "../../db/gmsh/3d/cylinder/curvature/0.80.msh";
 mesh[3]  = "../../db/gmsh/3d/cylinder/curvature/0.75.msh";
 mesh[4]  = "../../db/gmsh/3d/cylinder/curvature/0.70.msh";
 mesh[5]  = "../../db/gmsh/3d/cylinder/curvature/0.65.msh";
 mesh[6]  = "../../db/gmsh/3d/cylinder/curvature/0.60.msh";
 mesh[7]  = "../../db/gmsh/3d/cylinder/curvature/0.55.msh";
 mesh[8]  = "../../db/gmsh/3d/cylinder/curvature/0.50.msh";
 mesh[9]  = "../../db/gmsh/3d/cylinder/curvature/0.45.msh";
 mesh[10]  = "../../db/gmsh/3d/cylinder/curvature/0.40.msh";
 mesh[11]  = "../../db/gmsh/3d/cylinder/curvature/0.35.msh";
 mesh[12]  = "../../db/gmsh/3d/cylinder/curvature/0.30.msh";
 mesh[13]  = "../../db/gmsh/3d/cylinder/curvature/0.25.msh";
 mesh[14]  = "../../db/gmsh/3d/cylinder/curvature/0.24.msh";
 mesh[15]  = "../../db/gmsh/3d/cylinder/curvature/0.23.msh";
 mesh[16]  = "../../db/gmsh/3d/cylinder/curvature/0.22.msh";
 mesh[17]  = "../../db/gmsh/3d/cylinder/curvature/0.21.msh";
 mesh[18]  = "../../db/gmsh/3d/cylinder/curvature/0.20.msh";
 mesh[19]  = "../../db/gmsh/3d/cylinder/curvature/0.19.msh";
 mesh[20]  = "../../db/gmsh/3d/cylinder/curvature/0.18.msh";
 mesh[21]  = "../../db/gmsh/3d/cylinder/curvature/0.17.msh";
 mesh[22]  = "../../db/gmsh/3d/cylinder/curvature/0.16.msh";
 mesh[23] = "../../db/gmsh/3d/cylinder/curvature/0.15.msh";
 mesh[24] = "../../db/gmsh/3d/cylinder/curvature/0.14.msh";
 mesh[25] = "../../db/gmsh/3d/cylinder/curvature/0.13.msh";
 mesh[26] = "../../db/gmsh/3d/cylinder/curvature/0.12.msh";
 mesh[27] = "../../db/gmsh/3d/cylinder/curvature/0.11.msh";
 mesh[28] = "../../db/gmsh/3d/cylinder/curvature/0.10.msh";
 mesh[29] = "../../db/gmsh/3d/cylinder/curvature/0.095.msh";
 mesh[30] = "../../db/gmsh/3d/cylinder/curvature/0.09.msh";
 mesh[31] = "../../db/gmsh/3d/cylinder/curvature/0.085.msh";
 mesh[32] = "../../db/gmsh/3d/cylinder/curvature/0.08.msh";
 mesh[33] = "../../db/gmsh/3d/cylinder/curvature/0.075.msh";
 mesh[34] = "../../db/gmsh/3d/cylinder/curvature/0.07.msh";
 mesh[35] = "../../db/gmsh/3d/cylinder/curvature/0.065.msh";
 mesh[36] = "../../db/gmsh/3d/cylinder/curvature/0.06.msh";
 mesh[37] = "../../db/gmsh/3d/cylinder/curvature/0.055.msh";
 mesh[38] = "../../db/gmsh/3d/cylinder/curvature/0.05.msh";
 mesh[39] = "../../db/gmsh/3d/cylinder/curvature/0.045.msh";
 mesh[40] = "../../db/gmsh/3d/cylinder/curvature/0.04.msh";
 mesh[41] = "../../db/gmsh/3d/cylinder/curvature/0.035.msh";
 mesh[42] = "../../db/gmsh/3d/cylinder/curvature/0.03.msh";

 for( int i=0;i<(int) mesh.size();i++ )
 {
  cout << color(none,magenta,black);
  cout << "____________________________________ Iteration: " 
       << i << endl << endl;
  cout << resetColor();

  // set each bubble length
  vector< real > triEdge;
  triEdge.resize(2);
  triEdge[0] = 0.6; // wall
  triEdge[1] = 0.90-0.01*i*2; // bubble 1 

  Model3D m1;
  Simulator3D s1;

  m1.readMSH(mesh[i]);
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
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setWallBC();

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
  s1.setBeta(beta);
  s1.setSigma(sigma);
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.setCfl(cfl);
  s1.init();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  //s1.stepLagrangian();
  //s1.stepALE();
  s1.stepALEVel();
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
  save.saveVTKHalf(vtkFolder,"simCutPlane",i);
  save.saveVTKSurface(vtkFolder,"sim",i);
  save.saveBubbleInfo(datFolder);

  cout << color(none,magenta,black);
  cout << "________________________________________ END of " 
       << i << endl << endl;;
  cout << resetColor();
 }

 PetscFinalize();
 return 0;
}



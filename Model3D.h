/**=========================================================================
 *\file Id:Model3D.h, GESAR
 *\author Gustavo Rabello dos Anjos
 *\date   18-jun-2007
 *\comment Classe Model 3D para o modulo HYDRO
 *\references 
 *\version 1.0
 *\updates
 	version    date         author             comment
	1.0        23/08/2007   Gustavo            inicio de implementacao
 =========================================================================**/

#ifndef Model3D_H
#define Model3D_H

#include <list>
#include <vector>
#include <string.h>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"
#include "TElement.h"
#include "math.h"
#include "compare.h"
#include "colors.h"
#include "structs.h"
#include "geometry.h"
#include "tetgen.h"

/**
 * @brief classe responsavel pela preparacao da malha para entrada no
 * simulador. Definicao das condicoes de contorno, leitura de arquivos,
 * inclusao de vertices auxiliares nos elementos etc.
 *
 * @return 
 **/
class Model3D
{
 public:
  Model3D();                      // construtor padrao
  virtual ~Model3D();             // destrutor padrao

  // reading files
  void readVTK( const char* filename );
  void readVTKCC( const char* filename );
  void readVTKSurface( const char* filename );
  void readMSH( const char* filename );
  void readBC( const char* filename );
  void readBaseStateNu(const char* _filename);

  void setSphere(real _xC,real _yC,real _zC,real _r,real _eps);
  void setCube(real _lim1,real _lim2,real _eps);
  void setCube(real _xlimInf,real _xlimSup,
               real _ylimInf,real _ylimSup,
	       	   real _zlimInf,real _zlimSup,real _eps);

  // surface points treatment
  void setSurface();
  void setSurfaceFace();
  void setSurfaceFace2();
  void setSurfaceTri();
  void setConvexTri();
  void buildSurfMesh();
  SurfaceMesh arrangeMesh(SurfaceMesh _mesh,int _nVerts,int _begin);
  void setInOutVert();
  void setInOutElem();
  void computeAverageTriangleEdge();
  void insertPointsByLength();
  void insertPointsByInterfaceDistance();
  void removePointsByLength();
  void insertPointsByArea();
  void surfaceTriangulator(int _v);
  void surfaceTriangulatorEarClipping(int _v);
  void surfaceTriangulatorQualityEarClipping(int _v);
  void deleteSurfacePoint(int _v);
  void markSurfElemForDeletion(int _elem);
  void deleteSurfaceElements();
  void insertPoint(int _edge);
  void insertPointWithCurvature(int _edge);
  void insertPointsBetweenBubblesByPosition();
  void deletePoint(int _v);
  void setPolyhedron(int _v);
  void flipTriangleEdge();
  void contractEdgeByLength();
  int findEdge(int _v1,int _v2);
  void removePointsByInterfaceDistance();
  void remove3dMeshPointsByDistance();
  void breakup();
  clVector triangleQuality(int _v);
  void setMapEdgeTri();
  void setSurfaceConfig();
  bool checkNormal(int _surfaceNode,int _v1,int _v2,int _vIn);
  void saveVTK( const char* _dir,const char* _filename, int _iter );
  void saveVTKConvex( const char* _dir,const char* _filename, int _iter );
  void saveVTKSurface( const char* _dir,const char* _filename, int _iter );
  bool testFace(int v1, int v2, int v3, int v4);

  // meshing with TETGEN
  void setMeshStep(int nX,int nY,int nZ);
  void setMeshDisk(int nLados1Poli,int nCircMax,int nZ);
  void mesh2Dto3D();
  void mesh2Dto3DOriginal();
  void mesh3DPoints();
  bool checkMeshQuality(tetgenio &_tetmesh);
  tetgenio convertSurfaceMeshToTetGen(SurfaceMesh _mesh,tetgenio &_tetmesh);
  Mesh3D convertTetgenToMesh3d(tetgenio &_tetmesh);
  void convertTetgenToModel3D(tetgenio &_tetmesh);
  void convertModel3DtoTetgen(tetgenio &_tetmesh);

  // boundary condition settings
  void setNuCteDiskBC();
  void setNuCDiskBC();
  void setNuZDiskBC();
  void readAndSetVelocityDiskBC(const char* _dir,const char* _filename);
  void readAndSetPressureDiskBC(const char* _dir,const char* _filename);
  void setCDiskBC();
  void setDiskFSBC();
  void setDiskCouetteBC();
  void setInterfaceBC();
  void setBubbleBubbleBC();
  void setBubbleBC2();
  void setBubble3DBC();
  void setCubeBC();
  void setCubeBC2();
  void setWallBC();
  void setWallAnnularBC();
  void set2BubbleBC();
  void setStepBC();
  void setCStepBC();
  void setStepReservBC();
  void setStepReservInvBC();
  void setCouetteBC();

  // adimensionalisation
  void setAdimenDiskCouette();
  void setAdimenDisk();
  void setAdimenStep();

  void setPerturbSurf();
  void setPerturbSurf2();
  void setPerturbSurfSquare();

  // misc
  void moveXPoints(clVector &_vec,real _dt);
  void moveYPoints(clVector &_vec,real _dt);
  void moveZPoints(clVector &_vec,real _dt);
  void setMiniElement();            
  void setQuadElement();             
  void setNeighbour();
  void setNeighbourSurface();
  void setVertNeighbour();
  void setOFace();
  void printMeshReport(tetgenio &_tetmesh);
  void clearBC();
  void reAllocStruct();
  void computeSurfaceNormal();
  void computeSurfaceAverageNormal();
  void setKappaSurface();
  void setKappaSurface(clVector &_kappa);
  void setCloser();
  void setInterfaceDistance();
  void computeKappaGeo();

  // get and set methods
  clVector* getX();
  clVector* getXVert();
  real getMaxX();
  real getMinX();
  void setX(clVector _X);
  clVector* getY();
  clVector* getYVert();
  real getMaxY();
  real getMinY();
  void setY(clVector _Y);
  real getMaxZ();
  real getMinZ();
  clVector* getZ();
  clVector* getZVert();
  void setZ(clVector _Z);

  clVector* getUC();
  clVector* getVC();
  clVector* getWC();
  clVector* getPC();
  clVector* getCC();
  clVector* getCCOriginal();
  clVector* getOutflow();
  clVector* getIdbcu();
  clVector* getIdbcv();
  clVector* getIdbcw();
  clVector* getIdbcp();
  clVector* getIdbcc();
  clMatrix* getIEN();
  clVector* getInterfaceDistance();
  clDMatrix* getCurvature();
  SurfaceMesh* getSurfMesh();
  SurfaceMesh* getInterfaceMesh();
  SurfaceMesh* getConvexMesh();
  Mesh3D* getMesh3d();
  void setMeshX(clVector _X); 
  void setMeshY(clVector _Y);
  void setMeshZ(clVector _Z);
  int getNumVerts();
  int getNumNodes();
  int getNumElems();
  real getVolume(int _v1,int _v2,int _v3,int _v4);
  real getVolume(int _elem);
  real getAreaVert(int _v1,int _v2,int _v3);
  real getAreaElem(int _elem);
  real getAreaHeron(int _elem);
  real getLength(int _v1,int _v2);
  clMatrix* getOFace();
  clVector* getSurface();
  clVector* getNonSurface();
  real getXCenter();
  real getYCenter();
  real getZCenter();
  real getBubbleRadius();
  vector< list<int> >* getNeighbourElem();
  vector< list<int> >* getNeighbourVert();
  vector< list<int> >* getNeighbourFace();
  vector< list<int> >* getNeighbourFaceVert();
  vector< list<int> >* getElemSurface();
  vector< list<int> >* getSurfaceViz();
  vector< list<int> >* getFaceIEN();
  list<int>* getOutVert();
  list<int>* getInVert();
  list<int>* getOutElem();
  list<int>* getInElem();
  real getTriEdge();
  void setTriEdge(real _triEdge);

  void operator=(Model3D &_mRight);

 private:
  tetgenio in,mid,out;
  clVector uc,vc,wc,pc,cc;
  clMatrix IEN;
  clDMatrix curvature;
  Mesh3D mesh3d;
  SurfaceMesh surfMesh,interfaceMesh,convexMesh;
  clVector X,Y,Z;
  clVector xConvex,yConvex,zConvex;
  clVector outflow,idbcu,idbcv,idbcw,idbcp,idbcc;
  clMatrix faceFace,freeFace,mapViz;
  clMatrix oFace;
  clVector V; // vetor de volumes dos elementos de malha 
  clVector idRegion;
  clVector surface,nonSurface;
  clMatrix mapEdgeTri;
  clVector xSurface,ySurface,zSurface;
  clVector closer,xCloser,yCloser,zCloser,closerViz;
  clVector interfaceDistance;

  int numVerts;                   // numero total de vertices da malha
  int numElems;                   // numero total de elementos da malha
  int numNodes;                   // numero total de nos da malha
  int dVerts;                     // delta vertices (v_new-v_old)
  int numTriangles;
  real rMax;                      // tamanho max do raio do disco
  real xCenter,yCenter,zCenter;
  real bubbleRadius;
  real triEdge,averageTriEdge;
  int isp;                        // isp: num of inserted surface points
  int rsp;                        // rsp: num of removed surface points
  int csp;                        // csp: num of contracted surface points
  int ip;                         // ip: num of inserted 3d mesh points
  int rp;                         // rp: num of removed 3d mesh points
  int rpi;                        // rpi: by interface distance
  int flip;
  int badtet;                     // num of shit tetrahedrons

  vector< list<int> > neighbourElem;  // lista de elementos de cada no
  vector< list<int> > neighbourVert;  // lista de vizinhos de cada no
  vector< list<int> > neighbourFace;  // lista de vizinhos de cada no
  vector< list<int> > faceIEN;
  vector< list<int> > neighbourFaceVert;
  vector< list<int> > elemSurface;
  vector< list<int> > surfaceViz;  // lista de vizinhos na interface
  vector< list<int> > neighbourSurfaceElem; // lista de elementos triangulares
  vector< list<int> > neighbourPoint;  // lista de pontos vizinhos da superficie
  list<int> outVert,inVert; // lista de elementos do interior 
  list<int> outElem,inElem; // lista de elementos do interior
};

#endif

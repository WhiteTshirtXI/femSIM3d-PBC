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
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
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
  Model3D(const Model3D &_mRight);
  virtual ~Model3D();             // destrutor padrao

  // reading files
  void readVTK( const char* filename );
  void readVTKHeaviside( const char* filename );
  void readVTKSurface( const char* filename );
  void readMSH( const char* filename );
  void readBC( const char* filename );

  void setSphere(real _xC,real _yC,real _zC,real _r,real _eps);
  void setCube(real _lim1,real _lim2,real _eps);
  void setCube(real _xlimInf,real _xlimSup,
               real _ylimInf,real _ylimSup,
	       	   real _zlimInf,real _zlimSup,real _eps);

  // surface points treatment
  void setSurface();
  void setSurfaceTri();
  void setConvexTri();
  void buildSurfMesh();
  SurfaceMesh arrangeMesh(SurfaceMesh _mesh,int _nVerts,int _begin);
  void setInOutVert();
  void setInOutElem();
  void insertPointsByLength(const char* _interpolation);
  void insertPointsByLength(const char* _interpolation,real _param);
  void insertPointsByCurvature(const char* _interpolation);
  void removePointsByCurvature();
  void smoothPoint(int _node);
  void smoothPointsByCurvature();
  void checkAngleBetweenPlanes();
  void insertPointsByInterfaceDistance(const char* _interpolation);
  void removePointsByLength();
  void removePointsByLength(real _param);
  void removePointsByNeighbourCheck();
  void removePointByNeighbourCheck(int _node);
  void insertPointsByArea();
  void surfaceTriangulator(int _v);
  void surfaceTriangulatorEarClipping(int _v,
	                                  list<int> _list,
	                                  const char* _mode);
  void deleteSurfacePoint(int _v);
  void markSurfElemForDeletion(int _elem);
  void deleteSurfaceElements();
  void insertSurfacePoint(int _edge,const char* _interpolation);
  void removeSurfacePoint(int _node);
  void insertPointsBetweenBubblesByPosition();
  list<int> setPolyhedron(list<int> _myList);
  void flipTriangleEdges();
  void contractEdgesByLength(const char* _interpolation);
  void contractEdgesByLength(const char* _interpolation,real _param);
  void contractEdgesByLength2(const char* _interpolation,real _param);
  int findEdge(int _v1,int _v2);

  // 3D points treatment
  void mark3DPointForDeletion(int _vert);
  void delete3DPoints();
  void removePointsByInterfaceDistance();
  void remove3dMeshPointsByDistance();
  void remove3dMeshPointsByDiffusion();
  void insert3dMeshPointsByDiffusion();
  void insert3dMeshPointsByVolume();
  void remove3dMeshPointsByVolume();
  void remove3dMeshPointsByHeight();

  void breakup();
  clVector triangleQuality(int _v);
  void setMapEdge();
  void setMapEdgeTri();
  void setSurfaceConfig();
  void saveVTK( const char* _dir,const char* _filename, int _iter );
  void saveVTKConvex( const char* _dir,const char* _filename, int _iter );
  void saveVTKSurface( const char* _dir,const char* _filename, int _iter );
  bool testFace(int v1, int v2, int v3, int v4);

  // Volume calculation
  void setInitSurfaceVolume();
  void setSurfaceVolume();
  real getSurfaceVolume(int _region);
  real getSurfaceVolumeTET(int _region);

  // Area calculation
  void setInitSurfaceArea();
  void setSurfaceArea();
  real getSurfaceArea(int _region);

  void applyBubbleVolumeCorrection();
  clVector computeConvexRegionCentroid(int _region);
  clVector computeConvexRegionCentroid2D(real _zPlane);

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

  vector<int> getOPER();
  vector<int> getOPERSURF();

  vector<int> getISP();
  vector<int> getISPC();
  vector<int> getRSP();
  vector<int> getRSPN();
  vector<int> getRSPC();
  vector<int> getFLIP();
  vector<int> getSPC();
  vector<int> getSPP();
  vector<int> getINTET();
  vector<real> getMinArea();
  vector<real> getMaxArea();
  vector<real> getMinLength();
  vector<real> getMaxLength();
  vector<int> getIdMinArea();
  vector<int> getIdMaxArea();
  vector<int> getNumSurfElems();
  vector<int> getNumSurfVerts();

  vector<int> getIP();
  vector<int> getIPD();
  vector<int> getRP();
  vector<int> getRPI();
  vector<int> getRPD();
  vector<int> getRPDist();
  vector<int> getRPH();
  vector<int> getRPV();
  vector<int> getCSP();
  vector<real> getMinVolume();
  vector<real> getMaxVolume();
  vector<int> getIdMinVolume();
  vector<int> getIdMaxVolume();

  // boundary condition settings
  void setInfiniteDiskBC(real _F,real _G,real _H);
  void setFiniteDiskBC();
  void setCDiskBC();
  void setDiskFSBC();
  void setDiskCFSBC();
  void setDiskCouetteBC();
  void setInterfaceBC();
  void setCubeBC();
  void setGenericBC();
  void setGenericBC(real _vel);
  void setWallBC();
  void setWallBC(real _vel);
  void setMovingWallBC();
  void setMovingWallBC(real _vel);
  void setMicroWallBC();
  void setCircularWallBC();
  void setWallCouetteBC();
  void setWallAnnularBC();
  void setWallInterfaceBC();
  void set2BubblesBC();
  void set2AxiBubblesBC();
  void setStepBC();
  void setWallStepBC();
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
  void setSphereToEllipsoid(real _factor);
  void setBiggerSphere(real _factor);

  // misc
  void moveXPoints(clVector &_vec,real _dt);
  void moveYPoints(clVector &_vec,real _dt);
  void moveZPoints(clVector &_vec,real _dt);
  void setMiniElement();            
  void setCentroid();
  void centroidPositionCorrection();
  void edgeMidPointPositionCorrection();
  void checkTetrahedronOrientation();
  void checkTriangleOrientation();
  void setQuadElement();             
  void setNeighbour();
  void setNeighbourSurfaceElem();
  list<int> getNeighbourSurfacePoint(int _node);
  void setNeighbourSurfacePoint();
  void setVertNeighbour();
  void setOFace();
  void triMeshStats();
  void tetMeshStats();
  void clearBC();
  void reAllocStruct();
  void computeSurfaceNormal();
  void computeSurfaceAverageNormal();
  void setKappaSurface();
  void setKappaSurface(clVector &_kappa);
  void setCloser();
  void setInterfaceDistance();
  void restoreMappingArrays();
  clVector getNormalAndKappa(int _node,list<int> _myList);
  //clVector getNormalAndKappaByGauss(int _node,list<int> _myList);
  clVector getNormalAndKappaByDesbrun(int _node,list<int> _myList);
  void setNormalAndKappa();
  clVector getNormalElem(int _elem);
  clVector getNormalElem(int _v1,int _v2,int _v3);

  // annular boundary conditions
  void checkLineOrientation();
  void setNormalAndKappa2D();

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
  clVector* getEdgeSize();
  void setEdgeSize(clVector _edgeSize);
  vector<real> getAverageTriLength();
  vector<real> getAverageTriArea();
  vector<real> getAverageTetVolume();

  clVector* getUC();
  clVector* getVC();
  clVector* getWC();
  clVector* getPC();
  clVector* getCC();
  clVector* getHeaviside();
  clVector* getOutflow();
  clVector* getIdbcu();
  clVector* getIdbcv();
  clVector* getIdbcw();
  clVector* getIdbcp();
  clVector* getIdbcc();
  clMatrix* getIEN();
  clMatrix* getMapEdge();
  clMatrix* getMapEdgeTri();
  clVector* getInterfaceDistance();
  clDMatrix* getCurvature();
  clVector* getElemIdRegion();
  clVector* getVertIdRegion();
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
  real getLength(int _v1,int _v2);
  clMatrix* getOFace();
  clVector* getSurface();
  clVector* getNonSurface();
  real getXCenter();
  real getYCenter();
  real getZCenter();
  vector< list<int> >* getNeighbourElem();
  vector< list<int> >* getNeighbourVert();
  vector< list<int> >* getNeighbourFace();
  vector< list<int> >* getNeighbourPoint();
  vector< list<int> >* getFaceIEN();
  list<int>* getBoundaryVert();
  list<int>* getInVert();
  list<int>* getOutElem();
  list<int>* getInElem();
  real getMinEdge();
  void setTriEdge();
  void setTriEdge(vector< real > _triEdge);
  void initMeshParameters();
  vector<real> getTriEdge();
  void setTetVol(vector< real > _tetVol);
  vector<real> getTetVol();
  vector<real> getInitSurfaceArea();
  vector<real> getInitSurfaceVolume();
  vector<real> getSurfaceArea();
  vector<real> getDArea();
  vector<real> getErrorArea();
  vector<real> getSurfaceVolume();
  vector<real> getDVolume();
  vector<real> getErrorVolume();
  void setSingleElement();
  void setTwoElements();
  void setThreeElements();
  void setFourElements();
  void integralParabolic();
  clVector* getCloser();

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
  clVector heaviside;
  clVector surface,nonSurface;
  clMatrix mapEdgeTri,mapEdge;  
  clVector xSurface,ySurface,zSurface;
  clVector closer,xCloser,yCloser,zCloser;
  clVector interfaceDistance;
  clVector vertIdRegion,elemIdRegion;
  clVector edgeSize; // edge lenght vector :: Dim = number of vertices

  int numVerts;                   // numero total de vertices da malha
  int numElems;                   // numero total de elementos da malha
  int numNodes;                   // numero total de nos da malha
  int dVerts;                     // delta vertices (v_new-v_old)
  real rMax;                      // tamanho max do raio do disco
  real xCenter,yCenter,zCenter;
  real minEdge;
  real minEdgeTri;

  vector<int> oper,opersurf; // oper: num of operations of time step
  vector<int> isp;           // isp: num of inserted surface points by length
  vector<int> ispc;          // ispc: num of inserted surface points by curv
  vector<int> rsp;           // rsp: num of removed surface points by length
  vector<int> rspn;          // rspn: num of removed surface points by neigh check
  vector<int> rspc;          // rspc: num of removed surface points by curv
  vector<int> csp;           // csp: num of contracted surface points
  vector<int> flip;          // flip: flipping operations
  vector<int> spc;           // spc: smoothing operations
  vector<int> spp;           // spp: smoothing operations
                            
  vector<int> ip;            // ip: num of inserted 3d mesh points
  vector<int> ipd;           // ipd: by diffusion 
  vector<int> rp;            // rp: num of removed 3d mesh points
  vector<int> rpi;           // rpi: by interface distance
  vector<int> rpv;           // rpv: by volume 
  vector<int> rpd;           // rpd: by diffusion 
  vector<int> rpdist;        // rpd: by distance 
  vector<int> rph;           // rph: by height 
  vector<int> badtet;        // num of shit tetrahedrons

  vector<int> intet;         // csp: num of surface tetrahedrons 
  vector<int> idMinVolume;   // ID of min tet volume
  vector<int> idMaxVolume;   // ID of max tet volume
  vector<real> minVolume;    // min tet volume
  vector<real> maxVolume;    // max tet volume
  vector<int> idMinArea;     // ID of min tri area 
  vector<int> idMaxArea;     // ID of max tri area 
  vector<real> minArea;      // min triangle area 
  vector<real> maxArea;      // max triangle area
  vector<real> minLength;    // min triangle length
  vector<real> maxLength;    // max triangle length
  vector<int> numSurfElems;  // number of surface elements
  vector<int> numSurfVerts;  // number of surface points                            

  vector<real> initSurfaceVolume,surfaceVolume; // vector de volumes iniciais
  vector<real> dVolume,errorVolume;  
  vector<real> initSurfaceArea,surfaceArea; // vector de areas iniciais
  vector<real> dArea,errorArea;  
  vector<real> triEdge;               // vector de tamanho de aresta 
  vector<real> averageTriLength;    // vector de tamanho de aresta medio
  vector<real> averageTriArea;      // vector de tamanho de aresta medio
  vector<real> averageTetVolume;    // vector de volume medio
  vector<real> tetVol;                // vector de volumes 
  vector< list<int> > neighbourElem;  // lista de elementos de cada no
  vector< list<int> > neighbourVert;  // lista de vizinhos de cada no
  vector< list<int> > neighbourFace;  // lista de vizinhos de cada no
  vector< list<int> > faceIEN;
  vector< list<int> > neighbourSurfaceElem; // lista de elementos triangulares
  vector< list<int> > neighbourPoint;  // lista de pontos vizinhos da superficie
  vector< list<int> > neighbourEdge; // elems que compart. a face do triangulo 
  list<int> boundaryVert,inVert; // lista de elementos do interior 
  list<int> outElem,inElem; // lista de elementos do interior
};

#endif

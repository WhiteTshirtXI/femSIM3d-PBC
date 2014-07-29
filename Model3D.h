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
#include <set>
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
#include "MathTools.h"

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

  // surface points treatment
  void setSurface();
  void setSurfaceTri();
  SurfaceMesh arrangeMesh(SurfaceMesh _mesh,int _nVerts,int _begin);
  void setInOutVert();
  void setInOutElem();
  void insertPointsByLength(const char* _interpolation);
  void insertPointsByLength(const char* _interpolation,double _param);
  void insertPointsByCurvature(const char* _interpolation);
  void removePointsByCurvature();
  void smoothPoint(int _node);
  void smoothPointsByCurvature();
  void checkAngleBetweenPlanes();
  void insertPointsByInterfaceDistance(const char* _interpolation);
  void removePointsByLength();
  void removePointsByLength(double _param);
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
  void contractEdgesByLength(const char* _interpolation,double _param);
  void contractEdgesByLength2(const char* _interpolation,double _param);
  int findEdge(int _v1,int _v2);

  // 3D points treatment
  void mark3DPointForDeletion(int _vert);
  void delete3DPoints();
  void removePointsByInterfaceDistance();
  void remove3dMeshPointsByDistance();
  void remove3dMeshPointsByDiffusion();
  void remove3dMeshPointsByDiffusion(double _param);
  void insert3dMeshPointsByDiffusion();
  void insert3dMeshPointsByDiffusion(double _param);
  void insert3dMeshPointsByVolume();
  void remove3dMeshPointsByVolume();
  void remove3dMeshPointsByHeight();

  void breakup();
  clVector triangleQuality(int _v);
  void setMapping();
  void setMapEdgeTri();
  void setSurfaceConfig();
  void saveVTK( const char* _dir,const char* _filename, int _iter );
  void saveVTKSurface( const char* _dir,const char* _filename, int _vertID, int _iter );
  bool testFace(int v1, int v2, int v3, int v4);

  // Volume calculation
  void setInitSurfaceVolume();
  void setSurfaceVolume();
  double getSurfaceVolume(int _region);
  double getSurfaceVolumeTET(int _region);

  // Area calculation
  void setInitSurfaceArea();
  void setSurfaceArea();
  double getSurfaceArea(int _region);

  void applyBubbleVolumeCorrection();
  clVector computeConvexRegionCentroid(int _region);
  clVector computeConvexRegionCentroid2D(double _zPlane);

  // meshing with TETGEN
  void setMeshStep(int nX,int nY,int nZ);
  void setMeshStep(int nX,int nY,int nZ,const char* _param);
  void setMeshDisk(int nLados1Poli,int nCircMax,int nZ);
  void setMeshDisk(int nLados1Poli,int nCircMax,int nZ,const char* _param);
  void transformDiskToSphere();
  void mesh2Dto3D();
  void mesh2Dto3D(const char* _param);
  void mesh2Dto3DOriginal();
  void mesh2Dto3DOriginal(const char* _param);
  void mesh3DPoints();
  void mesh3DPoints(const char* _param);
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
  vector<double> getMinArea();
  vector<double> getMaxArea();
  vector<double> getMinLength();
  vector<double> getMaxLength();
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
  vector<double> getMinVolume();
  vector<double> getMaxVolume();
  vector<int> getIdMinVolume();
  vector<int> getIdMaxVolume();

  // boundary condition settings
  void setInfiniteDiskBC(double _F,double _G,double _H);
  void setInfiniteSphereBC(double _F,double _G, double _H);
  void setFiniteDiskBC();
  void setCDiskBC();
  void setDiskFSBC();
  void setDiskCFSBC();
  void setDiskCouetteBC();
  void setInterfaceBC();
  void setGenericBC();
  void setGenericBC(double _vel);
  void setWallBC();
  void setWallBC(double _vel);
  void setMovingWallBC();
  void setMovingWallBC(double _vel);
  void setMicroWallBC();
  void setCircularWallBC();
  void setWallCouetteBC();
  void setWallAnnularBC();
  void setWallInterfaceBC();
  void set2BubblesBC();
  void set2AxiBubblesBC();
  void setStepBC();
  void setOnePointPressureBC(); // <<<
  void setWallNormalVWBC(); // << slip condition, except PBC walls
  void setWallSlipEdgesBC();
  void setCubeVortexBC(); // <<< TaylorGreen vortex w/ slip walls
  void setStretchJetMesh(); // geometrical transform
  void setUnstretchJetMesh(); // geometrical transform
  void setWallMovingPBC(double _velInf, double _velSup); 
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
  void setSphereToEllipsoid(double _factor);
  void setBiggerSphere(double _factor);

  // misc
  void moveXPoints(clVector &_vec,double _dt);
  void moveYPoints(clVector &_vec,double _dt);
  void moveZPoints(clVector &_vec,double _dt);
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
  double getMaxX();
  double getMinX();
  void setX(clVector _X);
  clVector* getY();
  clVector* getYVert();
  double getMaxY();
  double getMinY();
  void setY(clVector _Y);
  double getMaxZ();
  double getMinZ();
  clVector* getZ();
  clVector* getZVert();
  void setZ(clVector _Z);
  clVector* getEdgeSize();
  void setEdgeSize(clVector _edgeSize);
  vector<double> getAverageTriLength();
  vector<double> getAverageTriArea();
  vector<double> getAverageTetVolume();

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
  int getNumVerts();
  int getNumNodes();
  int getNumElems();
  double getVolume(int _v1,int _v2,int _v3,int _v4);
  double getVolume(int _elem);
  double getLength(int _v1,int _v2);
  clMatrix* getOFace();
  clVector* getSurface();
  double getXCenter();
  double getYCenter();
  double getZCenter();
  vector< list<int> >* getNeighbourElem();
  vector< list<int> >* getNeighbourVert();
  vector< list<int> >* getNeighbourFace();
  vector< list<int> >* getNeighbourPoint();
  vector< list<int> >* getFaceIEN();
  list<int>* getBoundaryVert();
  list<int>* getInVert();
  list<int>* getOutElem();
  list<int>* getInElem();
  vector<int>* getPbcIndicesLeft(); // PBC
  vector<int>* getPbcIndicesRight(); // PBC
  vector< vector<int> >* getPbcIndicesMaster(); // PBC
  vector< vector<int> >* getPbcIndicesSlave(); // PBC
  double getMinEdge();
  void setTriEdge();
  void setTriEdge(vector< double > _triEdge);
  void initMeshParameters();
  vector<double> getTriEdge();
  void setTetVol(vector< double > _tetVol);
  vector<double> getTetVol();
  vector<double> getInitSurfaceArea();
  vector<double> getInitSurfaceVolume();
  vector<double> getSurfaceArea();
  vector<double> getDArea();
  vector<double> getErrorArea();
  vector<double> getSurfaceVolume();
  vector<double> getDVolume();
  vector<double> getErrorVolume();
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
  clMatrix oFace;
  clVector V; // vetor de volumes dos elementos de malha 
  clVector heaviside;
  clVector surface;
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
  double rMax;                      // tamanho max do raio do disco
  double xCenter,yCenter,zCenter;
  double minEdge;
  double minEdgeTri;

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
  vector<double> minVolume;    // min tet volume
  vector<double> maxVolume;    // max tet volume
  vector<int> idMinArea;     // ID of min tri area 
  vector<int> idMaxArea;     // ID of max tri area 
  vector<double> minArea;      // min triangle area 
  vector<double> maxArea;      // max triangle area
  vector<double> minLength;    // min triangle length
  vector<double> maxLength;    // max triangle length
  vector<int> numSurfElems;  // number of surface elements
  vector<int> numSurfVerts;  // number of surface points                            

  vector<double> initSurfaceVolume,surfaceVolume; // vector de volumes iniciais
  vector<double> dVolume,errorVolume;  
  vector<double> initSurfaceArea,surfaceArea; // vector de areas iniciais
  vector<double> dArea,errorArea;  
  vector<double> triEdge;               // vector de tamanho de aresta 
  vector<double> averageTriLength;    // vector de tamanho de aresta medio
  vector<double> averageTriArea;      // vector de tamanho de aresta medio
  vector<double> averageTetVolume;    // vector de volume medio
  vector<double> tetVol;                // vector de volumes 
  vector< list<int> > neighbourElem;  // lista de elementos de cada no
  vector< list<int> > neighbourVert;  // lista de vizinhos de cada no
  vector< list<int> > neighbourFace;  // lista de vizinhos de cada no
  vector< list<int> > faceIEN;
  vector< list<int> > neighbourSurfaceElem; // lista de elementos triangulares
  vector< list<int> > neighbourPoint;  // lista de pontos vizinhos da superficie
  vector< list<int> > neighbourEdge; // elems que compart. a face do triangulo 
  list<int> boundaryVert,inVert; // lista de elementos do interior 
  list<int> outElem,inElem; // lista de elementos do interior

  /* PBC */
  vector<int> pbcIndicesLeft;
  vector<int> pbcIndicesRight;
  vector< vector<int> > pbcIndicesMaster;
  vector< vector<int> > pbcIndicesSlave; 

};


#endif

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
#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"
#include "math.h"
#include "compare.h"
#include "colors.h"

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

  void readVTK( const char* filename );
  void readVTKCC( const char* filename );
  void readVTKSurface( const char* filename );
  void readMSH( const char* filename );
  void readBC( const char* filename );
  void readBaseStateNu(const char* _filename);
  void setMeshStep(int nX,int nY,int nZ);
  void setStepBC();
  void setCStepBC();
  void setStepReservBC();
  void setStepReservInvBC();
  void setCouetteBC();
  void setAdimenStep();
  void setMeshDisk(int nLados1Poli,int nCircMax,int nZ);
  void setTriEdge();

  // surface points treatment
  void setSurface();
  void setSurfaceFace();
  void setSurfaceTri();
  void setOutTri();
  void setInOutVert();
  void setInOutElem();
  void setTriangleMinEdge();
  void insertPointsByLength();
  void removePointsByLength();
  void insertRemovePointsByLength();
  void insertPointsByArea();
  void surfaceTriangulator(int _v);
  void surfaceTriangulatorEarClipping(int _v);
  void surfaceTriangulatorQualityEarClipping(int _v);
  void deleteSurfacePoint(int _v);
  void deleteSurfaceElementByPoint(int _v);
  void insertPoint(int _v);
  void deletePoint(int _v);
  void setPolyedron(int _v);
  void flipTriangleEdge( int _edge );
  int findEdge(int _v1,int _v2);
  void removePointsByInterfaceDistance();
  void breakup();
  clVector triangleQuality(int _v);
  clVector dsearchn(clVector _X,clVector _Y,clVector _Z,
	                clVector &_XI,clVector &_YI,clVector &_ZI);

  void meshTest();
  void mesh2Dto3D();
  void mesh2Dto3DOriginal();
  void mesh3DPoints();
  void setNuCteDiskBC();
  void setNuCDiskBC();
  void setNuZDiskBC();
  void readAndSetVelocityDiskBC(const char* _dir,const char* _filename);
  void readAndSetPressureDiskBC(const char* _dir,const char* _filename);
  void setCDiskBC();
  void setDiskFSBC();
  void setDiskCouetteBC();
  void setSphere(real _xC,real _yC,real _zC,real _r,real _eps);
  void setCube(real _lim1,real _lim2,real _eps);
  void setCube(real _xlimInf,real _xlimSup,
               real _ylimInf,real _ylimSup,
	       	   real _zlimInf,real _zlimSup,real _eps);
  void setInterfaceBC();
  void setAdimenDiskCouette();
  void setAdimenDisk();
  void setPerturbSurf();
  void setPerturbSurf2();
  void setPerturbSurfSquare();
  void setMiniElement();            
  void setQuadElement();             
  void setNeighbour();
  void setNeighbourSurface();
  void setVertNeighbour();
  void setOFace();
  void setSurfaceConfig();
  bool testFace(int v1, int v2, int v3, int v4);
  void setBubbleBubbleBC();
  void setBubbleBC2();
  void setBubble3DBC();
  void setCubeBC();
  void setCubeBC2();
  void setWallBC();
  void setWallAnnularBC();
  void set2BubbleBC();

  real getMaxAbsUC();
  real getMinAbsUC();
  real getMaxAbsVC();
  real getMinAbsVC();
  real getMaxAbsWC();
  real getMinAbsWC();
  real getDeltaXMin();
  real getDeltaYMin();
  real getDeltaZMin();
  clVector* getX();
  real getMaxX();
  real getMinX();
  void setX(clVector _X);
  clVector* getY();
  real getMaxY();
  real getMinY();
  void setY(clVector _Y);
  real getMaxZ();
  real getMinZ();
  clVector* getZ();
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
  clMatrix* getIENTri();
  clMatrix* getIENConvexTri();
  clMatrix* getIENOriginal();
  int getNumVerts();
  int getNumVertsOriginal();
  int getNumNodes();
  int getNumElems();
  int getNumGLEU();
  int getNumGLEP();
  int getNumGLEC();
  real getVolume(int _elem);
  real getArea(int _v1,int _v2,int _v3);
  real getArea(int _elem);
  real getAreaHeron(int _elem);
  real getLength(int _v1,int _v2);
  void clearBC();
  void reAllocStruct();
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
  void operator=(Model3D &_mRight);
  bool checkNormal(int _surfaceNode,int _v1,int _v2,int _vIn);


 private:
  clVector uc,vc,wc,pc,cc;
  clMatrix IEN,IENTri,IENConvexTri;
  clVector X,Y,Z;
  clVector outflow,idbcu,idbcv,idbcw,idbcp,idbcc;
  clMatrix faceFace,freeFace,mapViz;
  clMatrix oFace;
  clVector V; // vetor de volumes dos elementos de malha 
  clVector idRegion;
  clVector surface,nonSurface;
  clMatrix mapEdgeTri;
  clMatrix IENOriginal;

  int numVerts;                   // numero total de vertices da malha
  int numVertsOriginal;
  int numElems;                   // numero total de elementos da malha
  int numNodes;                   // numero total de nos da malha
  int numTriangles;
  int numGLEU;                    // numero de graus de liberdade de vel.
  int numGLEP;                    // numero de graus de liberdade de P.
  int numGLEC;                    // numero de graus de liberdade de C.
  real rMax;                      // tamanho max do raio do disco
  real xCenter,yCenter,zCenter;
  real bubbleRadius;
  real minEdge;

  vector< list<int> > neighbourElem;  // lista de elementos de cada no
  vector< list<int> > neighbourVert;  // lista de vizinhos de cada no
  vector< list<int> > neighbourFace;  // lista de vizinhos de cada no
  vector< list<int> > faceIEN;
  vector< list<int> > neighbourFaceVert;
  vector< list<int> > elemSurface;
  vector< list<int> > surfaceViz;  // lista de vizinhos na interface
  vector< list<int> > xSurfaceViz; // lista de coords X de vizinhos na interface
  vector< list<int> > ySurfaceViz; // lista de coords Y de vizinhos na interface
  vector< list<int> > zSurfaceViz; // lista de coords Z de vizinhos na interface
  vector< list<int> > neighbourSurfaceElem; // lista de elementos triangulares
  vector< list<int> > neighbourPoint;  // lista de pontos vizinhos da superficie
  list<int> outVert,inVert; // lista de elementos do interior 
  list<int> outElem,inElem; // lista de elementos do interior

  void saveVTKSurface( const char* _dir,const char* _filename, int _iter );
};

#endif

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
  void clearBC();
  void setMeshStep(int nX,int nY,int nZ);
  void setStepBC();
  void setStepReservBC();
  void setStepReservInvBC();
  void setCouetteBC();
  void setAdimenStep();
  void setMeshDisk(int nLados1Poli,int nCircMax,int nZ);
  void mesh2Dto3D();
  void mesh2Dto3DOriginal();
  void mesh3DPoints();
  void setNuCteDiskBC();
  void setNuCDiskBC();
  void setNuZDiskBC();
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
  clVector* getOutflow();
  clVector* getIdbcu();
  clVector* getIdbcv();
  clVector* getIdbcw();
  clVector* getIdbcp();
  clVector* getIdbcc();
  clMatrix* getIEN();
  clMatrix* getIENTri();
  clMatrix* getIENConvexTri();
  int getNumVerts();
  int getNumNodes();
  int getNumElems();
  int getNumGLEU();
  int getNumGLEP();
  int getNumGLEC();
  real getVolume(int _elem);
  clMatrix* getOFace();
  clVector* getSurface();
  clVector* getNonSurface();
  real getXCenter();
  real getYCenter();
  real getZCenter();
  real getBubbleRadius();
  void setSurface();
  void setSurfaceFace();
  void setSurfaceTri();
  void setOutTri();
  void setInOutVert();
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


  clMatrix IENOriginal;
  int numVertsOriginal;
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

  int numVerts;                   // numero total de vertices da malha
  int numElems;                   // numero total de elementos da malha
  int numNodes;                   // numero total de nos da malha
  int numGLEU;                    // numero de graus de liberdade de vel.
  int numGLEP;                    // numero de graus de liberdade de P.
  int numGLEC;                    // numero de graus de liberdade de C.
  real rMax;                      // tamanho max do raio do disco
  real xCenter,yCenter,zCenter;
  real bubbleRadius;

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
  list<int> outVert,inVert; // lista de elementos do interior 
  list<int> outElem,inElem; // lista de elementos do interior
  //vector< string > idRegion; // lista de coords Z de vizinhos na interface

};

#endif

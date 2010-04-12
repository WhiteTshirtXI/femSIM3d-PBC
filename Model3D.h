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

  bool readVTK( const char* filename );
  bool readBC( const char* filename );
  void readBaseStateNu(const char* _filename);
  void clearBC();
  void setStep(int nX,int nY,int nZ);
  void setStepBC();
  void setStepReservBC();
  void setStepReservInvBC();
  void setCouetteBC();
  void setAdimenStep();
  void setDisk(int nLados1Poli,int nCircMax,int nZ);
  void setNuCteDiskBC();
  void setNuCDiskBC();
  void setNuZDiskBC();
  void setDiskFSBC();
  void setDiskCouetteBC();
  void setAdimenDiskCouette();
  void setAdimenDisk();
  void setPerturbSurf();
  void setPerturbSurf2();
  void setPerturbSurfSquare();
  void setCentroid();             // acrescenta centroide na malha
  void setNeighbour();
  void setVertNeighbour();
  void setOFace();
  bool testFace(int v1, int v2, int v3, int v4);
  //clVector sortFreeVector();
  void setBubbleBubbleBC();
  void setBubbleBC2();
  void setBubble3DBC();
  void setCubeCubeBC(real param);

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
  int getNumVerts();
  int getNumNodes();
  int getNumElems();
  int getNumGLEU();
  int getNumGLEP();
  int getNumGLEC();
  clMatrix* getOFace();
  clVector* getSurface();
  clVector* getNonSurface();
  real getXCenter();
  real getYCenter();
  real getZCenter();
  real getBubbleRadius();
  void setSurface();
  void setSurfaceFace();
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

 private:
  clVector uc,vc,wc,pc,cc;
  clMatrix IEN;
  clMatrix IENSort;
  clVector X,Y,Z;
  clVector outflow,idbcu,idbcv,idbcw,idbcp,idbcc;
  clMatrix faceFace,freeFace,mapViz;
  clMatrix oFace;
  clVector V; // vetor de volumes dos elementos de malha 
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
  list<int> outVert,inVert; // lista de vertices do interior e de contorno

};

#endif

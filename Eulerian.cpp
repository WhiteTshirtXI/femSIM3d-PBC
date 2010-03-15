// =================================================================== // 
// this is file Eulerian.cpp, created at 20-Ago-2009                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "Eulerian.h"

Eulerian::~Eulerian(){};

Eulerian::Eulerian(Model3D &_m)
{
 uc = _m.getPointerUC();
 vc = _m.getPointerVC();
 wc = _m.getPointerWC();
 pc = _m.getPointerPC();
 X = _m.getPointerX();
 Y = _m.getPointerY();
 Z = _m.getPointerZ();
 IEN = _m.getPointerIEN();
 idbcu = _m.getPointerIdbcu();
 idbcv = _m.getPointerIdbcv();
 idbcw = _m.getPointerIdbcw();
 idbcp = _m.getPointerIdbcp();
 numVerts = _m.getNumVerts();
 numNodes = _m.getNumNodes();
 numElems = _m.getNumElems();
 neighbourVert = _m.neighbourVert;
 surfaceViz = _m.surfaceViz;
 surface = _m.getPointerSurface();
 //nonSurface = _m.getPointerNonSurface();
 uSmooth.Dim(numNodes);
 vSmooth.Dim(numNodes);
 wSmooth.Dim(numNodes);
 uSmoothSurface.Dim(numNodes);
 vSmoothSurface.Dim(numNodes);
 wSmoothSurface.Dim(numNodes);
}

clVector Eulerian::compute(real _dt)
{
 dt = _dt;
 //setSurface();

 stepSmooth();
 stepSmoothTangent();
 setBC();
 setCentroid();

 clVector aux = uSmooth;
 aux.Append(vSmooth);
 aux.Append(wSmooth);
 return aux;
}

// calcula velocidade euleriana em todos os vertices
void Eulerian::stepSmooth()
{
 real aux;
 listElem plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum;
 real size; // numero de elementos da lista
 xAverage.Dim(numVerts);
 yAverage.Dim(numVerts);
 zAverage.Dim(numVerts);
 uSmooth.Dim(numNodes);
 vSmooth.Dim(numNodes);
 wSmooth.Dim(numNodes);

 // loop em todos os vertices 
 for( int i=0;i<numVerts;i++ )
 {
  plist = neighbourVert.at(i);
  size = plist.size();
  xSum = 0.0;
  ySum = 0.0;
  zSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   xSum = xSum + X->Get(*vert);
   ySum = ySum + Y->Get(*vert);
   zSum = zSum + Z->Get(*vert);
  }
  xAverage.Set( i,xSum/size ); // X medio
  aux = (xAverage.Get(i) - X->Get(i))/dt; // velocidade USmooth
  uSmooth.Set(i,aux);

  yAverage.Set( i,ySum/size ); // Y medio
  aux = (yAverage.Get(i) - Y->Get(i))/dt; // velocidade VSmooth
  vSmooth.Set(i,aux);

  zAverage.Set( i,zSum/size ); // Y medio
  aux = (zAverage.Get(i) - Z->Get(i))/dt; // velocidade WSmooth
  wSmooth.Set(i,aux);
//--------------------------------------------------
//   if( i == 0 )
//   {
//   cout << "-------- " << i << " --------" << endl;
//   cout << "uSmooth: " << uSmooth.Get(i) << endl;
//   cout << "vSmooth: " << vSmooth.Get(i) << endl;
//   cout << "wSmooth: " << wSmooth.Get(i) << endl;
//   cout << "xAverage: " << xAverage << endl;
//   cout << "yAverage: " << yAverage << endl;
//   cout << "zAverage: " << zAverage << endl;
//   cout << "X: " << X->Get(i) << endl;
//   cout << "Y: " << Y->Get(i) << endl;
//   cout << "Z: " << Z->Get(i) << endl;
//   cout << "dt: " << dt << endl;
//   cout << "--------------------" << endl;
//   }
//-------------------------------------------------- 
 }
} // fecha metodo stepSmooth

void Eulerian::setBC()
{
 int aux;
 // ---------------------------------------------------------------- //
 // impondo as condicoes de contorno para u, v e w
 // ---------------------------------------------------------------- //
 for( int i=0;i<idbcu->Dim();i++ )
 {
  aux = (int) idbcu->Get(i);
  uSmooth.Set( aux,uc->Get(aux));
 }

 for( int i=0;i<idbcv->Dim();i++ )
 {
  aux = (int) idbcv->Get(i);
  vSmooth.Set( aux,vc->Get(aux));
 }

 for( int i=0;i<idbcw->Dim();i++ )
 {
  aux = (int) idbcw->Get(i);
  wSmooth.Set( aux,wc->Get(aux));
 }
}

void Eulerian::setCentroid()
{  
 // calculando os valores de uSmooth e vSmooth nos centroides
 // atraves da media dos valores nos vertices
 int v[5];
 real aux;

 for( int mele=0;mele<numElems;mele++ )
 {
  v[0] = (int) IEN->Get(mele,0);
  v[1] = (int) IEN->Get(mele,1);
  v[2] = (int) IEN->Get(mele,2);
  v[3] = (int) IEN->Get(mele,3);
  v[4] = (int) IEN->Get(mele,4);

  aux = ( uSmooth.Get(v[0])+uSmooth.Get(v[1])+
	      uSmooth.Get(v[2])+uSmooth.Get(v[3]) )*0.25;
  uSmooth.Set( v[4],aux );
  aux = ( vSmooth.Get(v[0])+vSmooth.Get(v[1])+
	      vSmooth.Get(v[2])+vSmooth.Get(v[3]) )*0.25;
  vSmooth.Set( v[4],aux );
  aux = ( wSmooth.Get(v[0])+wSmooth.Get(v[1])+
	      wSmooth.Get(v[2])+wSmooth.Get(v[3]) )*0.25;
  wSmooth.Set( v[4],aux );
 }
}

// calcula velocidade tangencial euleriana nos vertices da interface
// atraves da media dos vizinhos.
void Eulerian::stepSmoothTangent()
{
 // velocidade da interface eh igual a velocidade do fluido
 // para isso precisa-se encontrar os vertices da interface e impor a
 // velocidade do fluido (calculada pelo Semi-lagrangeano)
 real aux,j;

 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  j = surface->Get(i);

  // ******************************************************** //
  // eh necessario pensar e refazer a estrategia!!!
  // ******************************************************** //
  real closer = search( i,xAverage.Get(j),yAverage.Get(j),zAverage.Get(j) );

  // criando variaveis para nomear pontos e calcular vetor unitario
  real P0x = X->Get( j );
  real P0y = Y->Get( j );
  real P0z = Z->Get( j );
  real P1x = X->Get( closer );
  real P1y = Y->Get( closer );
  real P1z = Z->Get( closer );

  // calculo do vetor unitario na direcao 0-1
  // x1Unit e y1Unit sao coordenadas do vetor unitario com origem em
  // 0,0. Para recoloca-lo em sua posicao original eh necessario 
  // somar P0x e P0y ao vetor unitario; 
  real dist01 = sqrt( (P1x-P0x)*(P1x-P0x)+
                      (P1y-P0y)*(P1y-P0y)+
					  (P1z-P0y)*(P1z-P0y) );
  real x1Unit = (P1x-P0x)/dist01;
  real y1Unit = (P1y-P0y)/dist01;
  real z1Unit = (P1z-P0z)/dist01;

  // velocidade suavizada ja calculada acima
  real velX = uSmooth.Get(j); // velocidade na direcao x
  real velY = vSmooth.Get(j); // velocidade na direcao y
  real velZ = wSmooth.Get(j); // velocidade na direcao y

  // calculo da projecao da velocidade v na direcao da aresta 0-1
  // proj = | v.u^ | . u^
  aux = velX*x1Unit + velY*y1Unit + velZ*z1Unit; // v.u^ (produto escalar)

  //real xProj =  aux*x1Unit; 
  //real yProj =  aux*y1Unit;
  //real zProj =  aux*z1Unit;

  uSmooth.Set(j,0.0);
  vSmooth.Set(j,0.0);
  wSmooth.Set(j,0.0);
  //uSmooth.Set(j,xProj);
  //vSmooth.Set(j,yProj);
  //wSmooth.Set(j,zProj);
 }
} // fecha metodo stepSmoothTangent

// a velocidade eh calculada tomando como base a media das posicoes dos
// vertices vizinhos pertencentes a interface. A distribuicao de pontos
// na interface fica mais uniforme que na do metodo stepSmoothTangent.
void Eulerian::stepSmoothTangent2()
{
 real j,aux;
 listElem plist;
 list<int>::iterator vert;
 real xSum,ySum;
 real size; // numero de elementos da lista


 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  j = surface->Get(i);

  listElem plist;
  list<int>::iterator vert;
  plist = surfaceViz.at(i);
  size = plist.size();
  xSum = 0.0;
  ySum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   xSum = xSum + X->Get(*vert);
   ySum = ySum + Y->Get(*vert);
  }
  xAverage.Set( j,xSum/size ); // X medio
  aux = (xAverage.Get(j) - X->Get(j))/dt; // velocidade USmooth
  uSmooth.Set(j,aux);

  yAverage.Set( j,ySum/size ); // Y medio
  aux = (yAverage.Get(j) - Y->Get(j))/dt; // velocidade VSmooth
  vSmooth.Set(j,aux);

  real closer = search( i,xAverage.Get(j),yAverage.Get(j),zAverage.Get(j) );

  // criando variaveis para nomear pontos e calcular vetor unitario
  real P0x = X->Get( j );
  real P0y = Y->Get( j );
  real P1x = X->Get( closer );
  real P1y = Y->Get( closer );

  // calculo do vetor unitario na direcao 0-1
  // x1Unit e y1Unit sao coordenadas do vetor unitario com origem em
  // 0,0. Para recoloca-lo em sua posicao original eh necessario 
  // somar P0x e P0y ao vetor unitario; 
  real dist01 = sqrt( (P1x-P0x)*(P1x-P0x)+(P1y-P0y)*(P1y-P0y) );
  real x1Unit = (P1x-P0x)/dist01;
  real y1Unit = (P1y-P0y)/dist01;

  // velocidade suavizada ja calculada acima
  real velX = uSmooth.Get(j); // velocidade na direcao x
  real velY = vSmooth.Get(j); // velocidade na direcao y

  // calculo da projecao da velocidade v na direcao da aresta 0-1
  // proj = | v.u^ | . u^
  aux = velX*x1Unit + velY*y1Unit; // v.u^ (produto escalar)
  real xProj =  aux*x1Unit; 
  real yProj =  aux*y1Unit;

  //uSmooth.Set(j,0.0);
  //vSmooth.Set(j,0.0);
  uSmooth.Set(j,xProj);
  vSmooth.Set(j,yProj);
 }
} // fecha metodo stepSmoothTangent2

void Eulerian::setSurface()
{
 int aux;

 // xSurface representa os valores de X dos nos da interface
 // ySurface representa os valores de Y dos nos da interface
 xSurface.Dim(surface->Dim());
 ySurface.Dim(surface->Dim());
 zSurface.Dim(surface->Dim());
 for( int i=0;i<surface->Dim();i++ )
 {
  aux = surface->Get(i);
  xSurface.Set(i,X->Get(aux));
  ySurface.Set(i,Y->Get(aux));
  zSurface.Set(i,Z->Get(aux));
 }
}

int Eulerian::search(int node,real _XI, real _YI, real _ZI)
{
 listElem plist;
 list<int>::iterator vert;
 real dist, dmin;
 plist = surfaceViz.at(node);

 // o primeiro elemento da lista eh a coordenada de trabalho
 vert=plist.begin();
 ++vert; // indo para o 2o. elemento da lista

 // pega o 2o. elemento da lista para comparar
 dmin = (_XI-X->Get(*vert))*(_XI-X->Get(*vert))+
        (_YI-Y->Get(*vert))*(_YI-Y->Get(*vert))+
        (_ZI-Z->Get(*vert))*(_ZI-Z->Get(*vert));

 int vmin=*vert; // atribuindo o vmin do primeiro elemento
 //cout << "vmin1: " << vmin << endl;

 ++vert; // indo para o 3o. elemento da lista
 //cout << "vmin2: " << *vert << endl;

 for( vert=vert; vert != plist.end(); ++vert )
 {
 dist = (_XI-X->Get(*vert))*(_XI-X->Get(*vert))+
        (_YI-Y->Get(*vert))*(_YI-Y->Get(*vert))+
        (_ZI-Z->Get(*vert))*(_ZI-Z->Get(*vert));

  if (dist<dmin)
  {
   dmin=dist;
   vmin=*vert; 
  }
 }
 return vmin;
}



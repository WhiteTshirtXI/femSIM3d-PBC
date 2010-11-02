// =================================================================== //
// this is file Interface3D.cpp, created at 26-Mar-2009                //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Interface3D.h"
#include "InOut.h"

Interface3D::Interface3D(Model3D &_m)
{
 m = &_m;
 cc = m->getCC();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 IEN = m->getIEN();
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 neighbourElem = m->getNeighbourElem();
 neighbourVert = m->getNeighbourVert(); 
 neighbourFace = m->getNeighbourFace(); 
 elemSurface = m->getElemSurface(); 
 neighbourFaceVert = m->getNeighbourFaceVert();
 surfaceViz = m->getSurfaceViz();
 surface = m->getSurface();
 closer.Dim(numNodes);
 xCloser.Dim(numNodes);
 yCloser.Dim(numNodes);
 zCloser.Dim(numNodes);
 closerViz.Dim(numNodes);
 distance.Dim(numVerts);
 xCenter = m->getXCenter();
 yCenter = m->getYCenter();
 zCenter = m->getZCenter();
 bubbleRadius = m->getBubbleRadius();
 setSolverSmooth(new PCGSolver());
}

Interface3D::~Interface3D()
{
 //delete solverC;
}

clVector Interface3D::curvature1()
{
 setCloser();
 real aux;

 // esta rotina retorna todos os valores em X,Y e Z dos nos da interface
 // vinculados pela menor distancia de todos os nos da malha.
 // calculo da distancia atraves da norma 2
 for( int i=0;i<numVerts;i++ )
 {
  aux = sqrt( (X->Get(i)-xCloser.Get(i))*(X->Get(i)-xCloser.Get(i))+
              (Y->Get(i)-yCloser.Get(i))*(Y->Get(i)-yCloser.Get(i))+
              (Z->Get(i)-zCloser.Get(i))*(Z->Get(i)-zCloser.Get(i)) );
  distance.Set(i,fabs(aux));
 }

 // sign distance function
 clVector ccAux = *cc;ccAux-0.5;ccAux=2*ccAux;

 // sign distance function (phi)
 distance = distance.MultVec(ccAux);

 return distance;
} // fecha metodo curvature1 

// imposicao da curvatura --> este metodo nao eh um calculo, so serve
// para testes simples do codigo
clVector Interface3D::curvature2()
{
 setCloser();
 real aux,aux1,aux2,aux3;
 clVector xAux(numNodes);
 clVector yAux(numNodes);
 clVector zAux(numNodes);

 //rotina para calculo da distancia dos nos da malha para seu
 // respectivo no da interface atraves do teorema de pitagoras
 for( int i=0;i<numVerts;i++ )
 {
  aux1 = (X->Get(i)-xCenter)*(X->Get(i)-xCenter);
  aux2 = (Y->Get(i)-yCenter)*(Y->Get(i)-yCenter);
  aux3 = (Z->Get(i)-zCenter)*(Z->Get(i)-zCenter);
  aux = sqrt( aux1+aux2+aux3 );
  aux = fabs(aux - bubbleRadius);
  distance.Set(i,aux);
 }

 // sign distance function
 clVector ccAux = *cc;ccAux-0.5;ccAux=2*ccAux;

 // sign distance function (phi)
 distance = distance.MultVec(ccAux);

 return distance;
} // fecha metodo curvature2

// calculo do Kappa geometricamente. Utiliza neighbourVert
// No calculo dos vetores da superficie (triangulo), este metodo
// encontra a metade do segmento da aresta do triangulo e depois 
// calcula a area para, em seguida, integrar a forca na aresta.
clVector Interface3D::computeKappa1()
{
 setCloser();
 int surfaceNode;
 real force,sumForce,sumArea,sumLength,fx,fy,fz;
 real xCross,yCross,zCross;
 list<int> plist,plist2;
 list<int>::iterator face,vert;

 // loop sobre todos os nos da superficie 
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  real P0x = X->Get(surfaceNode);
  real P0y = Y->Get(surfaceNode);
  real P0z = Z->Get(surfaceNode);

  //int c1 = 0;
  xCross = 0;
  yCross = 0;
  zCross = 0;
  fx = 0;
  fy = 0;
  fz = 0;
  force = 0;
  sumForce = 0;
  sumArea = 0;
  sumLength = 0;

  // loop sobre os vizinhos do vertice i
  plist = elemSurface->at (surfaceNode); 
  for( face=plist.begin();face!=plist.end();++face )
  {
   // 3D: 2 pontos pois a face em 3D pertencente a superficie contem 
   // 3 pontos (P0 - surfaceNode, P1 e P2 que sao pontos do triangulo)
   plist2 = neighbourFaceVert->at (*face);
   vert=plist2.begin();
   int v1 = *vert;++vert;
   int v2 = *vert;
   vert=plist2.end();

   real P1x = X->Get(v1);
   real P1y = Y->Get(v1);
   real P1z = Z->Get(v1);
   real P2x = X->Get(v2);
   real P2y = Y->Get(v2);
   real P2z = Z->Get(v2);

   // ponto medio aresta 01
   real Pm01x = P0x + (P1x-P0x)/2.0; 
   real Pm01y = P0y + (P1y-P0y)/2.0; 
   real Pm01z = P0z + (P1z-P0z)/2.0; 

   // ponto medio aresta 02
   real Pm02x = P0x + (P2x-P0x)/2.0; 
   real Pm02y = P0y + (P2y-P0y)/2.0; 
   real Pm02z = P0z + (P2z-P0z)/2.0; 

   // distance do ponto 0 ate metade do segmento 01
   real a = sqrt( (P0x-Pm01x)*(P0x-Pm01x)+
	              (P0y-Pm01y)*(P0y-Pm01y)+
				  (P0z-Pm01z)*(P0z-Pm01z) );

   // distance do ponto 0 ate metade do segmento 02
   real b = sqrt( (P0x-Pm02x)*(P0x-Pm02x)+
	              (P0y-Pm02y)*(P0y-Pm02y)+
			      (P0z-Pm02z)*(P0z-Pm02z) );

   // distance da metade do segmento 01 ate metade do segmento 02
   real c = sqrt( (Pm02x-Pm01x)*(Pm02x-Pm01x)+
	              (Pm02y-Pm01y)*(Pm02y-Pm01y)+
			      (Pm02z-Pm01z)*(Pm02z-Pm01z) );

   // vetores 
   real x1 = Pm01x-P0x;
   real y1 = Pm01y-P0y;
   real z1 = Pm01z-P0z;

   real x2 = Pm02x-P0x;
   real y2 = Pm02y-P0y;
   real z2 = Pm02z-P0z;

   real xReta = Pm02x-Pm01x;
   real yReta = Pm02y-Pm01y;
   real zReta = Pm02z-Pm01z;

   // vetores unitarios deslocados para origem do sistema (0,0,0)
   real x1Unit = x1/a;
   real y1Unit = y1/a;
   real z1Unit = z1/a;

   real x2Unit = x2/b;
   real y2Unit = y2/b;
   real z2Unit = z2/b;

   real xRetaUnit = xReta/c;
   real yRetaUnit = yReta/c;
   real zRetaUnit = zReta/c;

   // calculando o produto vetorial de cada elemento na superficie
   clVector cross = crossProd(x1,y1,z1,x2,y2,z2);
   xCross += cross.Get(0);
   yCross += cross.Get(1);
   zCross += cross.Get(2);

   // area do elemento triangular com vertices ( P0,P1m,P2m )
   real area = 0.5*sqrt( cross.Get(0)*cross.Get(0)+
                         cross.Get(1)*cross.Get(1)+
                         cross.Get(2)*cross.Get(2));

   // soma dos vetores 1Unit + 2Unit = resultante
   real xRes = x1Unit+x2Unit;
   real yRes = y1Unit+y2Unit;
   real zRes = z1Unit+z2Unit;

   // produto escalar --> projecao do vetor Unit no segmento de reta
   // | Unit.RetaUnit | . RetaUnit
   // resultado = vetor tangente a reta situado na superficie
   real prod = xRes*xRetaUnit + yRes*yRetaUnit + zRes*zRetaUnit;
   real xTang = xRetaUnit*prod;
   real yTang = yRetaUnit*prod;
   real zTang = zRetaUnit*prod;

   // subtraindo vetor tangente do vetor unitario para encontrar as
   // coordenadas do vetor normal ao lado 'c' situadas na superficie
   real xNormal = xRes - xTang;
   real yNormal = yRes - yTang;
   real zNormal = zRes - zTang;

   // tamanho do vetor considerando que tem origem em 0,0,0
   real len = sqrt( (xNormal*xNormal)+(yNormal*yNormal)+(zNormal*zNormal) ); 

   // Unitario do vetor resultante do plano do triangulo
   // combinacao linear dos vetores unitarios das arestas do triangulo
   real xNormalUnit = xNormal/len;
   real yNormalUnit = yNormal/len;
   real zNormalUnit = zNormal/len;

   // normal integrada na distancia (MOD) dos 2 vertices medianos
   // force = resultante das componentes * tamanho da aresta que sera
   // usada como referencia no calculo da area do triangulo
   fx += xNormalUnit*c;
   fy += yNormalUnit*c;
   fz += zNormalUnit*c;

   sumArea += area;
   sumLength += c;
   force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) ); // perdendo a direcao!

  //cout << surfaceNode << endl;
//--------------------------------------------------
//    if( surfaceNode == 29 )
//    {
// 	cout << "Triangulo: ------------------------- " << c1 << endl;
// 	cout << "v1 = " << v1 << " " << "v2 = " << v2 << endl;
// 	cout << "Unit = " << xRes << " " << yRes << " " << zRes << endl;
// 	cout << "RetaUnit = " << xRetaUnit << " " 
// 	                      << yRetaUnit << " " 
// 						  << zRetaUnit << endl;
// 	cout << "c = " << c << endl;
// 	cout << "Cross  = " << cross.Get(0) << " " 
// 	                    << cross.Get(1) << " " 
// 						<< cross.Get(2) << endl;
// 	cout << "xCrossUnit = " << xCrossUnit << endl;
// 	cout << "yCrossUnit = " << yCrossUnit << endl;
// 	cout << "zCrossUnit = " << zCrossUnit << endl;
// 	cout << "Normal = " << xNormal << " " 
// 	                    << yNormal << " " 
// 						<< zNormal << endl;
// 	cout << "area = " << area << endl;
// 	cout << "force = " << force << endl;
// 	cout << "sumArea = " << sumArea << endl;
// 	cout << "fx = " << fx << endl;
// 	cout << "fy = " << fy << endl;
// 	cout << "fz = " << fz << endl;
// 	//cout << "sumForce = " << sumForce << endl;
// 	//cout << "pressure = " << sumForce/sumArea << endl;
// 	c1++;
//    }
//-------------------------------------------------- 
   }

  // aplicando o teste para saber o sentido correto de aplicacao da
  // pressao no noh.
  if( (fx*xCross+fy*yCross+fz*zCross) < 0.0 )
   force = -force;

  real pressure = force/sumArea;
//--------------------------------------------------
// 
//   if( surfaceNode == 214 )
//   {
//   cout << surfaceNode << endl;
//   cout << "  " << "force = " << force << endl; 
//   cout << "  " << "fx = " << fx << endl; 
//   cout << "  " << "fy = " << fy << endl;
//   cout << "  " << "fz = " << fz << endl;
//   cout << "  " << "sumArea = " << sumArea << endl;
//   cout << "  " << "pressure = " << pressure << endl;
//   }
//-------------------------------------------------- 

  distance.Set( surfaceNode,pressure );
 }
 return distance;

} // fecha metodo computeKappa1

// calculo do Kappa geometricamente. Utiliza neighbourVert
// No calculo dos vetores da superficie (triangulo), este metodo
// encontra a metade da area do triangulo e depois o comprimento da base
// correspondente para integrar a forca na aresta.
clVector Interface3D::computeKappa2()
{
 setCloser();
 int surfaceNode;
 real force,sumForce,sumArea,sumLength,fx,fy,fz;
 real xCross,yCross,zCross;
 list<int> plist,plist2;
 list<int>::iterator face,vert;

 // loop sobre todos os nos da superficie 
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  real P0x = X->Get(surfaceNode);
  real P0y = Y->Get(surfaceNode);
  real P0z = Z->Get(surfaceNode);

  int c1 = 0;
  xCross = 0;
  yCross = 0;
  zCross = 0;
  fx = 0;
  fy = 0;
  fz = 0;
  force = 0;
  sumForce = 0;
  sumArea = 0;
  sumLength = 0;

  // loop sobre nos triangulos da interface/superficie vizinhos do vertice i
  plist = elemSurface->at (surfaceNode); 
  for( face=plist.begin();face!=plist.end();++face )
  {
   // 3D: 2 pontos pois a face em 3D pertencente a superficie contem 
   // 3 pontos (P0 - surfaceNode, P1 e P2 que sao pontos do triangulo)
   plist2 = neighbourFaceVert->at (*face);
   vert=plist2.begin();
   int v1 = *vert;++vert;
   int v2 = *vert;
   vert=plist2.end();

   real P1x = X->Get(v1);
   real P1y = Y->Get(v1);
   real P1z = Z->Get(v1);
   real P2x = X->Get(v2);
   real P2y = Y->Get(v2);
   real P2z = Z->Get(v2);

   // distance do ponto 0 ate metade do segmento 01
   real a = sqrt( (P0x-P1x)*(P0x-P1x)+
	              (P0y-P1y)*(P0y-P1y)+
				  (P0z-P1z)*(P0z-P1z) );

   // distance do ponto 0 ate metade do segmento 02
   real b = sqrt( (P0x-P2x)*(P0x-P2x)+
	              (P0y-P2y)*(P0y-P2y)+
			      (P0z-P2z)*(P0z-P2z) );

   // distance da metade do segmento 01 ate metade do segmento 02
   real c = sqrt( (P2x-P1x)*(P2x-P1x)+
	              (P2y-P1y)*(P2y-P1y)+
			      (P2z-P1z)*(P2z-P1z) );

   // vetores
   real x1 = P1x-P0x;
   real y1 = P1y-P0y;
   real z1 = P1z-P0z;
   real x2 = P2x-P0x;
   real y2 = P2y-P0y;
   real z2 = P2z-P0z;
   real xReta = P2x-P1x;
   real yReta = P2y-P1y;
   real zReta = P2z-P1z;

   // vetores unitarios
   real x1Unit = x1/a;
   real y1Unit = y1/a;
   real z1Unit = z1/a;

   real x2Unit = x2/b;
   real y2Unit = y2/b;
   real z2Unit = z2/b;

   real xRetaUnit = xReta/c;
   real yRetaUnit = yReta/c;
   real zRetaUnit = zReta/c;

   // calculando o produto vetorial de cada elemento triangular da superficie
   clVector cross = crossProd(x1,y1,z1,x2,y2,z2);
   xCross += cross.Get(0);
   yCross += cross.Get(1);
   zCross += cross.Get(2);

   real area = 0.5*sqrt( cross.Get(0)*cross.Get(0)+
                         cross.Get(1)*cross.Get(1)+
                         cross.Get(2)*cross.Get(2));

   // soma dos vetores 1Unit + 2Unit = resultante
   real xRes = x1Unit+x2Unit;
   real yRes = y1Unit+y2Unit;
   real zRes = z1Unit+z2Unit;

   // produto escalar --> projecao do vetor Res no segmento de reta
   // | Unit.RetaUnit | . RetaUnit
   // resultado = vetor tangente a reta situado na superficie
   real prod = xRes*xRetaUnit + yRes*yRetaUnit + zRes*zRetaUnit;
   real xTang = xRetaUnit*prod;
   real yTang = yRetaUnit*prod;
   real zTang = zRetaUnit*prod;

   // subtraindo vetor tangente do vetor unitario para encontrar as
   // coordenadas do vetor normal situada na superficie
   real xNormal = xRes - xTang;
   real yNormal = yRes - yTang;
   real zNormal = zRes - zTang;

   real len = sqrt( (xNormal*xNormal)+(yNormal*yNormal)+(zNormal*zNormal) ); 

   // Unitario do vetor resultante situado no plano do triangulo
   // combinacao linear dos vetores unitarios das arestas do triangulo
   real xNormalUnit = xNormal/len;
   real yNormalUnit = yNormal/len;
   real zNormalUnit = zNormal/len;

   // tamanho do comprimento da aresta do triangulo que tem a metade da
   // area do elemento
   // semelhanca de triangulos: area = bh/2
   // c'h' = A    ----> b'b' = Ac/h
   real h = 2*area/c;
   real e = sqrt( c*area/h );

   // normal integrada na distancia (MOD) dos 2 vertices medianos
   // force = resultante das componentes * tamanho da aresta que sera
   // usada como referencia no calculo da area do triangulo
   fx += xNormalUnit*e;
   fy += yNormalUnit*e;
   fz += zNormalUnit*e;

   sumArea += area/2.0;
   sumLength += c;
   force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) );

//--------------------------------------------------
//    if( surfaceNode == 531 || surfaceNode == 528 || surfaceNode == 525 )
//    {
// 	cout << "Triangulo: ------------------------- " << c1 << endl;
// 	cout << "v1 = " << v1 << " " << "v2 = " << v2 << endl;
// 	//cout << "Unit = " << xUnit << " " << yUnit << " " << zUnit << endl;
// 	cout << "RetaUnit = " << xRetaUnit << " " 
// 	                      << yRetaUnit << " " 
// 						  << zRetaUnit << endl;
// 	cout << "c = " << c << endl;
// 	cout << "Cross  = " << cross.Get(0) << " " 
// 	                    << cross.Get(1) << " " 
// 						<< cross.Get(2) << endl;
// 	cout << "xCross = " << xCross << endl;
// 	cout << "yCross = " << yCross << endl;
// 	cout << "zCross = " << zCross << endl;
// 	cout << "Normal = " << xNormal << " " 
// 	                    << yNormal << " " 
// 						<< zNormal << endl;
// 	cout << "area = " << area << endl;
// 	cout << "force = " << force << endl;
// 	cout << "sumArea = " << sumArea << endl;
// 	cout << "fx = " << fx << endl;
// 	cout << "fy = " << fy << endl;
// 	cout << "fz = " << fz << endl;
// 	//cout << "sumForce = " << sumForce << endl;
// 	//cout << "pressure = " << sumForce/sumArea << endl;
// 	c1++;
//    }
//-------------------------------------------------- 
  }
  // aplicando o teste para saber o sentido correto de aplicacao da
  // pressao no noh.
  if( (fx*xCross+fy*yCross+fz*zCross) < 0.0 )
   force = -force;

  real pressure = force/sumArea;

//--------------------------------------------------
//   if( surfaceNode == 1268 )
//   {
//   cout << surfaceNode << endl;
//   cout << "  " << "force = " << force << endl; 
//   cout << "  " << "fx = " << fx << endl; 
//   cout << "  " << "fy = " << fy << endl;
//   cout << "  " << "fz = " << fz << endl;
//   cout << "  " << "sumArea = " << sumArea << endl;
//   cout << "  " << "pressure = " << pressure << endl;
//   }
//-------------------------------------------------- 

  distance.Set( surfaceNode,pressure );
 }
 return distance;

} // fecha metodo computeKappa2

clVector Interface3D::smoothing(clMatrix &_AcTilde,clVector &_b1cTilde)
{
 clVector cTilde(_b1cTilde.Dim());
 cout << " ----> calculando escalar ----------- " << endl;
 solverC->solve(1E-08,_AcTilde,cTilde,_b1cTilde);
 cout << " ------------------------------------ " << endl;
 return cTilde;

}

void Interface3D::setCloser()
{
 int aux;

 // xSurface representa os valores de X dos nos da interface
 // ySurface representa os valores de Y dos nos da interface
 // zSurface representa os valores de Z dos nos da interface
 xSurface.Dim( surface->Dim() );
 ySurface.Dim( surface->Dim() );
 zSurface.Dim( surface->Dim() );
 for( int i=0;i<surface->Dim();i++ )
 {
  aux = surface->Get(i);
  xSurface.Set(i,X->Get( aux ));
  ySurface.Set(i,Y->Get( aux ));
  zSurface.Set(i,Z->Get( aux ));
 }

 // closer=surface(dsearchn(X(surface),Y(surface),X,Y));
 // esta funcao retorna o noh da interface (surface) mais 
 // proximo de cada noh da malha (vertices)
 closer = dsearchn(xSurface,ySurface,zSurface,*X,*Y,*Z);
 for( int i=0;i<closer.Dim();i++ )
 {
  aux = closer.Get(i);
  closer.Set(i,surface->Get(aux)); // alterando os valores de closer(i)
  xCloser.Set(i,X->Get( closer.Get(i) ));
  yCloser.Set(i,Y->Get( closer.Get(i) ));
  zCloser.Set(i,Z->Get( closer.Get(i) ));
 }
}

// espalhando capa calculado na superfice para todos os pontos
clDMatrix Interface3D::setKappaSurface(clVector &_kappaAux)
{
 int aux;
 clDMatrix kappa(3*numNodes);
 for( int i=0;i<numNodes;i++ )
 {
  aux = closer.Get(i);
  kappa.Set(i,_kappaAux.Get(aux));
  kappa.Set(i+numNodes,_kappaAux.Get(aux));
  kappa.Set(i+2*numNodes,_kappaAux.Get(aux));
 }
 return kappa;
}

// espalhando capa calculado na superfice para todos os pontos
clDMatrix Interface3D::setKappaSurface(clVector &_kappaNx,
                                       clVector &_kappaNy,
									   clVector &_kappaNz)
{
 int aux;
 clDMatrix kappa(3*numNodes);
 for( int i=0;i<numNodes;i++ )
 {
  aux = closer.Get(i);
  kappa.Set(i,_kappaNx.Get(aux));
  kappa.Set(i+numNodes,_kappaNy.Get(aux));
  kappa.Set(i+2*numNodes,_kappaNz.Get(aux));
 }
 return kappa;
}

clVector Interface3D::dsearchn(clVector _X,clVector _Y,clVector _Z,
                               clVector &_XI,clVector &_YI,clVector &_ZI)
{
 clVector vmin(_XI.Dim());
 vmin.SetAll(0);
 real dist, dmin;

 // loop sobre todos os pontos a serem procurados
 for( int i=0;i<_XI.Dim();i++ )
 {

  // distancia de um ponto qualquer x(0),y(o) para o ponto informado xi,yi
  dmin = (_XI.Get(i)-_X.Get(0))*(_XI.Get(i)-_X.Get(0))+
         (_YI.Get(i)-_Y.Get(0))*(_YI.Get(i)-_Y.Get(0))+
         (_ZI.Get(i)-_Z.Get(0))*(_ZI.Get(i)-_Z.Get(0));

  // loop sobre todos os pontos da malha (menos o ponto inicial) para
  // saber qual o ponto apresenta a menor distancia com relacao ao ponto
  // informado xi,yi
  for (int j=1;j<_X.Dim();j++)
  {
   //calculo da distancia de outro ponto da malha x(i),y(i)
   dist = (_XI.Get(i)-_X.Get(j))*(_XI.Get(i)-_X.Get(j))+
	      (_YI.Get(i)-_Y.Get(j))*(_YI.Get(i)-_Y.Get(j))+
	      (_ZI.Get(i)-_Z.Get(j))*(_ZI.Get(i)-_Z.Get(j));
   if (dist<dmin)
   {
	dmin=dist;
	vmin.Set(i,j); 
   }
  }
 }

 // retorna o indice do ponto da malha
 return vmin;
}

void Interface3D::setSolverSmooth(Solver *s){ solverC = s; }

clVector Interface3D::getCloser(){ return closer; }

clVector Interface3D::crossProd(real _x1,real _y1,real _z1,
                                real _x2,real _y2,real _z2)
{
 real crossX = (_y1*_z2)-(_z1*_y2);
 real crossY = -( (_x1*_z2)-(_z1*_x2) );
 real crossZ = (_x1*_y2)-(_y1*_x2);

 clVector cross(3);
 cross.Set(0,crossX);
 cross.Set(1,crossY);
 cross.Set(2,crossZ);

 return cross;
}

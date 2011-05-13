// =================================================================== // 
// this is file SemiLagrangean.cpp, created at 20-Ago-2007             //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "SemiLagrangean.h"

SemiLagrangean::SemiLagrangean( Model3D &_m )
{
 m = &_m;
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 IEN = m->getIEN();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 convLin.Dim(numVerts,numVerts);
 neighbourElem = m->getNeighbourElem(); 
 oFace = m->getOFace();
}

SemiLagrangean::SemiLagrangean(Model3D &_m,clVector &_uSol,
                                           clVector &_vSol,
										   clVector &_wSol,
										   clVector &_velU,
										   clVector &_velV,
										   clVector &_velW,
										   clVector &_cSol)
{
 m = &_m;
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 IEN = m->getIEN();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 uSol = _uSol;
 vSol = _vSol;
 wSol = _wSol;
 velU = _velU;
 velV = _velV;
 velW = _velW;
 cSol = _cSol;
 convLin.Dim(numVerts,numVerts);
 neighbourElem = m->getNeighbourElem(); 
 oFace = m->getOFace();
}

void SemiLagrangean::compute(real dt)
{
 getDepartElem(dt); // procura de elemento

 clVector uVert(numVerts);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 uSol.CopyTo(0,uVert);
 vSol.CopyTo(0,vVert);
 wSol.CopyTo(0,wVert);

 uParticle = convLin*uVert; 
 vParticle = convLin*vVert; 
 wParticle = convLin*wVert;
 cParticle = convLin*cSol;

 clVector zerosConv(numNodes-numVerts);
 uParticle.Append(zerosConv);
 vParticle.Append(zerosConv);
 wParticle.Append(zerosConv);

 setBC();
 setCentroid();

} // fim do metodo compute

void SemiLagrangean::computeFreeSurface(real dt)
{
 getDepartElem2(dt); // procura de elemento

 clVector uVert(numVerts);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 uSol.CopyTo(0,uVert);
 vSol.CopyTo(0,vVert);
 wSol.CopyTo(0,wVert);

 uParticle = convLin*uVert; 
 vParticle = convLin*vVert; 
 wParticle = convLin*wVert;
 cParticle = convLin*cSol;

 clVector zerosConv(numNodes-numVerts);
 uParticle.Append(zerosConv);
 vParticle.Append(zerosConv);
 wParticle.Append(zerosConv);

 setBC();
 setCentroid();

} // fim do metodo computeFreeSurface

void SemiLagrangean::setCentroid()
{
 int v[5];
 real aux;

 for( int mele=0;mele<numElems;mele++ )
 {
  v[0] = (int) IEN->Get(mele,0);
  v[1] = (int) IEN->Get(mele,1);
  v[2] = (int) IEN->Get(mele,2);
  v[3] = (int) IEN->Get(mele,3);
  v[4] = (int) IEN->Get(mele,4);

  aux = ( uParticle.Get(v[0])+
	      uParticle.Get(v[1])+
		  uParticle.Get(v[2])+
		  uParticle.Get(v[3]) )/4.0;
  uParticle.Set(v[4],aux);

  aux = ( vParticle.Get(v[0])+
          vParticle.Get(v[1])+
          vParticle.Get(v[2])+
	 	  vParticle.Get(v[3]) )/4.0;
  vParticle.Set(v[4],aux);

  aux = ( wParticle.Get(v[0])+
	      wParticle.Get(v[1])+
          wParticle.Get(v[2])+
		  wParticle.Get(v[3]) )/4.0;
  wParticle.Set(v[4],aux);
 }

}// fim do metodo compute -> setCentroid

void SemiLagrangean::getDepartElem(real dt)
{
 real xP,yP,zP;

 clVector xParticle = *X - velU*dt; 
 clVector yParticle = *Y - velV*dt;
 clVector zParticle = *Z - velW*dt;

 list<int> plist;
 list<int>::iterator mele;

 for (int ii = 0; ii < numVerts; ii++)
 {
  plist = neighbourElem->at(ii);
  mele=plist.begin(); // pega o primeiro elemento da lista, seja la qual for

  xP = xParticle.Get(ii);
  yP = yParticle.Get(ii);
  zP = zParticle.Get(ii);

  jumpToElem(*mele,ii,xP,yP,zP);
 } 
} // fim do metodo compute -> getDepartElem

// rotina usada em simulacao FREE SURFACE
void SemiLagrangean::getDepartElem2(real dt)
{
 real xP,yP,zP;

 //clVector xParticle = *X; 
 //clVector yParticle = *Y - vSol*dt;
 //clVector zParticle = *Z - wSol*dt;

 //clVector xParticle = *X - uSol*dt; 
 //clVector yParticle = *Y;
 //clVector zParticle = *Z - wSol*dt;

 clVector xParticle = *X - uSol*dt; 
 clVector yParticle = *Y - vSol*dt;
 clVector zParticle = *Z;

 list<int> plist;
 list<int>::iterator mele;

 for (int ii = 0; ii < numVerts; ii++)
 {
  plist = neighbourElem->at(ii);
  mele=plist.begin(); // pega o primeiro elemento da lista, seja la qual for

  xP = xParticle.Get(ii);
  yP = yParticle.Get(ii);
  zP = zParticle.Get(ii);

  //cout << "vertice de origem = " << ii << endl;
  jumpToElem(*mele,ii,xP,yP,zP);
 } 
} // fim do metodo compute -> getDepartElem2


// verifica se esta dentro do elemento, se estiver calcula os
// coeficientes de interpolacao (convLin). Se nao encontrar e tiver
// vizinho, pula para o novo vizinho. Se nao, interpola na face
void SemiLagrangean::jumpToElem(int destElem,int iiVert,real R2X,real R2Y,real R2Z)
{
 real l1,l2,l3,l4;
 real Bl1,Bl2,Bl3;
 int v[4],v1,v2,v3,v4,vjump;
 int ib1=0;
 int ib2=0;
 int ib3=0;
//--------------------------------------------------
//  cout << iiVert << " " << R2X << " " 
//                        << R2Y << " " 
// 					   << R2Z << endl;
//-------------------------------------------------- 
 if( testElement(destElem,iiVert,R2X,R2Y,R2Z,&l1,&l2,&l3,&l4) )
 {
  v[0] = v1 = (int) IEN->Get(destElem,0);
  v[1] = v2 = (int) IEN->Get(destElem,1);
  v[2] = v3 = (int) IEN->Get(destElem,2);
  v[3] = v4 = (int) IEN->Get(destElem,3);
  convLin.Set(iiVert,v[0],l1);
  convLin.Set(iiVert,v[1],l2);
  convLin.Set(iiVert,v[2],l3);
  convLin.Set(iiVert,v[3],l4);

  //cout << " Ponto localizado " << iiVert << endl;
 } 
 else 
 {
  v[0] = v1 = (int) IEN->Get(destElem,0);
  v[1] = v2 = (int) IEN->Get(destElem,1);
  v[2] = v3 = (int) IEN->Get(destElem,2);
  v[3] = v4 = (int) IEN->Get(destElem,3);

  if((l1<=l2) && (l1<=l3) && (l1<=l4)){ vjump=0; ib1=v2; ib2=v3; ib3=v4; };
  if((l2<=l1) && (l2<=l3) && (l2<=l4)){ vjump=1; ib1=v1; ib2=v3; ib3=v4; };
  if((l3<=l1) && (l3<=l2) && (l3<=l4)){ vjump=2; ib1=v1; ib2=v2; ib3=v4; };
  if((l4<=l1) && (l4<=l2) && (l4<=l3)){ vjump=3; ib1=v1; ib2=v2; ib3=v3; };

  if ( oFace->Get(destElem,vjump)!=-1)
  {
   jumpToElem( (int) oFace->Get(destElem,vjump), iiVert,R2X,R2Y,R2Z );
  }
  else
  {
   computeIntercept( iiVert,R2X,R2Y,R2Z,ib1,ib2,ib3,&Bl1,&Bl2,&Bl3 );
   convLin.Set( iiVert,ib1,Bl1 );
   convLin.Set( iiVert,ib2,Bl2 );  
   convLin.Set( iiVert,ib3,Bl3 );  
  }
 }
 return;
}

// calcula ponto de intersecao resolvendo um sistema linear 3x3 atraves
// da regra de Cramer
//
// | a1 b1 c1 |   | x1 | = | d1 |
// | a2 b2 c2 | . | x2 | = | d2 |
// | a3 b3 c3 |   | x3 | = | d3 |
//
//       detx          dety
// x1 = ------   x2 = ------   x3 = 1.0-x1-x2
//       det            det
void SemiLagrangean::computeIntercept(int ii,real R2X,real R2Y,real R2Z,int ib1,int ib2,int ib3,real *Bl1,real *Bl2,real *Bl3)
{
 real R1X = X->Get(ii); real R1Y = Y->Get(ii); real R1Z = Z->Get(ii);

 real B1X = X->Get(ib1); real B1Y = Y->Get(ib1); real B1Z = Z->Get(ib1);
 real B2X = X->Get(ib2); real B2Y = Y->Get(ib2); real B2Z = Z->Get(ib2);
 real B3X = X->Get(ib3); real B3Y = Y->Get(ib3); real B3Z = Z->Get(ib3);

 real a1 = B1X-B3X; real b1 = B2X-B3X; real c1 = R1X-R2X; real d1 = R1X-B3X;
 real a2 = B1Y-B3Y; real b2 = B2Y-B3Y; real c2 = R1Y-R2Y; real d2 = R1Y-B3Y;
 real a3 = B1Z-B3Z; real b3 = B2Z-B3Z; real c3 = R1Z-R2Z; real d3 = R1Z-B3Z;

 real det  = (a1*b2*c3)+(a3*b1*c2)+(a2*b3*c1)-
             (a3*b2*c1)-(a1*b3*c2)-(a2*b1*c3);
 real detx = (d1*b2*c3)+(d3*b1*c2)+(d2*b3*c1)-
             (d3*b2*c1)-(d1*b3*c2)-(d2*b1*c3);
 real dety = (a1*d2*c3)+(a3*d1*c2)+(a2*d3*c1)-
             (a3*d2*c1)-(a1*d3*c2)-(a2*d1*c3);

 real x1 = detx/det;
 real x2 = dety/det;

 *Bl1 = x1;
 *Bl2 = x2;
 *Bl3 = 1.0-*Bl1-*Bl2;

 return;
}

bool SemiLagrangean::testElement(int mele,int ii,real xP,real yP,real zP, real *l1,real *l2,real *l3,real *l4)
{
 int v[4],v1,v2,v3,v4;
 real V,V1,V2,V3;
 real EPSlocal = 10e-5;

 v[0] = v1 = (int) IEN->Get(mele,0);
 v[1] = v2 = (int) IEN->Get(mele,1);
 v[2] = v3 = (int) IEN->Get(mele,2);
 v[3] = v4 = (int) IEN->Get(mele,3);

 V = (1.0/6.0) * (+1*( (X->Get(v2)*Y->Get(v3)*Z->Get(v4)) 
	                    +(Y->Get(v2)*Z->Get(v3)*X->Get(v4)) 
	                    +(Z->Get(v2)*X->Get(v3)*Y->Get(v4)) 
	                    -(Y->Get(v2)*X->Get(v3)*Z->Get(v4)) 
	                    -(X->Get(v2)*Z->Get(v3)*Y->Get(v4)) 
	                    -(Z->Get(v2)*Y->Get(v3)*X->Get(v4)) )
	        -X->Get(v1)*( +Y->Get(v3)*Z->Get(v4)
		                 +Y->Get(v2)*Z->Get(v3) 
		                 +Z->Get(v2)*Y->Get(v4)
		                 -Y->Get(v2)*Z->Get(v4)
		 	  		     -Z->Get(v3)*Y->Get(v4) 
			 		     -Z->Get(v2)*Y->Get(v3) )
	        +Y->Get(v1)*( +X->Get(v3)*Z->Get(v4)
	 	 			     +X->Get(v2)*Z->Get(v3)
	 				     +Z->Get(v2)*X->Get(v4)
		                 -X->Get(v2)*Z->Get(v4)
					     -Z->Get(v3)*X->Get(v4) 
					     -Z->Get(v2)*X->Get(v3) )
		    -Z->Get(v1)*( +X->Get(v3)*Y->Get(v4)
			             +X->Get(v2)*Y->Get(v3) 
					     +Y->Get(v2)*X->Get(v4)
		                 -X->Get(v2)*Y->Get(v4)
				         -Y->Get(v3)*X->Get(v4) 
					     -Y->Get(v2)*X->Get(v3) ) );

 V1 = (1.0/6.0) * (+1*( (X->Get(v2)*Y->Get(v3)*Z->Get(v4)) 
	                     +(Y->Get(v2)*Z->Get(v3)*X->Get(v4)) 
	                     +(Z->Get(v2)*X->Get(v3)*Y->Get(v4)) 
	                     -(Y->Get(v2)*X->Get(v3)*Z->Get(v4)) 
	                     -(X->Get(v2)*Z->Get(v3)*Y->Get(v4)) 
	                     -(Z->Get(v2)*Y->Get(v3)*X->Get(v4)) )
               -xP*( +Y->Get(v3)*Z->Get(v4)
		                  +Y->Get(v2)*Z->Get(v3) 
			              +Z->Get(v2)*Y->Get(v4)
		                  -Y->Get(v2)*Z->Get(v4)
				   	      -Z->Get(v3)*Y->Get(v4) 
					      -Z->Get(v2)*Y->Get(v3) )
	           +yP*( +X->Get(v3)*Z->Get(v4)
			 		      +X->Get(v2)*Z->Get(v3)
					      +Z->Get(v2)*X->Get(v4)
		                  -X->Get(v2)*Z->Get(v4)
					      -Z->Get(v3)*X->Get(v4) 
					      -Z->Get(v2)*X->Get(v3) )
		       -zP*( +X->Get(v3)*Y->Get(v4)
			              +X->Get(v2)*Y->Get(v3) 
					      +Y->Get(v2)*X->Get(v4)
		                  -X->Get(v2)*Y->Get(v4)
				          -Y->Get(v3)*X->Get(v4) 
					      -Y->Get(v2)*X->Get(v3) ) );

 V2 = (1.0/6.0) * (+1*( (xP*Y->Get(v3)*Z->Get(v4)) 
	                    +(yP*Z->Get(v3)*X->Get(v4)) 
	                    +(zP*X->Get(v3)*Y->Get(v4)) 
	                    -(yP*X->Get(v3)*Z->Get(v4)) 
	                    -(xP*Z->Get(v3)*Y->Get(v4)) 
	                    -(zP*Y->Get(v3)*X->Get(v4)) )
	                -X->Get(v1)*( +Y->Get(v3)*Z->Get(v4)
		                 +yP*Z->Get(v3) 
			             +zP*Y->Get(v4)
		                 -yP*Z->Get(v4)
					             -Z->Get(v3)*Y->Get(v4) 
					     -zP*Y->Get(v3) )
	                +Y->Get(v1)*( +X->Get(v3)*Z->Get(v4)
					     +xP*Z->Get(v3)
					     +zP*X->Get(v4)
		                 -xP*Z->Get(v4)
					             -Z->Get(v3)*X->Get(v4) 
					     -zP*X->Get(v3) )
		            -Z->Get(v1)*( +X->Get(v3)*Y->Get(v4)
			             +xP*Y->Get(v3) 
					     +yP*X->Get(v4)
		                 -xP*Y->Get(v4)
				                 -Y->Get(v3)*X->Get(v4) 
					     -yP*X->Get(v3) ) );

 V3 = (1.0/6.0) * (+1*( (X->Get(v2)*yP*Z->Get(v4)) 
	                    +(Y->Get(v2)*zP*X->Get(v4)) 
	                    +(Z->Get(v2)*xP*Y->Get(v4)) 
	                    -(Y->Get(v2)*xP*Z->Get(v4)) 
	                    -(X->Get(v2)*zP*Y->Get(v4)) 
	                    -(Z->Get(v2)*yP*X->Get(v4)) )
	        -X->Get(v1)*( +yP*Z->Get(v4)
	                     +Y->Get(v2)*zP 
	                     +Z->Get(v2)*Y->Get(v4)
	                     -Y->Get(v2)*Z->Get(v4)
						 -zP*Y->Get(v4) 
						 -Z->Get(v2)*yP )
			+Y->Get(v1)*( +xP*Z->Get(v4)
			             +X->Get(v2)*zP
						 +Z->Get(v2)*X->Get(v4)
						 -X->Get(v2)*Z->Get(v4)
						 -zP*X->Get(v4) 
						 -Z->Get(v2)*xP )
			-Z->Get(v1)*( +xP*Y->Get(v4)
			             +X->Get(v2)*yP 
						 +Y->Get(v2)*X->Get(v4)
						 -X->Get(v2)*Y->Get(v4)
						 -yP*X->Get(v4) 
						 -Y->Get(v2)*xP ) );


 *l1 = V1/V;
 *l2 = V2/V;
 *l3 = V3/V;
 *l4 = 1.0 - *l3 - *l2 - *l1;

 // retorna verdadeiro quando a particula cai dentro do elemento
 return ( (*l1>=0.0-EPSlocal) && (*l1<=1.0+EPSlocal) && 
          (*l2>=0.0-EPSlocal) && (*l2<=1.0+EPSlocal) && 
		  (*l3>=0.0-EPSlocal) && (*l3<=1.0+EPSlocal) && 
		  (*l4>=0.0-EPSlocal) && (*l4<=1.0+EPSlocal) ); 
}

void SemiLagrangean::setBC()
{
 int aux;
 for( int i=0;i<idbcu->Dim();i++ )
 {
  aux = (int) idbcu->Get(i);
  uParticle.Set( aux,uc->Get(aux) ); 
 }

 for( int i=0;i<idbcv->Dim();i++ )
 {
  aux = (int) idbcv->Get(i);
  vParticle.Set( aux,vc->Get(aux) ); 
 }

 for( int i=0;i<idbcw->Dim();i++ )
 {
  aux = (int) idbcw->Get(i);
  wParticle.Set( aux,wc->Get(aux) ); 
 }

 /* APENAS UTILIZADO NAS SIMULACOES DISK (NUCTE, NUC, NUZ) */
 /* ------------- caso com c.c. livre em w --------------- */
//--------------------------------------------------
//  for( int i=0;i<idbcv->Dim();i++ )
//  {
//   aux = (int) idbcv->Get(i);
//   if( Z->Get(aux) > 1.0 )
//    wParticle.Set( aux,wSol.Get(aux) ); 
//  }
//-------------------------------------------------- 
 /* ------------------------------------------------------ */

 for( int i=0;i<idbcc->Dim();i++ )
 {
  aux = (int) idbcc->Get(i);
  cParticle.Set( aux,cc->Get(aux) ); 
 }
}

void SemiLagrangean::meshInterp(clVector &_X,clVector &_Y,clVector &_Z)
{
 real xP,yP,zP;
 list<int> plist;
 list<int>::iterator mele;
 interpLin.Dim(_X.Dim(),numVerts);

 for (int ii = 0; ii < _X.Dim(); ii++)
 {
  if( ii<numVerts )
   plist = neighbourElem->at(ii);
  else
   plist = neighbourElem->at(0);

  mele=plist.begin(); // pega o primeiro elemento da lista, seja la qual for

  xP = _X.Get(ii);
  yP = _Y.Get(ii);
  zP = _Z.Get(ii);

  //cout << "vertice de origem = " << ii << endl;
  jumpToElem2(*mele,ii,xP,yP,zP);
 } 
} // fim do metodo compute -> getDepartElem

void SemiLagrangean::jumpToElem2(int destElem,int iiVert,real R2X,
                                 real R2Y,real R2Z)
{
 real l1,l2,l3,l4;
 real Bl1,Bl2,Bl3;
 int v[4],v1,v2,v3,v4,vjump;
 int ib1=0;
 int ib2=0;
 int ib3=0;
 if( testElement(destElem,iiVert,R2X,R2Y,R2Z,&l1,&l2,&l3,&l4) )
 {
  v[0] = v1 = (int) IEN->Get(destElem,0);
  v[1] = v2 = (int) IEN->Get(destElem,1);
  v[2] = v3 = (int) IEN->Get(destElem,2);
  v[3] = v4 = (int) IEN->Get(destElem,3);
  interpLin.Set(iiVert,v[0],l1);
  interpLin.Set(iiVert,v[1],l2);
  interpLin.Set(iiVert,v[2],l3);
  interpLin.Set(iiVert,v[3],l4);

  //cout << " Ponto localizado " << iiVert << endl;
 } 
 else 
 {
  v[0] = v1 = (int) IEN->Get(destElem,0);
  v[1] = v2 = (int) IEN->Get(destElem,1);
  v[2] = v3 = (int) IEN->Get(destElem,2);
  v[3] = v4 = (int) IEN->Get(destElem,3);

  if((l1<=l2) && (l1<=l3) && (l1<=l4)){ vjump=0; ib1=v2; ib2=v3; ib3=v4; };
  if((l2<=l1) && (l2<=l3) && (l2<=l4)){ vjump=1; ib1=v1; ib2=v3; ib3=v4; };
  if((l3<=l1) && (l3<=l2) && (l3<=l4)){ vjump=2; ib1=v1; ib2=v2; ib3=v4; };
  if((l4<=l1) && (l4<=l2) && (l4<=l3)){ vjump=3; ib1=v1; ib2=v2; ib3=v3; };

  if ( oFace->Get(destElem,vjump)!=-1)
  {
   jumpToElem2( (int) oFace->Get(destElem,vjump), iiVert,R2X,R2Y,R2Z );
  }
  else
  {
   computeIntercept( iiVert,R2X,R2Y,R2Z,ib1,ib2,ib3,&Bl1,&Bl2,&Bl3 );
   interpLin.Set( iiVert,ib1,Bl1 );
   interpLin.Set( iiVert,ib2,Bl2 );  
   interpLin.Set( iiVert,ib3,Bl3 );  
  }
 }
 return;
}

clMatrix* SemiLagrangean::getInterpLin(){ return &interpLin; }
clVector* SemiLagrangean::getUSL(){ return &uParticle; }
clVector* SemiLagrangean::getVSL(){ return &vParticle; }
clVector* SemiLagrangean::getWSL(){ return &wParticle; }
clVector* SemiLagrangean::getCSL(){ return &cParticle; }


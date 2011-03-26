/**=========================================================================
 *\file Id:TElement3D.h, GESAR
 *\author Gustavo Rabello dos Anjos
 *\date   20-jun-2007
 *\comment Classe TElement 3D para o modulo HYDRO
 *\references
 *\version 1.0
 *\updates
 	version    date         author             comment
	1.0        08/07/2007   Gustavo            finalizacao da classe
 =========================================================================**/

#ifndef TElement_H
#define TElement_H

/* defining element size and rule */
#define NUMGLEU 5  // number of velocity nodes per element
#define NUMGLEP 4  // number of pressure nodes per element
#define NUMRULE 45 // number of quadrature rule
#define NUMGLEC 4  // number of scalar nodes per element
#define NUMRULEC 5 // number of quadrature rule

#include "clVector.h"

class TElement
{
 public:
  
  TElement(); 
  virtual ~TElement();

  double massele[NUMGLEU][NUMGLEU];  ///< matriz de massa do elemento           
  double kxx[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 11
  double kxy[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 12
  double kxz[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 13
  double kyx[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 21
  double kyy[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 22
  double kyz[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 23
  double kzx[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 31
  double kzy[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 32
  double kzz[NUMGLEU][NUMGLEU];      ///< matriz do laplaciano do bloco 33
  double gxele[NUMGLEU][NUMGLEP];    ///< matriz do oper. grad do bloco 1
  double gyele[NUMGLEU][NUMGLEP];    ///< matriz do oper. grad do bloco 2
  double gzele[NUMGLEU][NUMGLEP];    ///< matriz do oper. grad do bloco 3
  double dxele[NUMGLEP][NUMGLEU];    ///< matriz do oper. div do bloco 1
  double dyele[NUMGLEP][NUMGLEU];    ///< matriz do oper. div do bloco 2
  double dzele[NUMGLEP][NUMGLEU];    ///< matriz do oper. div do bloco 3
  double masselec[NUMGLEC][NUMGLEC]; ///< matriz de massa do elemento           
  double kelec[NUMGLEC][NUMGLEC];    ///< matriz do laplaciano do bloco 11
  double kxxc[NUMGLEC][NUMGLEC];     ///< matriz do laplaciano do bloco 11
  double kyyc[NUMGLEC][NUMGLEC];     ///< matriz do laplaciano do bloco 11
  double kzzc[NUMGLEC][NUMGLEC];     ///< matriz do laplaciano do bloco 11
  double gxelec[NUMGLEC][NUMGLEC];   ///< matriz do oper. grad do bloco 1
  double gyelec[NUMGLEC][NUMGLEC];   ///< matriz do oper. grad do bloco 2
  double gzelec[NUMGLEC][NUMGLEC];   ///< matriz do oper. grad do bloco 3
  double dxelec[NUMGLEC][NUMGLEC];   ///< matriz do oper. div do bloco 1
  double dyelec[NUMGLEC][NUMGLEC];   ///< matriz do oper. div do bloco 2
  double dzelec[NUMGLEC][NUMGLEC];   ///< matriz do oper. div do bloco 3

};

#endif


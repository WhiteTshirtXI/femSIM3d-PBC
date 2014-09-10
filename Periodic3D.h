/**
 * \file:   Periodic3D.h
 * \author: Gustavo Charles P. de Oliveira (tavolesliv@gmail.com)
 * \author: Norberto Mangiavacchi (norberto@uerj.br)
 *
 * \brief   Class liable to set PBC (Periodic Boundary Conditions)
 *          over boundaries of domains of simulation.
 *
 * \version 1.0
 *
 * --------------------- CLASS HISTORY ------------------------------
 * \date Created on April 9th, 2013, 17:43 PM
 *
 * ------------------------------------------------------------------
 * \warning DISCLAIMER: This code is under gradual development 
 * 			in collaboration with the following universities/labs: 
 * 			UERJ/GESAR-Brazil, UFRJ/PEMM-Brazil, EPFL/LTCM-Switzerland. 
 * 			As a component of femSIM3D code, is intended to academic 
 * 			purposes only.
 * ------------------------------------------------------------------
 */

#ifndef PERIODIC3D_H
#define	PERIODIC3D_H

// C++ general headers 
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>

// Code headers 
#include "clVector.h"
#include "clMatrix.h"
#include "colors.h"
#include "Model3D.h"
#include "LibTypes.h"
#include "MathTools.h"

class Periodic3D
{
    
public:
    
    Periodic3D(); // default constructor
    Periodic3D(Model3D &_M3D); // overloaded constructor
    virtual ~Periodic3D(); // destructor
    
	/* INT functions */
    int GetNyPointsL();
    int GetNyPointsR();
    
	/* VOID functions */
    void MountPeriodicVectors(string _printIPs);
    void MountPeriodicVectorsNew(string _printIPs);
    void SetVelocityPBC(clVector &_uVelocity, clVector &_vVelocity,
                        clVector &_wVelocity, clVector &_VecXMin,
                        clVector &_VecXMax, int L, string direction);
    void SetVelocityPBCNew(clVector &_uVelocity, clVector &_vVelocity,
                        clVector &_wVelocity, vector<int>* _master,
                        vector<int>* _slave, int _L, string _direction);
    void SetPurePressurePBC(clVector &_Pressure, clVector &_VecXMin,
                            clVector &_VecXMax, int L, string direction);
    void SetPurePressurePBCNew(clVector &_Pressure, vector<int>* _master,
                            vector<int>* _slave, int _L, string _direction);
    void SetJumpPressurePBC(clVector &_Pressure, clVector &_VecXMin,
                            clVector &_VecXMax, int L, double jump);
    void SetPureScalarPBC(clVector &_Scalar, clVector &_VecXMin,
                            clVector &_VecXMax, int L, string direction);
    void SetPureScalarPBCNew(clVector &_Scalar, vector<int>* _master,
                            vector<int>* _slave, int _L, string _direction);
	void SetIndicesVector(vector<int>* _master, vector<int>* _slave);

	/* POINTER functions */
	clVector* GetVecXMin();
    clVector* GetVecXMax();

	vector<int>* GetMasterIndices();
	vector<int>* GetSlaveIndices();
	vector<int>* GetInterfaceIndices();

	list<int>* GetIL();
	list<int>* GetIR();
    
private:
    
	/* INT objects */
    int NumVerts; // number of mesh vertices
    int NumNodes; // number of mesh vertices + centroids
    int nyPointsL; // number of y points over \Gamma_{left}
    int nyPointsR; // number of y points over \Gamma_{right}
    
	/* POINTER objects */
    Model3D *M3DPtr; // pointer to class Model3D
    clVector *XPtr, *YPtr, *ZPtr; // pointers to X,Y,Z vectors from class Model3D
	vector<int> *MasterIndicesPtr;
	vector<int> *SlaveIndicesPtr;
    
	/* CLVECTOR objects */
	clVector VecXMin, VecXMax; // stores indices of points over \Gammas_{left,right}
	
	/* VECTOR OBJECTS */
	vector<int> MasterIndices;
	vector<int> SlaveIndices;
	vector<int> InterfaceIndices;

	list<int> IL;
	list<int> IR;
	
};

#endif	/* PERIODIC3D_H */

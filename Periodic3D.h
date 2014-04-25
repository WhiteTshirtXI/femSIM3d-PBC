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
    int GetNyPointsM();
    int GetNumVertsMid();
    
	/* BOOL functions */
    bool VerifyParallelismY(clVector _u, clVector _v);
    bool VerifyParallelismZ(clVector _u, clVector _v);
    
	/* VOID functions */
    void MountPeriodicVectors(Model3D &_M3D);
    void MountPeriodicVectorsNew(Model3D &_M3D);
    void SetVelocityPBC(clVector &_uVelocity, clVector &_vVelocity,
                        clVector &_wVelocity, clVector &_VecXMin,
                        clVector &_VecXMax, int L, string direction);
    void SetVelocityPBCVector(clVector &_uVelocity, clVector &_vVelocity,
                        clVector &_wVelocity, vector<int> &_master,
                        vector<int> &_slave, int L, string direction);
    void SetPurePressurePBC(clVector &_Pressure, clVector &_VecXMin,
                            clVector &_VecXMax, int L, string direction);
    void SetPurePressurePBCVector(clVector &_Pressure, vector<int> &_master,
                            vector<int> &_slave, int L, string direction);
    void SetJumpPressurePBC(clVector &_Pressure, clVector &_VecXMin,
                            clVector &_VecXMax, int L, double jump);
    void SetPureScalarPBC(clVector &_Scalar, clVector &_VecXMin,
                            clVector &_VecXMax, int L, string direction);
    void SetPureScalarPBCVector(clVector &_Scalar, vector<int> &_master,
                            vector<int> &_slave, int L, string direction);
	void ExtractMiddleVerts(const char *MeshFileName);
	void ForcedParallelismY(int j, clVector _u, clVector _v);
	void ForcedParallelismZ(int j, clVector _u, clVector _v);

	void SetIndicesVector(vector<int>* _master, vector<int>* _slave);
	void SetGeometricalShape(string _shape);

	/* POINTER functions */
	clVector* GetVecXMin();
    clVector* GetVecXMax();
    clVector* GetVecXMid();
    clVector* GetVecXMidVerts();
	vector<int> GetMasterIndices();
	vector<int> GetSlaveIndices();

    
private:
    
	/* INT objects */
    int NumVerts; // number of mesh vertices
    int NumVertsMid; // number of interior mesh vertices without centroids
    int NumNodes; // number of mesh vertices + centroids
    int nyPointsL; // number of y points over \Gamma_{left}
    int nyPointsR; // number of y points over \Gamma_{right}
    int nyPointsM; // number of points over \Omega - ( \Gammas_{left,right} )
    
	string GEOMETRICAL_SHAPE;
    
	/* POINTER objects */
    Model3D *M3DPtr; // pointer to class Model3D
    clVector *XPtr, *YPtr, *ZPtr; // pointers to X,Y,Z vectors from class Model3D
	vector<int> *MasterIndicesPtr;
	vector<int> *SlaveIndicesPtr;
    
	/* CLVECTOR objects */
	clVector VecXMin, VecXMax; // stores indices of points over \Gammas_{left,right}
    clVector VecXMid; // stores indices of points over \Omega - ( \Gammas_{left,right} )
    clVector VecXMidVerts; // store indices according to NumVertsMid
    clVector YLeftBoundaryVector; // stores the y-components of \Gamma_{left}
    clVector YRightBoundaryVector; // stores the y-components of \Gamma_{right}
    clVector ZLeftBoundaryVector; // stores the z-components of \Gamma_{left}
    clVector ZRightBoundaryVector; // stores the z-components of \Gamma_{right}
	
	/* VECTOR OBJECTS */
	vector<int> MasterIndices;
	vector<int> SlaveIndices;
	
};

#endif	/* PERIODIC3D_H */

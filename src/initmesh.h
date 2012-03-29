#ifndef ___INITMESH_H__
#define ___INITMESH_H__
#include <string>

// this data structure can be used both in 2d and 3d situation;

// Define Point type
struct point3d // 3D point 
{
	double x[3];
};
typedef point3d Point;

typedef struct __rbc_type
{
	int* mapping;
	int  rnodenum;
	double* u;
	double* q;
} rbc;

typedef struct __bc_type
{
	int freedom;
	int* bctags;
	double* bcvalues;
} bc;

typedef struct __solu_type
{
	double *u;
	double *q;
} solu;

typedef struct __mesh_type
{
	int nodenum;           // number of boundary nodes 
	int fnodenum;          // number of field nodes
	double** node;         // boundary node (G)
	int* nodeinfo;         // check if the node is the boundary node

	int** elem;            // boundary element
	int** relem;           // radiative boundary element
	int** felem;           // field element

	double** elemnorm;     // boundary element‘ outter norm
	int* eleminfo;         // boundary elements' tags
	int* matinfo;          // tags for mat number 

	int elemnum;           // number of boundary elements
	int relemnum;					 // radiative boundary element
	int felemnum;          // number of field element 
	bc boundaryCondition;
} Mesh;

int initMesh(Mesh& bmesh,std::string filename);

#endif

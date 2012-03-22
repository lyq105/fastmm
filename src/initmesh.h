#ifndef ___INITMESH_H__
#define ___INITMESH_H__
#include <string>

// this data structure can be used both in 2d and 3d situation;
typedef struct __rbc_type
{
	int* mapping;
	int  rnodenum;
	double* u;
	double* q;
} rbc;

struct mesh
{
	int nodenum;           // number of boundary nodes 
	int fnodenum;          // number of field nodes
	double** node;         // boundary node (G)
	int* nodeinfo;         // check if the node is the boundary node

	int** elem;            // boundary element
	int** relem;           // radiative boundary element
	int** felem;           // field element

	double** elemnorm;        // boundary element‘ outter norm
	int* eleminfo;         // boundary elements' tags
	int* matinfo;          // tags for mat number 

	int elemnum;           // number of boundary elements
	int relemnum;					 // radiative boundary element
	int felemnum;          // number of field element 
};

int initMesh(mesh& bmesh,std::string filename);



#endif

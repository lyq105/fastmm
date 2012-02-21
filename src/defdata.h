#ifndef DEF_DATA_H
#define DEF_DATA_H

//======================================================
//             Define some data structures.
//======================================================

enum ElementType 
{
	link1d = 1,
	link2d = 2,
	tri2d  = 3,
	tri3d  = 4,
	quad2d = 5,
	quad3d = 6,
	tet3d  = 7,
	cube3d = 8,
};
enum Efileformat
{
	longformat = 8,
	shortformat = 6,
};

enum boolean
{
	Ture = 1,
	False = 0,
};
// Define Point type
struct point3d // 3D point 
{
	double x[3];
};
typedef point3d Point;

struct freedomstatus
{
	boolean valueisgiven;
	boolean fluxisgiven;
};

typedef freedomstatus FreedomStatus;


// 
// Because of both value and flux cannot be given at same time. So
// boundary condition is defined as below.
//

struct boundarycondition
{
	int totlefreedom;
	int* freedomnumber;
	FreedomStatus* freedomstatlist;
	double* values;
};

typedef boundarycondition BoundaryCondition;

struct cell
{
	int* cellnumber;
};
typedef cell Cell;

struct mesh
{
	int dim;
	ElementType et;
	int nodepercell;
	int TotleNodeNumber;
	int TotleElementNumber;
	Point* MeshPointList;
	Cell* CellList; 
	BoundaryCondition bc;
};

typedef mesh Mesh;



//==============================================================
//                  Define some global data;
//==============================================================

// Definitions of Utility Value 
//      File formats of data file use in reading mesh data; 

const double triref[3][2] = {{0,0},{1,0},{0,1}};
const double tetref[4][3] = {{0,0,0},{0,1,0},{0,0,1},{0,0,1}};

const double PI = 3.14159265;
#endif //DEF_DATA_H

#ifndef __LAPACE2D_0_H_
#define __LAPACE2D_0_H_
#include "initmesh.h"

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

void GAUSSELIMINATE(double **E,int M,double *RHS);//适用于全带宽存储阵
int laplace2d_0(mesh bmesh, bc bc1, solu& solution);
int print_temp(mesh bmesh, solu solution,string resfname);
double u_boundary(int number, mesh bmesh, solu solution);
double u_field(double sx0, double sx1, mesh bmesh, solu solution);

#endif

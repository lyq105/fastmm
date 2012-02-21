//
// A:  liyiqiang(lyq105@163.com)
// P:  test if mesh function is ok.
// C:  2011-07-20 10:23:47
// M:  2011-07-20 10:23:47
// 
#include <stdio.h>
#include <stdlib.h>
#include "defdata.h"
#include "mesh.h"
#include "quadtree.h"
extern int maxTreeDepth;
extern int maxEleminCell;

int main(int argc, const char *argv[])
{
/*
	if(argc<=1)
	{
		printf("Please input you finite element mesh file names!!!\n");
		return 0;
	}
	*/
	char nodefilename []= "annulus.n";
	char elemfilename []= "annulus.e";
	char outfilename []= "link2d.plt";

	Mesh mesh;
	Quadtree qtree;

	if(!InitialMesh(nodefilename,elemfilename,&mesh,link2d))
	{
		printf("Failed to reading files!\n");
		return 0;
	};
  //maxTreeDepth = 20;
  //maxEleminCell = 5;
	quadtree_creat(qtree,mesh);	

	Mesh2tecplot(&mesh,link2d,outfilename);
	FreeMeshMem(&mesh);
	return 1;
}
	
/**************************************************************
 * check boundary element's norm.
 * check if the norm of boundary element is point to the 
 * interior of domain.
 **************************************************************/
int checkBoundaryNorm(Mesh& mesh)
{
	int i,j,k;
	int ten = mesh.TotleElementNumber;
	for (i = 0; i < ten; i++) 
	{
	}
	return ten;
}

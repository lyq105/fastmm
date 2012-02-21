/*------------------------------------------------------------
A:  liyiqiang(lyq105@163.com)
P:  Read data from ansys mesh generator.
C:  2011-03-10 15:33:20
M:  2011-03-10 15:33:20
------------------------------------------------------------*/
//#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include "defdata.h"
#include <string.h>

//#define _DEBUG_
//#define _USE_CPP_



int InitialMesh(char* nodefilename,char* elemfilename,Mesh* mesh,ElementType et,Efileformat efmt = shortformat)
{
	FILE* fp;
	FILE* efp;

	char buffer[200];
	char x[3][21]={'\0'};
	char elem[8][10]={'\0'};
	char ch;

	switch(et)
	{
		case link1d:
			mesh->dim=1;
			mesh->nodepercell=2;
			break;
		case link2d:
			mesh->dim=2;
			mesh->nodepercell=2;
			break;
		case tri2d:
			mesh->dim=2;
			mesh->nodepercell=3;
			break;
		case tri3d:
			mesh->dim=3;
			mesh->nodepercell=3;
			break;
		case quad2d:
			mesh->dim=2;
			mesh->nodepercell=4;
			break;
		case quad3d:
			mesh->dim=3;
			mesh->nodepercell=4;
			break;
		case tet3d:
			mesh->dim=3;
			mesh->nodepercell=4;
			break;
		case cube3d:
			mesh->dim=3;
			mesh->nodepercell=8;
			break;
	}

	fp=fopen(nodefilename,"r");
	if (fp == NULL)
	{
		printf("Cannot open file %s!!\n",nodefilename);
		return 0;
	}
	efp=fopen(elemfilename,"r");
	if (efp == NULL)
	{
		printf("Cannot open file %s!!\n",elemfilename);
		return 0;
	}


	int nnum=0;
	while((ch=fgetc(fp))!=EOF)
	{
		if(ch == '\n')nnum++;
	}
	rewind(fp);


	int elemnum=0;
	while((ch=fgetc(efp))!=EOF)
	{
		if(ch=='\n')elemnum++;
	}
	rewind(efp);

	mesh->TotleNodeNumber = nnum;
	mesh->TotleElementNumber = elemnum;

#ifdef _DEBUG_
	system("pause");
	printf("There are %d nodes in this mesh!\n",mesh->TotleNodeNumber);
	printf("there are %d element in this mesh!!\n",mesh->TotleElementNumber);
	system("pause");
#endif


#ifdef _USE_CPP_
	mesh->MeshPointList = new Point[mesh->TotleNodeNumber];
	mesh->CellList = new Cell[mesh->TotleElementNumber];
#else  // using Ansi C.
	mesh->MeshPointList = (Point*)malloc(sizeof(Point)*mesh->TotleNodeNumber);
	mesh->CellList = (Cell*)malloc(sizeof(Cell)*mesh->TotleElementNumber);
#endif


	int pindex = 0;
	while(fgets(buffer,200,fp)!=NULL)
	{
		for (int j=0;j<mesh->dim;j++)
		{
			x[j][0]='\0';
			if (strlen(buffer) == 8) continue;
			for(int i=0;i<20;i++)
			{
				int strid=j*20+i+8;
				if(buffer[strid] != '\0')
				{
					x[j][i] = buffer[strid];
					x[j][i+1]='\0';
				}
				else
				{
					break;
				}
			}
			mesh->MeshPointList[pindex].x[j]=atof(x[j]);
		}

		pindex++;

#ifdef _DEBUG_
		//system("pause");
		printf("%8d ",pindex-1);
		for (int j=0;j<mesh->dim;j++)
			printf("%12.8f ",atof(x[j]));
		printf("\n");
#endif
	}

	fclose(fp);

	int eindex=0;
	while(fgets(buffer,200,efp)!=NULL)
	{

		// Allocate memory of every cell.

#ifdef _USE_CPP_
		mesh->CellList[eindex].cellnumber = new int[mesh->nodepercell];
#else  // using Ansi C.
		mesh->CellList[eindex].cellnumber = (int*) malloc(sizeof(int)*mesh->nodepercell);
#endif
		for (int j=0;j < mesh->nodepercell;j++)
		{
			elem[j][0]='\0';
			for(int i=0;i<efmt;i++)
			{
				int strid=j*efmt+i;
				if(et == tet3d && j == 3)
				{
					strid = 4*efmt+i;
				}
				if(buffer[strid] != '\0')
				{
					elem[j][i] = buffer[strid];
					elem[j][i+1]='\0';
				}
			}
			mesh->CellList[eindex].cellnumber[j]=atoi(elem[j])-1;
		}

		eindex++; // element index ++;

#ifdef _DEBUG_
		for (int j=0;j< nodenum;j++)
			printf("%s ",(elem[j]));
		printf("\n");
#endif
	}

	fclose(efp);
	return 1;
}
int InitBoundaryCondition(BoundaryCondition* bc,char* bcfilename)
{
	char buffer[200];
	char x[2][20];
	char ch;
	FILE* fp;
	fp=fopen(bcfilename,"r");
	if(!fp)
	{
		printf("Cann't open file %s\n",bcfilename);
		return 0;
	}

	int nnum=0;
	while((ch=fgetc(fp))!=EOF)
	{
		if(ch == '\n')nnum++;
	}
	rewind(fp);

	bc->totlefreedom = nnum;


#ifdef _USE_CPP_
	bc->freedomnumber = new int[bc->totlefreedom];
	bc->freedomstatlist = new FreedomStatus[bc->totlefreedom];
	bc->values = new double[bc->totlefreedom];
#else  // using Ansi C.
	bc->freedomnumber = (int*)malloc(sizeof(int)*bc->totlefreedom);
	bc->freedomstatlist =(FreedomStatus*)malloc(sizeof(FreedomStatus)*bc->totlefreedom);
	bc->values = (double*)malloc(sizeof(double)*bc->totlefreedom);
#endif
	int bcindex = 0;
	while(fgets(buffer,200,fp)!=NULL)
	{
		for (int j=0;j<2;j++)
		{
			x[j][0]='\0';
			for(int i=0;i<8;i++)
			{
				int strid=j*8+i;
				if(buffer[strid] != '\0')
				{
					x[j][i] = buffer[strid];
					x[j][i+1]='\0';
				}
			}
		}
		bc->freedomnumber[bcindex]=atoi(x[0]);
		if (atoi(x[1]) == 1)
		{
			bc->freedomstatlist[bcindex].valueisgiven = False;
			bc->freedomstatlist[bcindex].fluxisgiven = Ture;
		}
		else if(atoi(x[1]) == 2)
		{
			bc->freedomstatlist[bcindex].valueisgiven = Ture;
			bc->freedomstatlist[bcindex].fluxisgiven = False;
		}
		bcindex++;
	}
	return 0;
}

int FreeMeshMem(Mesh* mesh)
{
#ifdef _USE_CPP_
	delete [] mesh->MeshPointList;

	for (int eindex=0;eindex<mesh->TotleElementNumber;eindex++)
	{
		delete [] mesh->CellList[eindex].cellnumber;
	}

	delete [] mesh->CellList;
#else  // using Ansi C.
	free(mesh->MeshPointList);
	mesh->MeshPointList = NULL;
	for (int eindex=0;eindex<mesh->TotleElementNumber;eindex++)
	{
		free(mesh->CellList[eindex].cellnumber);
		mesh->CellList[eindex].cellnumber = NULL;
	}
	free(mesh->CellList);
	mesh->CellList = NULL;
#endif
	return 0;
}

int FreeBCMem(BoundaryCondition* bc)
{
#ifdef _USE_CPP_
	delete [] bc->freedomnumber;
	bc->freedomnumber=NULL;
	delete [] bc->freedomstatlist;
	bc->freedomstatlist=NULL;
	delete [] bc->values;
	bc->values=NULL;

#else  // using Ansi C.
	free(bc->freedomnumber);
	bc->freedomnumber=NULL;
	free(bc->freedomstatlist);
	bc->freedomstatlist=NULL;
	free(bc->values);
	bc->values=NULL;
#endif
	return 0;
}
int Mesh2tecplot(Mesh* mesh,ElementType et,char* outfilename)
{

	const char quadheader[] = "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n";
	const char triheader[] = "ZONE N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n";
	const char tetheader[] = "ZONE N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n";
	const char cubeheader[] = "ZONE N=%d, E=%d, F=FEPOINT, ET=BRICK\n";
	FILE* tecplot;
	tecplot=fopen(outfilename,"w");
	if (tecplot == NULL)
	{
		printf("Cannot open file %s!!\n",outfilename);
		return 0;
	}

	char* header;
	switch(et)
	{
		case link1d:
		case link2d:
		case quad2d:
		case quad3d:
			fprintf(tecplot,quadheader,mesh->TotleNodeNumber,mesh->TotleElementNumber);
			break;
		case tri2d:
		case tri3d:
			fprintf(tecplot,triheader,mesh->TotleNodeNumber,mesh->TotleElementNumber);
			break;
		case tet3d:
			fprintf(tecplot,tetheader,mesh->TotleNodeNumber,mesh->TotleElementNumber);
			break;
		case cube3d:
			fprintf(tecplot,cubeheader,mesh->TotleNodeNumber,mesh->TotleElementNumber);
			break;
	}
	for(int nindex=0;nindex < mesh->TotleNodeNumber;nindex++)
	{
		for (int ndindex=0;ndindex<mesh->dim;ndindex++)
		{
			fprintf(tecplot,"%12.8lf ",mesh->MeshPointList[nindex].x[ndindex]);
		}
		fprintf(tecplot,"\n");
	}
	for(int eindex=0;eindex< mesh->TotleElementNumber;eindex++)
	{
		for (int ndindex=0;ndindex< mesh->nodepercell;ndindex++)
		{
			fprintf(tecplot,"%9d ",mesh->CellList[eindex].cellnumber[ndindex]+1);
			if (et == link2d)
				fprintf(tecplot,"%9d ",mesh->CellList[eindex].cellnumber[ndindex]+1);
		}
		fprintf(tecplot,"\n");
	}
	return 1;
}



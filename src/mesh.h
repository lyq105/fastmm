#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <stdlib.h>
#include "defdata.h"

int InitialMesh(char* nodefilename,char* elemfilename,Mesh* mesh,ElementType et,Efileformat efmt = shortformat);
int Mesh2tecplot(Mesh* mesh,ElementType et,char* outfilename);
int FreeMeshMem(Mesh* mesh);

#endif // _MESH_H_

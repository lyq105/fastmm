#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "initmesh.h"
using namespace std;
// NOTE: this code may cause memory leaks.

template <typename T>
int allocatearray( T** &a,int fdim,int sdim)
{
	int i,j;
	a = new T*[fdim];
	for (i = 0; i < fdim; i++) {
		a[i] = new T[sdim];
		for (j = 0; j < sdim; j++) {
			a[i][j] = 0;
		}
	}
	return 0;
}
template <typename T>
int allocatevector(T* &v,int dim)
{
	int i;
	v = new T[dim];
	for (i = 0; i < dim; i++) {
		v[i] = 0;
	}
}

int initMesh(Mesh& bmesh,string filename)
{
	string buffer;
	istringstream sPeaser;
	ifstream ifile(filename.c_str(),ios::in);

	getline(ifile,buffer);
	sPeaser.str(buffer);
	sPeaser >> bmesh.nodenum 	>> bmesh.fnodenum >> bmesh.elemnum	>> bmesh.relemnum	>> bmesh.felemnum;
  //cout << bmesh.elemnum ;

  allocatearray(bmesh.node,bmesh.nodenum,2);
	allocatearray(bmesh.relem,bmesh.relemnum,2);
	allocatearray(bmesh.elem,bmesh.elemnum,2);
	allocatearray(bmesh.elemnorm,bmesh.elemnum,2);
	allocatearray(bmesh.felem,bmesh.felemnum,3);
	allocatevector(bmesh.eleminfo,bmesh.elemnum + bmesh.relemnum);
	allocatevector(bmesh.matinfo,bmesh.elemnum + bmesh.relemnum);
	allocatevector(bmesh.nodeinfo,bmesh.nodenum);


	for (int i = 0; i < bmesh.nodenum ; i++) 
	{
		int number,info;
		double x,y;
		getline(ifile,buffer);
		istringstream speaser(buffer);
		speaser >> number >> x >> y >> info;
    //cout << number <<" "<< x <<" "<< y << endl;
		bmesh.node[number-1][0] = x;
		bmesh.node[number-1][1] = y;
		bmesh.nodeinfo[number-1] = info;
	}

	for (int i = 0; i < bmesh.elemnum; i++) 
	{
		int number,e1,e2,mat,real;
		double normx,normy;
		getline(ifile,buffer);
		istringstream speaser(buffer);
		speaser >> real >> e1 >> e2 >> normx >>normy>> mat;
    //cout << real << " "<< e1 << " "<< e2 << " "<< normx << " "<< normy << " "<< mat << " "<< endl;
		bmesh.elem[i][0] = e1-1;
		bmesh.elem[i][1] = e2-1;
   
    bmesh.elemnorm[i][0] = normx;
		bmesh.elemnorm[i][1] = normy;
    //cout << bmesh.elemnorm[i][0] << " == " << bmesh.elemnorm[i][1]<<endl;
		bmesh.eleminfo[i] = real;
		bmesh.matinfo[i] = mat;
	}
	

	for (int i = 0; i < bmesh.felemnum; i++) 
	{
		int number,e1,e2,e3;
		getline(ifile,buffer);
		istringstream speaser(buffer);
		speaser >> number >> e1 >> e2 >> e3;
		//cout << number << " "<< e1<< " " << e2 << " "<< e3 << endl;
		bmesh.felem[i][0] = e1 - 1;
		bmesh.felem[i][1] = e2 - 1;
		bmesh.felem[i][2] = e3 - 1;
	}

	ofstream tec("tec.plt");
  
	tec << "zone " << "n=" << bmesh.nodenum 
								 << ",e= " << bmesh.elemnum
								 << ",et= triangle"
								 << ",f = fepoint" << endl;

	for (int i = 0; i< bmesh.nodenum; i++)
	{
		tec << bmesh.node[i][0] << " " << bmesh.node[i][1] << endl;
	}
	for (int i = 0; i< bmesh.elemnum; i++)
	{
		tec << bmesh.elem[i][0] +1 << " " << bmesh.elem[i][1]+1 <<" " << bmesh.elem[i][1]+1<< endl;
	}

	return 0;
}


//int main(int argc, const char *argv[])
//{
  //string datafile("bem2d.dat");
  //mesh me;
  //initMesh(me,datafile);
  //return 0;
//}



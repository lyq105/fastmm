#include <iostream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "initmesh.h"
#include "laplace2d_0.h"

const double pi = 3.14159265;
const double sigma = 5.670373E-8 ;
const double epsilong = 0.05;


// describe the 1st and 2nd boundary condition
int initialBC(bc& bc1,Mesh bmesh)
{
	bc1.freedom = bmesh.elemnum;

	bc1.bctags = new int [bc1.freedom];
	bc1.bcvalues = new double [bc1.freedom];

	for (int i = 0; i < bc1.freedom; i++)
	{	
   /* 
	for (int i = 0; i < bc1.freedom; i++)
	{		
		if(bmesh.eleminfo[i] == 1)    // first kind boundary condition
		{
			bc1.bctags[i] = 1;
			bc1.bcvalues[i] = 100;
		}
		if(bmesh.eleminfo[i] == 2) // second kind boundary condition
		{
			bc1.bctags[i] = 2;
			bc1.bcvalues[i] = 200;
		}
	}
  */
  
     int	nx0 = bmesh.elem[i][0];
		int nx1 = bmesh.elem[i][1];

    //cout << "imap " << imap  << " "<< nx0 <<" " <<nx1<< endl;
		double source_x0 = 0.5*(bmesh.node[nx0][0] + bmesh.node[nx1][0]);
		double source_x1 = 0.5*(bmesh.node[nx0][1] + bmesh.node[nx1][1]); 
		bc1.bctags[i] = 2;
		bc1.bcvalues[i] = 0;
    //if (fabs(source_x1+1)<1e-19 || fabs(source_x0+1)<1e-19)
    if (fabs(source_x1+1)<1e-19)
    {
      cout << i <<" "<< nx0 <<" " <<nx1<< " apply 1" <<endl;
			bc1.bctags[i] = 1;
			bc1.bcvalues[i] = 300;
    }

    //if (fabs(source_x1-1)<1e-19 || fabs(source_x0 -1)<1e-19)
    if (fabs(source_x1-1)<1e-19 )
    {
      cout << i <<" "<< nx0 <<" " <<nx1<< " apply 2" <<endl;
			bc1.bctags[i] = 1;
			bc1.bcvalues[i] = 100;
		}
    
  }
  return 0;
}

int assignRBCvalue(bc& bc1, rbc rbc1 )
{
  for (int i = 0; i < rbc1.rnodenum; i++) 
  {
			bc1.bctags[rbc1.mapping[i]] = 2;
			bc1.bcvalues[rbc1.mapping[i]] = rbc1.q[i];
  }
  return 0;
}

int getRBCvalue(solu solution, rbc& rbc1 )
{
  for (int i = 0; i < rbc1.rnodenum; i++) 
  {
		rbc1.u[i] = solution.u[rbc1.mapping[i]];
  }
  return 0;
}

int initRbc(Mesh bmesh, rbc* rbc1, int rbcnum)
{
	for (int rindex = 0; rindex < rbcnum; rindex++) 
	{
		int k = 0;
		for (int i = 0; i<bmesh.elemnum; i++)
		{
			//cout << "mat  " <<bmesh.matinfo[i] <<endl;
			//cout << "real  " <<bmesh.eleminfo[i] <<endl;
			if (bmesh.matinfo[i] == rindex + 1 && bmesh.eleminfo[i] == 1)
			{
				k++;
			}
		}
		rbc1[rindex].rnodenum = k;
		rbc1[rindex].mapping = new int[rbc1[rindex].rnodenum];
		rbc1[rindex].u = new double[rbc1[rindex].rnodenum];
		rbc1[rindex].q = new double[rbc1[rindex].rnodenum];
		k = 0;
		for (int i = 0; i<bmesh.elemnum; i++)
		{
			if (bmesh.matinfo[i] == rindex + 1 && bmesh.eleminfo[i] == 1)
			{
				rbc1[rindex].mapping[k] = i;
				//cout << i <<end;
				k++;
			}
		}
	}
	return 0;
}

int qr_boundary(Mesh bmesh,rbc& rbc1)
{
	int i,j;
	double hij,gij;
	double x,y;
	double unx,uny;

	double vf=0;

	double gaussnum = 7;
	double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
	double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};
	double** matrix;
	double* b;

	matrix = new double* [rbc1.rnodenum];
	for (int i = 0; i < rbc1.rnodenum; i++) 
	{
		matrix[i]= new double[rbc1.rnodenum];
		for (int j = 0; j < rbc1.rnodenum; j++) {
			matrix[i][j] = 0;
		}
	}


	b = new double[rbc1.rnodenum];

	double unx1,uny1;
	double unx2,uny2;
	for(i = 0; i< rbc1.rnodenum; i++) 
	{
		// find out the source point
		int imap = rbc1.mapping[i];
		int	nx0 = bmesh.elem[imap][0];
		int nx1 = bmesh.elem[imap][1];

    //cout << "imap " << imap  << " "<< nx0 <<" " <<nx1<< endl;
		double source_x0 = 0.5*(bmesh.node[nx0][0] + bmesh.node[nx1][0]);
		double source_x1 = 0.5*(bmesh.node[nx0][1] + bmesh.node[nx1][1]);

		unx1 = bmesh.elemnorm[imap][0];
		uny1 = bmesh.elemnorm[imap][1];
//cout << "(" << unx1 << "," << uny1 << ")"<<endl;
		double btemp = 0;
		vf =0;
		for (j = 0; j < rbc1.rnodenum; j++) 
		{
			// find out the number of 2 ends of a cell.
			int jmap = rbc1.mapping[j];
			int	ny0 = bmesh.elem[jmap][0];
			int ny1 = bmesh.elem[jmap][1];
			double x0[2],x1[2];	
			x0[0] = bmesh.node[ny0][0];
			x0[1] = bmesh.node[ny0][1];
			x1[0] = bmesh.node[ny1][0];
			x1[1] = bmesh.node[ny1][1];

			// calculate the length of a cell.
			double gama_j = (x0[0] - x1[0])*(x0[0] - x1[0])+(x0[1] - x1[1])*(x0[1] - x1[1]);
			gama_j = sqrt(gama_j);
			//cout << gama_j <<endl;
			unx2 = bmesh.elemnorm[jmap][0];
			uny2 = bmesh.elemnorm[jmap][1];
      //if ( i == 1)
        //cout << "(" << unx2 << "," << uny2 << ")"<<endl;
			if (i == j)
			{
				matrix[i][j]=1;
				btemp += epsilong*sigma*pow(rbc1.u[j],4);
			}
			else
			{
				double mtemp=0;
				for(int gindex = 0; gindex< gaussnum; gindex++)
				{
					double yy0 = 0.5*(x0[0]-x1[0])*gausspt[gindex] + 0.5*(x0[0] + x1[0]);
					double yy1 = 0.5*(x0[1]-x1[1])*gausspt[gindex] + 0.5*(x0[1] + x1[1]);
					double fz = (unx1*(source_x0 - yy0)+uny1*(source_x1 - yy1))*(unx2*(yy0-source_x0)+uny2*(yy1-source_x1));
          //cout << "fz = " << fz << endl;
          //if (fz < 0 ) cout << fz <<endl;
					double rr = ((yy0-source_x0)*(yy0-source_x0) + (yy1-source_x1)*(yy1-source_x1));
					//cout << "rr" << rr<<endl;
					mtemp += fabs(fz)/(2*rr*sqrt(rr))*gausswt[gindex];
				}
				vf += mtemp*gama_j*0.5;
				matrix[i][j] = -epsilong*(1/epsilong -1)*gama_j*mtemp*0.5;
        btemp -= epsilong*sigma*pow(rbc1.u[j],4)*(gama_j*0.5*mtemp);
				//btemp -= epsilong*sigma*bc1.bcvalues[j]*(gama_j*0.5*mtemp);
			}
		}	
    //cout <<"1 - vf = " <<1 - vf <<endl;
		//cout <<"vf = " << vf <<endl;
		b[i] = btemp;
	} // loop for element.

	/// step 3. solving the linear system.
	GAUSSELIMINATE(matrix,rbc1.rnodenum,b);	
	/// step 4. post processing of this problem, such as output the results.
	for (i = 0; i < rbc1.rnodenum; i++) 
	{
		rbc1.q[i] = -b[i];
	}
	return 0;
}


int main(int argc, const char *argv[])
{
	string filename = "bem2d.dat";
	if (argc >=2)
	{
		filename = argv[1];
	}

	Mesh bmesh;
	bc bc1;
	rbc rbc1;
	solu solution,solu0;

  initMesh(bmesh,filename);
	//handleDof();
	initialBC(bc1,bmesh);
	initRbc(bmesh,&rbc1,1);

	//cout << "return " << rbc1.rnodenum << endl;


	for (int i = 0; i < rbc1.rnodenum/2; i++) {
		//rbc1.u[i] = 10;
	}
	for (int i = rbc1.rnodenum/2; i < rbc1.rnodenum; i++) {
		//rbc1.u[i] = 20;
	}
	//qr_boundary(bmesh,rbc1);

  //assignRBCvalue(bc& bc1, rbc rbc1);
	solution.q = new double [bmesh.elemnum];
	solution.u = new double [bmesh.elemnum];
	solu0.q = new double [bmesh.elemnum];
	solu0.u = new double [bmesh.elemnum];

	laplace2d_0(bmesh,bc1,solu0);

  print_temp(bmesh,solu0,filename + "_nn.plt");
	char buffer[10];
	for(int i=0; i<80; i++)
	{
		sprintf(buffer,"%d",i);
		laplace2d_0(bmesh,bc1,solution);
		print_temp(bmesh,solution,filename + buffer + "_1.plt");
		getRBCvalue(solution,rbc1);
		qr_boundary(bmesh,rbc1);
		for (int i = 0; i < rbc1.rnodenum; i++) 
		{
			//cout << rbc1.q[i] <<endl;
		}
		assignRBCvalue(bc1,rbc1);
	}
	laplace2d_0(bmesh,bc1,solution);
  print_temp(bmesh,solution,filename + ".plt");
	return 0;
}


//=============================================================================
//   Project name:  constant element bem for solving 2D heat conduct problem.
//       Filename:  
//         Author:  liyiqiang(lyq105@163.com)
//         Create:  2011-06-22 21:16:26
//      Last edit:  2011-06-22 21:16:26
//=============================================================================

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
using namespace std;
#include "initmesh.h"
#include "laplace2d_0.h"





double gaussnum = 4;
double gausspt[4] = { 0.86113631, -0.86113631, 0.33998014, -0.33998014 };
double gausswt[4] = {0.34785485,0.34785485,0.65214515,0.65214515};

//double gaussnum = 7;
//double gausspt[7] = {0, 0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515,-0.40584515 };
//double gausswt[7] = {0.41795918,0.12948497,0.12948497,0.27970539,0.27970539,0.38183005,0.38183005};

const double pi = 3.14159265;


int laplace2d_0(mesh bmesh, bc bc1, solu& solution)
//int solveBem2d(int argc, const char *argv[])
{
	int i,j;
	double hij,gij;
	double x,y;
	double unx,uny;
	double *sx0,*sx1;

	double** matrix;
	double*  b;
	double*  u;
	//// substep 1 

	matrix = new double* [bmesh.elemnum];
	for (int i = 0; i < bmesh.elemnum; i++) 
	{
		matrix[i]= new double[bmesh.elemnum];
		for (int j = 0; j < bmesh.elemnum; j++) {
			matrix[i][j] = 0;
		}
	}

	b = new double[bmesh.elemnum];
	u = new double[bmesh.elemnum];
	sx0 = new double[bmesh.elemnum];
	sx1 = new double[bmesh.elemnum];


	for (i = 0; i < bmesh.elemnum; i++) 
	{
		b[i] = 0;
		u[i] = 0;
		sx0[i] = 0; 
		sx1[i] = 0; 
	}

	for (i = 0; i < bmesh.elemnum; i++) 
	{
		int	nindex0 = bmesh.elem[i][0];
		int nindex1 = bmesh.elem[i][1];

		// find out the source point.
		sx0[i] = 0.5*(bmesh.node[nindex0][0] + bmesh.node[nindex1][0]);
		sx1[i] = 0.5*(bmesh.node[nindex0][1] + bmesh.node[nindex1][1]);
	}

	for(i = 0; i< bmesh.elemnum; i++) 
	{
		// find out the source point
		double source_x0 = sx0[i];
		double source_x1 = sx1[i];
		double b_temp = 0;
		for (j = 0; j < bmesh.elemnum; j++) 
		{
			// find out the number of 2 ends of a cell.
			int	nindex0 = bmesh.elem[j][0];
			int nindex1 = bmesh.elem[j][1];
			double x0[2],x1[2];	
			x0[0] = bmesh.node[nindex0][0];
			x0[1] = bmesh.node[nindex0][1];
			x1[0] = bmesh.node[nindex1][0];
			x1[1] = bmesh.node[nindex1][1];

			// calculate the length of a cell.
			double gama_j = (x0[0] - x1[0])*(x0[0] - x1[0])+(x0[1] - x1[1])*(x0[1] - x1[1]);

			gama_j = sqrt(gama_j);
			// calculate the outter unit normal vector of the cell.
      unx = bmesh.elemnorm[j][0];
      uny = bmesh.elemnorm[j][1];
      //unx =  (bmesh.node[nindex1][1] - bmesh.node[nindex0][1])/gama_j;
      //uny =  -(bmesh.node[nindex1][0] - bmesh.node[nindex0][0])/gama_j;

			if(i == j) // diagonal element of matrix.
			{
				hij = 1.0/2;
        gij = gama_j / (2*pi) * (log(2.0/gama_j) + 1);
			}
			else
			{
				// use gauss integral 
				hij = 0;
				gij = 0;
				for(int gindex = 0; gindex< gaussnum; gindex++)
				{
					x = 0.5*(x0[0]-x1[0])*gausspt[gindex] + 0.5*(x0[0] + x1[0]);
					y = 0.5*(x0[1]-x1[1])*gausspt[gindex] + 0.5*(x0[1] + x1[1]);
					double rr = ((x-source_x0)*(x-source_x0) + (y-source_x1)*(y-source_x1));
					hij += 1.0/ rr * gausswt[gindex];
					gij += -log(sqrt(rr)) *gausswt[gindex];
				}
				double d = (x-source_x0)*unx + (y-source_x1)*uny;

				hij *= -1.0/(4*pi)*gama_j*d;
				gij *= 1.0/(4*pi)*gama_j;
			}
			// assemble the matrix and right hand side vector;

			if(bc1.bctags[j] == 1)
			{
				matrix[i][j] += -gij;
				b_temp += -hij*bc1.bcvalues[j];
			}
			else if (bc1.bctags[j] == 2)
			{
				matrix[i][j] += hij;
				b_temp += gij*bc1.bcvalues[j];
			}
		}//loop for source point.
		b[i] = b_temp;
	} // loop for element.

	/// step 3. solving the linear system.
	GAUSSELIMINATE(matrix,bmesh.elemnum,b);	
	/// step 4. post processing of this problem, such as output the results.
 
	for (i = 0; i < bmesh.elemnum; i++) 
	{
		if (bc1.bctags[i] == 1)
		{
			solution.u[i] = bc1.bcvalues[i];
			solution.q[i] = b[i];
		}
		if (bc1.bctags[i] == 2)
		{
			solution.u[i] = b[i];
			solution.q[i] = bc1.bcvalues[i];
		}
	}

  //print_temp(bmesh,bc1,b);
	return 0;
}

int print_temp(mesh bmesh, solu solution,string resfname)
{
	FILE* fp_fnode;
  fp_fnode = fopen(resfname.c_str(),"w");

	fprintf(fp_fnode, "zone f=fepoint,e=%d,n=%d,et=triangle\n",bmesh.felemnum,bmesh.nodenum);

	for (int i = 0; i< bmesh.nodenum; i++)
	{
		double source_x0 = bmesh.node[i][0];
		double source_x1 = bmesh.node[i][1];

		if(bmesh.nodeinfo[i] == 0)   // boundary point
		{
      double u = u_boundary(i,bmesh,solution);
      fprintf(fp_fnode, "%12.8f %12.8f %12.8f\n",source_x0,source_x1,u);
		}
    else
    {
      double u = u_field(source_x0,source_x1,bmesh,solution);
      fprintf(fp_fnode, "%12.8f %12.8f %12.8f\n",source_x0,source_x1,u);
    }
  }
  for (int i = 0; i < bmesh.felemnum; i++) 
  {
    fprintf(fp_fnode, "%8d %8d %8d\n",bmesh.felem[i][0]+1,bmesh.felem[i][1]+1,bmesh.felem[i][2]+1);
  }
	fclose(fp_fnode);
	return 0;
}

double u_boundary(int number, mesh bmesh, solu solution)
{
  double k = 0;
  double sum = 0;
  for(int i = 0; i<bmesh.elemnum;i++)
  {
    if(number == bmesh.elem[i][0] || number == bmesh.elem[i][1])
    {
      k += 1;
			sum += solution.u[i];
    }
  }
  return sum/k;
}
double u_field(double sx0, double sx1, mesh bmesh, solu solution)
{
	int j;
  double unx,uny;
  double hij,gij;
  double rr,x,y,d;
	double ret = 0;
	for (j = 0; j < bmesh.elemnum; j++) 
	{
		// find out the number of 2 ends of a cell.
		int	nindex0 = bmesh.elem[j][0];
		int nindex1 = bmesh.elem[j][1];
		double x0[2],x1[2];	
		x0[0] = bmesh.node[nindex0][0];
		x0[1] = bmesh.node[nindex0][1];
		x1[0] = bmesh.node[nindex1][0];
		x1[1] = bmesh.node[nindex1][1];

		// calculate the length of a cell.
		double gama_j = (x0[0] - x1[0])*(x0[0] - x1[0])+(x0[1] - x1[1])*(x0[1] - x1[1]);

		gama_j = sqrt(gama_j);
    unx = bmesh.elemnorm[j][0];
    uny = bmesh.elemnorm[j][1];
		//unx =  (bmesh.node[nindex1][1] - bmesh.node[nindex0][1])/gama_j;
		//uny =  -(bmesh.node[nindex1][0] - bmesh.node[nindex0][0])/gama_j;
		//cout << "(" << unx << "," << uny << ")\n";
    // use gauss integral 
    hij = 0;
    gij = 0;
    for(int gindex = 0; gindex< gaussnum; gindex++)
    {
       x = 0.5*(x0[0]-x1[0])*gausspt[gindex] + 0.5*(x0[0] + x1[0]);
       y = 0.5*(x0[1]-x1[1])*gausspt[gindex] + 0.5*(x0[1] + x1[1]);
       rr = ((x-sx0)*(x-sx0) + (y-sx1)*(y-sx1));
			 hij += 1.0/ rr * gausswt[gindex];
			 gij += -log(sqrt(rr)) *gausswt[gindex];
    }
    d = (x-sx0)*unx + (y-sx1)*uny;
    //cout <<"d = " <<d <<endl;
    hij *= -1.0/(4*pi)*gama_j*d;
    gij *= 1.0/(4*pi)*gama_j;

		ret += -hij*solution.u[j] + gij*solution.q[j];
  }

  return ret;
}


void GAUSSELIMINATE(double **E,int M,double *RHS)//适用于全带宽存储阵
{
  int i;
  int k;
  for(k=0;k<M-1;k++)
  {
    // 选主元
    double bmax=0.0;
    int ik;
    for(i=k;i<M;i++)
    {
      if(bmax<fabs(E[i][k]))
      {
        bmax=fabs(E[i][k]);
        ik=i;
      }
    }
    if(bmax<1.0e-10) 
    {
      cout<<"主元太小"<<endl;
      //system("pause"); // 主元太小
    }
    // 交换第ik行和第k行的元素
    if(ik!=k)
    {
      double t;
      for(i=k;i<M;i++)
      {
        t=E[ik][i];
        E[ik][i]=E[k][i];
        E[k][i]=t;
      }
      t=RHS[ik];
      RHS[ik]=RHS[k];
      RHS[k]=t;
    }
    // 消元
    for(i=k+1;i<M;i++)
    {
      if(E[i][k]!=0)
      {
        double lk=E[i][k]/E[k][k];
        int j;
        for(j=k+1;j<M;j++)
        {
          E[i][j]=E[i][j]-lk*E[k][j];
        }
        RHS[i]=RHS[i]-lk*RHS[k];
      }
    }
  }
  if(fabs(E[M-1][M-1])<1.0e-10)
  {
    cout<<"主元太小"<<endl;
    //system("pause"); // 主元太小
  }

  // 消元法结束后开始回代
  RHS[M-1]=RHS[M-1]/E[M-1][M-1];
  for(i=M-2;i>=0;i--)
  {
    double s=0.0;
    int j;
    for(j=i+1;j<M;j++)
    {
      s=s+E[i][j]*RHS[j];
    }
    RHS[i]=(RHS[i]-s)/E[i][i];
  }
}

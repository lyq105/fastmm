#include "moment.h"
#include "initmesh.h"
#include <complex.h>
#include <math.h>

//#define _Complex complex;


template <typename _T>
_T max (_T a, _T b)
{
	return a>b?a:b;
}	
template <typename _T>
_T min (_T a, _T b)
{
	return a<b?a:b;
}	


// 计算叶子节点的多极矩
void cal_multipole_moment(Quadtree& qtree,int cellindex,Mesh mesh,double* data)
{
	//  计算多极矩
	//  遍历叶子节点的单元，计算积分，并累加到多极系数中
	int eindex;
	double x1,x2,y1,y2;
	double phi,q;
	_Complex double z1,z2,zwbar,zp1,zp2, znorm;

	//  查询单元序号并开始遍历
	QuadtreeNode& qtnode = *qtree.treeNodeList[cellindex];
	int nexp = qtnode.mpcoeff.terms;

	int startIndex = qtnode.startIndex;
	int endIndex = qtnode.maxElem - 1;
	
	for (int i = startIndex; i <= endIndex; ++i)
	{
		// 取出单元号
		eindex = qtree.elemList[i];

		x1 = mesh.node[ mesh.elem[eindex][0] ] [0];
		x2 = mesh.node[ mesh.elem[eindex][0] ] [1];
		y1 = mesh.node[ mesh.elem[eindex][1] ] [0];
		y2 = mesh.node[ mesh.elem[eindex][1] ] [1];
		
		z1 = (x1 - qtnode.center.x[0]) + (x2 - qtnode.center.x[1]);
		z2 = (y1 - qtnode.center.x[0]) + (y2 - qtnode.center.x[1]);
		zwbar = conj(z2 - z1);             
		zwbar = zwbar/cabs(zwbar);
		zp1   = z1*zwbar;                   
		zp2   = z2*zwbar;                   
		znorm = (mesh.elemnorm[eindex][0]) + (mesh.elemnorm[eindex][1])*I;
		//znorm = cmplx(dnorm(1,nelm),dnorm(2,nelm));
		//
		// 指定迭代变量 
    
		if (mesh.boundaryCondition.bctags[i] == 1)
		{
			phi = 0;
			q = data[i];
		}
		if (mesh.boundaryCondition.bctags[i] == 2)
		{
			phi = data[i];
			q = 0;
		}

		// G kernel
		qtnode.mpcoeff.mp_data[0] = qtnode.mpcoeff.mp_data[0] - (zp2-zp1)*q;     
		for(int k = 1; k < nexp; k++)
		{
	 		// F kernel
			qtnode.mpcoeff.mp_data[k] = qtnode.mpcoeff.mp_data[k] + (zp2-zp1)*znorm*phi; 
			zp1  = zp1*z1/(k+1);
			zp2  = zp2*z2/(k+1);
			// G kernel
			qtnode.mpcoeff.mp_data[k] = qtnode.mpcoeff.mp_data[k] - (zp2-zp1)*q;
		}
	}// 遍历单元
} 
// 传递多极矩到父节点
void transfer_mm_to_its_father_by_m2m(Quadtree& qtree,int cellindex)
{
	_Complex double z0,z1;
	QuadtreeNode& qnode = *qtree.treeNodeList[cellindex]; 
	int fatherIndex = qnode.father;
	int nexp = qnode.mpcoeff.terms;

	QuadtreeNode& fnode = *qtree.treeNodeList[fatherIndex]; 
	// 计算节点的积分
	Point& qcenter = qnode.center;
	Point& fcenter = fnode.center;

	z0 = (fcenter.x[0] - qcenter.x[0]) + (fcenter.x[1] - qcenter.x[1])*I;

	z1 = 1; 
	for (int k = 0; k < nexp; k++)
	{
		for (int m = k; m< nexp; k++)
		{
			fnode.mpcoeff.mp_data[m] = fnode.mpcoeff.mp_data[m] + z1*qnode.mpcoeff.mp_data[m-k];
		}
		z1 = z1*z0/(k+1);
	}
}
// 将交互节点的多极矩转化为中心节点的局部矩 M2L
void transfer_lm_by_m2l(Quadtree& qtree,int cIndex, int iIndex)
{
	_Complex double z0,z1;
	signed int sgn;
	QuadtreeNode& qnode = *qtree.treeNodeList[cIndex];
	QuadtreeNode& inode = *qtree.treeNodeList[iIndex];
 	Point& qcenter = qnode.center;	// center cell center point. 
 	Point& icenter = inode.center;	// interact cell center point. 
	int nexp = qnode.mpcoeff.terms;
	int ntylr = qnode.lccoeff.terms;

	z0 = (qcenter.x[0] - icenter.x[0]) + (qcenter.x[1] - icenter.x[1])*I;

	qnode.lccoeff.lc_data[0] -= clog(z0)* inode.mpcoeff.mp_data[0]; 
	z1 = 1;

	int length = qnode.mpcoeff.terms + qnode.lccoeff.terms;
	for (int lIndex = 1; lIndex < length; lIndex++)
	{
		z1 = z1/z0;
		int kmin = max (0, lIndex - nexp);
		int kmax = min (lIndex, ntylr);

		sgn = pow(-1,kmin);

		for (int kIndex = kmin; kIndex < kmax; kIndex++)
		{
			qnode.lccoeff.lc_data[lIndex] = qnode.lccoeff.lc_data[lIndex]
			 + sgn*z1*inode.mpcoeff.mp_data[lIndex - kIndex];
			sgn*=-1;
		}	
		z1 = z1 * lIndex; 
	}
}	
// 将父节点的局部矩转换为当前节点的局部矩 L2L
void transfer_lm_by_l2l(Quadtree& qtree, int cIndex)
{
	_Complex double z0,z1;
	QuadtreeNode& qnode = *qtree.treeNodeList[cIndex];
	QuadtreeNode& fnode = *qtree.treeNodeList[qnode.father];
 	Point& qcenter = qnode.center;	// center cell center point. 
 	Point& fcenter = fnode.center;	// interact cell center point. 

	int ntylr = qnode.lccoeff.terms;

	z0 = (qcenter.x[0] - fcenter.x[0]) + (qcenter.x[1] - fcenter.x[1])*I;
	z1 = 1;
 	
	for (int k = 0; k < ntylr; k++)
	{
		for (int m = 0; m < ntylr - k; m++)
		{
			qnode.lccoeff.lc_data[m] = qnode.lccoeff.lc_data[m] + z1 * fnode.lccoeff.lc_data[m];
		}
		z1 = z1*z0/(k+1);
	}	
}
//  将局部矩传递到计算点pt处
void transfer_lm_to_local_point(Quadtree& qtree, int cIndex, Point& pt)
{
	_Complex double z0,z1;
	QuadtreeNode& qnode = *qtree.treeNodeList[cIndex];
 	Point& qcenter = qnode.center;	// center cell center point. 
	int ntylr = qnode.lccoeff.terms;

	z0 = (pt.x[0] - qcenter.x[0]) + ( pt.x[1] - qcenter.x[1])*I;
	z1 = 1;
 	
	for (int k = 0; k < ntylr; k++)
	{
		for (int m = 0; m < ntylr - k; m++)
		{
			//qnode.lccoeff.lc_data[m] = qnode.lccoeff.lc_data[m] + z1 * fnode.lccoeff.lc_data[m];
		}
		z1 = z1*z0/(k+1);
	}		
}

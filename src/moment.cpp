#include "moment.h"
#include <complex.h>

#define _Complex complex;
// 计算叶子节点的多极矩
void cal_multipole_moment(Quadtree& qtree,int cellindex,Mesh mesh,double* data)
{
	//  计算多极矩
	//  遍历叶子节点的单元，计算积分，并累加到多极系数中
	int eindex;
	double x1,x2,y1,y2;
	complex double z1,z2,zwbar,zp1,zp2;

	//  查询单元序号并开始遍历
	QuadtreeNode& qtnode = *qtree.treeNodeList[cellindex];

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
		zwbar = zwbar/abs(zwbar);
		zp1   = z1*zwbar;                   
		zp2   = z2*zwbar;                   
		znorm = mesh.elemnorm[0] + mesh.elemnorm[1];
		//znorm = cmplx(dnorm(1,nelm),dnorm(2,nelm)); 
		if(bc(1,nelm) .eq. 1.d0) then      
			phi = 0.D0
				q   = u(i)
		else if(bc(1,nelm) .eq. 2.d0) then
			phi = u(i)
				q   = 0.D0
				endif


				a(0) = a(0) - (zp2-zp1)*q;             // G kernel
				for(int k=1; k< nexp;k++)
					a(k) = a(k) + (zp2-zp1)*znorm*phi    // F kernel
						zp1  = zp1*z1/(k+1)
						zp2  = zp2*z2/(k+1)
						a(k) = a(k) - (zp2-zp1)*q;          // G kernel
						enddo





	}// 遍历单元
} 
// 传递多极矩到父节点
void transfer_mm_to_its_father_by_m2m(Quadtree& qtree,int cellindex)
{
	int fatherIndex = qtree.treeNodeList[cellIndex] -> father;
	// 计算节点的积分
}
// 将交互节点的多极矩转化为中心节点的局部矩
void transfer_lm_by_m2l(Quadtree& qtree,int cellIndex)
{
}	
// 将父节点的局部矩转换为当前节点的局部矩
void transfer_lm_by_l2l(Quadtree& qtree, int cellIndex)
{
}
//  将局部矩传递到计算点处
void transfer_lm_to_local_point()
{

}

#include "moment.h"

// 计算叶子节点的多极矩
void cal_multipole_moment(Quadtree& qtree,int cellindex,Mesh mesh,double* data)
{
	//  计算多极矩
	//  遍历叶子节点的单元，计算积分，并累加到多极系数中
	
	//  查询单元序号并开始遍历
	QuadtreeNode& qtnode = *qtree.treeNodeList[cellindex];

	int startIndex = qtnode.startIndex;
	int endIndex = qtnode.maxElem - 1;
	
	for (int i = startIndex; i <= endIndex; ++i)
	{
		// 取出单元的两个节点列表
		mesh.CellList[i]
		int n =
		z1 = cmplx(y(1,n1)-cx, y(2,n1)-cy) 
		z2 = cmplx(y(1,n2)-cx, y(2,n2)-cy) 
		zwbar = conjg(z2 - z1)             
		zwbar = zwbar/abs(zwbar)           ! omega bar
		zp1   = z1*zwbar                   
		zp2   = z2*zwbar                   
		znorm = cmplx(dnorm(1,nelm),dnorm(2,nelm)) ! complex normal n
		if(bc(1,nelm) .eq. 1.d0) then      ! Assign values to phi and q
			phi = 0.D0
				q   = u(i)
		else if(bc(1,nelm) .eq. 2.d0) then
			phi = u(i)
				q   = 0.D0
				endif


				a(0) = a(0) - (zp2-zp1)*q              ! G kernel
				do k=1,nexp
					a(k) = a(k) + (zp2-zp1)*znorm*phi    ! F kernel
						zp1  = zp1*z1/(k+1)
						zp2  = zp2*z2/(k+1)
						a(k) = a(k) - (zp2-zp1)*q            ! G kernel
						enddo

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

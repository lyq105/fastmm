/// This file is preform matrix vector products of Fast mulitipole boundary method.
#include "matvec.h"
#include <complex.h>

const double pi2 = 2 * 3.14159265; 

int find_elem(Quadtree& qtree,int cIndex, unsigned int* list)
{
	QuadtreeNode& qnode = *qtree.treeNodeList[cIndex];
	
	// 计算向量长度
	int len = qnode.maxElem;
	for (int i = 0; i <qnode.lenInterList; i++) 
	{
		len += qtree.treeNodeList[ qnode.interList[i] ].maxElem;
	}
	
	list = new unsigned int [len + 1];

	list[0] = len;
	int index = 1;
	for (int i = 0; i <qnode.lenInterList; i++) 
	{
		int maxelem = qtree.treeNodeList[ qnode.interList[i] ].maxElem;
		int start = qtree.treeNodeList[ qnode.interList[i] ].startIndex;
		int end = start + qtree.treeNodeList[ qnode.interList[i] ].maxElem - 1;	
		for (int j = start; j <= end; j++)
	 	{
			list[index] = qtree.elemList[j];
			index ++;
		}
	}

	return 0;
}




void matvec(Quadtree& qtree,Mesh mesh,double* itervec, double* data)
{
	// 构造矩阵向量乘法主要来自于两个部分
	//
	//  
	//  这一个函数的输出是一个向量，这一个向量的每一个分量是对应于配置点的
	//  在四叉树的前提下，搜索的过程是基于遍历叶子节点的
	// for (int i =0; i< vec_length; i++) // 对于每一个点计算
	// 第一个部分是临近单元和本身所在的单元的贡献，采用直接数值积分
	//
	// 数值积分计算积分贡献的地方是在
	// 算法过程是这样的。
	// 首先是遍历树的叶子节点。 在每一个叶子节点中遍历配置点 
	// 针对每一个配置点计算积分

	// 更新多极和局部系数
	_Complex double z0,z1;

	quadtree_upward(qtree,mesh,data); //多极系数
	quadtree_downward(qtree); // 局部系数

	for (int i = 0; i < qtree.numberTreenode; i++) 
	{
		if (qtree.treeNodeList[i] ->isLeaf == 1)
		{
			QuadtreeNode& qnode = qtree.treeNodeList[i];

			int start = qnode.startIndex; 
			int end = qnode.startIndex + qnode.maxElem - 1;

			for (int index = start; index <= end; index++)  // d单元循环，实际上是配置点的循环
			{
				// 第一部分是直接积分
				// 以单元的中点为配置点计算积分
				// 1---  计算配置点坐标
				int eindex = qtree.elemList[index];
				int nindex0 = mesh.elem[eindex][0];
				int nindex1 = mesh.elem[eindex][1];
				double source_x0 = 0.5*(bmesh.node[nindex0][0] + bmesh.node[nindex1][0]);
				double source_x1 = 0.5*(bmesh.node[nindex0][1] + bmesh.node[nindex1][1]);

				// 2--- 查找需要直接计算单元积分的单元号。

				unsigned int* elemNumberList;

				find_elem(qtree,i,elemNumberList);	
				for (int j = 1; j < elemNumberList[0]+1; j++) // 遍历当前树节点中的单元
				{
					// find out the number of 2 ends of a cell.
					int cindex = elemNumberList[j]
					int	nindex0 = bmesh.elem[cindex][0];
					int nindex1 = bmesh.elem[cindex][1];
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

					if(eindex == cindex) // diagonal element of matrix.
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

					if(mesh.boundaryCondition.bctags[cindex] == 1)
					{
						itervec[eindex] += -gij*data[eindex];
					}
					else if (mesh.boundaryCondition.bctags[cindex] == 2)
					{
						itervec[eindex] += hij*data[eindex];
					}
				}//loop for source point.


				// 第二部分是临近单元和远程单元的贡献，采用局部展开式和局部系数计算

				double fact = 1;
				for(int itylr=0; itylr < ntylr;itylr ++)
				{
					fact = fact/itylr;
					qnode.lccoeff.lc_data[itylr] *= fact;
				}
				zp = qnode.lccoeff.lc_data[itylr];
				z0 = (source_x0 - qnode.center[0]) + (source_x1 - qnode.center[1])*I;
				for(int itylr=ntylr-1; itylr >= 0; itylr--)
				{
					zp = zp*z0 + b(itylr,icell);
				}	
				zp = zp/pi2;
				itervec[eindex] += creal(zp);
			}
		}
	}

}

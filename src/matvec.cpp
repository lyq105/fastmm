/// This file is preform matrix vector products of Fast mulitipole boundary method.
#include "matvec.h"

void matvec()
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
	Quadtree& qtree;
	
	// 更新多极和局部系数
	
	quadtree_upward(qtree,mesh,data); //多极系数
	quadtree_downward(qtree); // 局部系数

	for (int i = 0; i < qtree.numberTreenode; i++) 
	{
		if (qtree.treeNodeList[i] ->isLeaf == 1)
		{
			QuadtreeNode& qnode = qtree.treeNodeList[i];

			int start = qnode.startIndex; 
		 	int end = qnode.startIndex + qnode.maxElem - 1;

			for (int eindex = start; eindex < end; eindex++)
		 	{
			// 第一部分是直接积分
			// 以单元的中点为配置点计算积分
			
			//// 当前树节点中的单元积分

			//// 邻接树节点的单元积分	


			// 第二部分是临近单元和远程单元的贡献，采用局部展开式和局部系数计算
			
			// 第一步计算当前配置点处的局部展开系数
		

			// 第二步利用多极系数计算远程单元和单元的贡献	
			
				
			}


			

		}
	}

}

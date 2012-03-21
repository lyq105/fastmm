#include "moment.h"

// 计算叶子节点的多极矩
void cal_multipole_moment(Quadtree& qtree,int cellindex,Mesh mesh,double* data)
{
	//  计算多极矩
	//  遍历叶子节点的单元，计算积分，并累加到多极系数中
	
	//  查询单元序号并开始遍历
	
	for ()
} 
// 传递多极矩到父节点
void transfer_mm_to_its_father_by_m2m(Quadtree& qtree,int cellindex)
{
	int fatherIndex = qtree.treeNodeList[cellIndex] -> father;
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

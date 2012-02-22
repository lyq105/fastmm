#include "moment.h"

// 计算多极矩
void cal_multipole_moment(Quadtree& qtree,int cellindex,Mesh mesh,double* data)
{
} 
// 传递多极矩到父节点
void transfer_mm_to_its_father_by_m2m(Quadtree& qtree,int cellindex)
{
}
// 将交互节点的多极矩转化为中心节点的局部矩
void transfer_lm_by_m2l(Quadtree& qtree,int cellIndex)
{
}	
// 将父节点的局部矩转换为当前节点的局部矩
void transfer_lm_by_l2l(Quadtree& qtree, int cellIndex)
{
}


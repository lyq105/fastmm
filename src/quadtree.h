//=============================================================================
//   Project name:  head of quadtree.cpp
//       Filename:  quadtree.h
//   Descriptions:  四叉树的头文件
//         Author:  liyiqiang(lyq105@163.com)
//         Create:  2011-06-22 21:16:26
//      Last edit:  2011-06-22 21:16:26
//=============================================================================

#ifndef __QUADTREE_H__ 
#define __QUADTREE_H__ 
#include "defdata.h"

// data type's Defination

typedef struct __multipole_coeff__
{
	int treeNodenumber;  // 树节点的数目
	int terms;           // 展开项数
	double** mp_data;    // 多极系数  
} MPcoeff;  // 多极矩系数


typedef struct __local_coeff__
{
	int treeNodenumber;  // 树节点的数目
	int terms;           // 展开项数
	double** lc_Data;    // 多极系数  p*num_treenode 
} LCcoeff;  // 局部矩系数

typedef struct _quadtreenode_type
{
	int isLeaf;     // if the node is a leaf;
	int level;      // tree node level.  
	int coord;      // coordinate in background tree;
	int father;     // tree node's father;
	int number;     // numbering of tree node;
	int startIndex; // start index of element number in the element list
	int maxElem;    // max number of boundary elements in this tree node; 
  int interList[27]; // numbering of tree node's interact list;
  int lenInterList; // number of interact node.
  int adjacentList[8];// 
  int lenAdjList;
	double length;  // tree node length;
	Point center;   // center of tree node;
	LCcoeff* lccoeff;
	MPcoeff* mpcoeff; // 
}QuadtreeNode;

typedef struct _quadtree_type
{
	int numberTreenode;   // Number of tree node; 
	int numElem;          // Number of element;
	int treeDepth;        // Depth of tree;
	int* levelCell;    // position of first level l cell; 
	int* numlevelCell; // how many non empty cell of level l;
	int* elemList;    // boundary elements in each tree level, max tree level is 20; 
	QuadtreeNode** treeNodeList; 
}Quadtree;


int quadtree_creat(Quadtree& qtree, Mesh& mesh); 

int quadtree_creat_childs(Mesh mesh,Quadtree& qtree,QuadtreeNode* ftnode);
int quadtree_init(Quadtree& qtree, int numelem); 
int quadtree_insert(Quadtree& qtree, QuadtreeNode* tnode);
int quadtree_destory(Quadtree& qtree);
int quadtree_upward(Quadtree& qtree);
int quadtree_downward(Quadtree& qtree);

int find_interact_list(Quadtree& qtree);
int is_adjacent(QuadtreeNode& qna, QuadtreeNode& qnb);

int plot_tree_and_list(char* filename, Quadtree qtree, int number);
int cal_convex_hall(Mesh& mesh,Point& center,double& length);
int plot_quadtree(Quadtree qtree,char filename[]);
int print_quadtree_info(Quadtree qtree,char filename[]);
int animate_it(Quadtree qtree);
int print_interact_list(Quadtree qtree,int number);

#endif

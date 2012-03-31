//=============================================================================
//   Project name:  libtree
//       Filename:  quadtree.cpp
//   Descriptions:  constrction a quad tree code.
//         Author:  liyiqiang(lyq105@163.com)
//         Create:  2011-06-22 21:16:26
//      Last edit:  2011-06-22 21:16:26
//=============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "quadtree.h"
//#include "defdata.h"
#include "moment.h"
#include <iostream>
#include <math.h>
using namespace std;

// 声明为static，使得变量为全局变量可以在程序中更改；

static int maxTreeDepth = 50; // 数的最大深度
static int maxEleminCell = 2; // 树节点中的最大单元数
static int main_number = 0;   // 树节点的初始编号
static int mmExpantionTerms = 10; //  多极展开项树


int animate_it(Quadtree qtree)
{
  char filename[200]; 
  char filename2[200]; 

  FILE* fp;
  for (int i = 4; i < qtree.numberTreenode; i++) 
  {
    sprintf(filename,"inter%d.dat",i);
    //cout << filename << endl;
    plot_tree_and_list(filename,qtree,i);
    
    sprintf(filename2,"inter%d.asy",i);
    fp = fopen(filename2,"w");
    fprintf(fp, "import plot_frame;\n");
    fprintf(fp, "plot_frame(\"%s\");",filename);
    fclose(fp);
  }
  return 0;
}
/// generater a tree graph using dot language;
int draw_tree_graph(Quadtree& qtree)
{
	FILE* fp;
	fp = fopen("tree_graph.dot","w");
	fprintf(fp, "digraph fmm_bem\{\n");
	fprintf(fp, "ranksep = 2;\n");
	fprintf(fp, "node[shape = box]\n");  // tree node's style
	// write tree node infomation.
	fprintf(fp, "node0 [label= \"root\", width = 2, height = 1.7,fontsize = 28]\n");
	for (int nIndex = 1; nIndex < qtree.numberTreenode; nIndex++)
 	{
		if (qtree.treeNodeList[nIndex] ->isLeaf ==1)
			fprintf(fp, "node%d [label= \"%d\", shape = circle, fillcolor = lightgrey, style = filled]\n",nIndex,nIndex);
		fprintf(fp, "node%d [label= \"%d\"]\n",nIndex,nIndex);
	}

	// specify the edge style.
	fprintf(fp, "edge[color = black]\n");
	
	// write edge infomation.
	for (int nIndex = 1; nIndex < qtree.numberTreenode; nIndex++)
	{
		int father = qtree.treeNodeList[nIndex] ->father;
		fprintf(fp, "node%d -> node%d\n",father,nIndex);		
	}

	fprintf(fp, "\}\n");

	fclose(fp);

	system("dot tree_graph.dot -Tpdf -o tree_graph.pdf");

	return 0;
}

int plot_tree_and_list(char* filename, Quadtree qtree, int number)
{
  FILE* fp;
	double coef_x[4]={-0.5,0.5,0.5,-0.5};
	double coef_y[4]={-0.5,-0.5,0.5,0.5};
  fp =fopen(filename,"w");
  int lenil = qtree.treeNodeList[number] ->lenInterList;
  fprintf(fp,"%d\n",lenil + 1);
  // draw central cell
  for(int j =0; j<4;j++)
  {
    double x = qtree.treeNodeList[number]->center.x[0] + coef_x[j] * qtree.treeNodeList[number]->length;
    double y = qtree.treeNodeList[number]->center.x[1] + coef_y[j] * qtree.treeNodeList[number]->length;
    fprintf(fp,"%12.8f %12.8f\n",x,y);
  }
  for (int inlenal = 0; inlenal < lenil; inlenal++) 
  {
    int index = qtree.treeNodeList[number] -> interList[inlenal]; 
    for(int j =0; j<4;j++)
    {
      double x = qtree.treeNodeList[index]->center.x[0] + coef_x[j] * qtree.treeNodeList[index]->length;
      double y = qtree.treeNodeList[index]->center.x[1] + coef_y[j] * qtree.treeNodeList[index]->length;
      fprintf(fp,"%12.8f %12.8f\n",x,y);
    }
  }
  // draw adj cell
  int lenadj = qtree.treeNodeList[number] ->lenAdjList;
  fprintf(fp,"%d\n",lenadj);
  printf("%d\n",lenadj);

  for (int inlenal = 0; inlenal < lenadj; inlenal++) 
  {
    int index = qtree.treeNodeList[number] ->adjacentList [inlenal]; 
    for(int j =0; j<4;j++)
    {
      double x = qtree.treeNodeList[index]->center.x[0] + coef_x[j] * qtree.treeNodeList[index]->length;
      double y = qtree.treeNodeList[index]->center.x[1] + coef_y[j] * qtree.treeNodeList[index]->length;
      fprintf(fp,"%12.8f %12.8f\n",x,y);
    }
  }
  // draw background quad tree.
  fprintf(fp,"%d\n",qtree.numberTreenode);
  for(int i=0; i<qtree.numberTreenode; i++)
	{
		for(int j =0; j<4;j++)
		{
			double x = qtree.treeNodeList[i]->center.x[0] + coef_x[j] * qtree.treeNodeList[i]->length;
			double y = qtree.treeNodeList[i]->center.x[1] + coef_y[j] * qtree.treeNodeList[i]->length;
			fprintf(fp,"%12.8f %12.8f\n",x,y);
		}
	}
  fclose(fp);

  return 0;
}

int print_interact_list(Quadtree qtree,int number)
{
  FILE* fp;
	double coef_x[4]={-0.5,0.5,0.5,-0.5};
	double coef_y[4]={-0.5,-0.5,0.5,0.5};
  char filename[100] = "bbbb.plt";
	fp =fopen(filename,"w");
  int lenil = qtree.treeNodeList[number] ->lenInterList;
  //int lenil = qtree.treeNodeList[number] ->lenAdjList;
	fprintf(fp,"zone N=%d,E=%d,T=\"quadtree2\",F=FEPOINT,ET=QUADRILATERAL\n",4*qtree.numberTreenode,lenil+1);
	for(int i=0; i<qtree.numberTreenode; i++)
	{
		for(int j =0; j<4;j++)
		{
			double x = qtree.treeNodeList[i]->center.x[0] + coef_x[j] * qtree.treeNodeList[i]->length;
			double y = qtree.treeNodeList[i]->center.x[1] + coef_y[j] * qtree.treeNodeList[i]->length;
			fprintf(fp,"%12.8f %12.8f\n",x,y);
		}
	}
	//for(int i=0; i<qtree.numberTreenode; i++)
  for (int inlenal = 0; inlenal < lenil; inlenal++) 
  {
    //int index = qtree.treeNodeList[number] -> adjacentList[inlenal]; 
    int index = qtree.treeNodeList[number] -> interList[inlenal]; 
    fprintf(fp, "%d %d %d %d\n",4*index+1,4*index+2,4*index+3,4*index+4);
  }
  fprintf(fp, "%d %d %d %d\n",4*number+1,4*number+2,4*number+3,4*number+4);
  fclose(fp);


  //int boolsm = is_adjacent(*(qtree.treeNodeList[241]),*(qtree.treeNodeList[240]));
  //cout << "isisi = " << boolsm << endl;
  return 0;
}

// 对一个树节点生成交互节点列表,和相邻节点列表，从同层次里面查找。
//      TODO: 当前算法只是遍历一遍同层次的树节点，可以优化一下。
int find_interact_list(Quadtree& qtree)
{
	//print_quadtree_info(qtree,"treeinfo2.txt");
  // 树节点的循环
  int firstLevel2cell = qtree.levelCell[2];
  
  for (int i = firstLevel2cell; i < qtree.numberTreenode; i++) 
  {
    // 遍历同层次的树节点
    int lev = qtree.treeNodeList[i] -> level;  // get depth of current node.

    int startIndex = qtree.levelCell[lev];
    int endIndex = startIndex + qtree.numlevelCell[lev] - 1;

    //cout << i << " === "
         //<< qtree.treeNodeList[i] -> level << " -- "
         //<< startIndex << " -- "
         //<< endIndex  << " -- "
         //<< endl;
    int fatherIndex = qtree.treeNodeList[i] -> father;
    int adjNumber = 0;
    int interactNumber = 0;
    for (int j = startIndex; j <= endIndex; ++j)  // 同层树节点的循环
    {
      // 计算树节点之间的距离来判断是否是交互节点
      // 
      // 先判断父节点是否相邻，在判断距离。  TODO TODO
      // 判断两个节点的父节点的距离。
      int testfatherIndex = qtree.treeNodeList[j] -> father;

      // 是否与当前节点相邻，如果相邻则列入相邻节点列表，否则加入交互节点列表。
      if (is_adjacent(*(qtree.treeNodeList[i]),*(qtree.treeNodeList[j])) )
      {
        qtree.treeNodeList[i] -> adjacentList[adjNumber] = j;
        adjNumber ++ ; 
      }
      else
      {
        // 如果父节点相邻
        if (is_adjacent(*(qtree.treeNodeList[fatherIndex]), *(qtree.treeNodeList[testfatherIndex])))
        {
          qtree.treeNodeList[i] -> interList[interactNumber] = j;
          interactNumber ++;
        }
      }// 是否相邻
    } // 遍历同层次的树节点完毕
    qtree.treeNodeList[i] -> lenInterList = interactNumber;
    qtree.treeNodeList[i] -> lenAdjList = adjNumber;
  } // 所有树节点计算完毕
  //print_interact_list(qtree,241);
  //plot_tree_and_list(qtree,155);
  //animate_it(qtree);
  return 0;
}

int is_adjacent(QuadtreeNode& qna, QuadtreeNode& qnb)
{
  double distant = (qna.center.x[0] - qnb.center.x[0])*(qna.center.x[0] - qnb.center.x[0])
    + (qna.center.x[1] - qnb.center.x[1])*(qna.center.x[1] - qnb.center.x[1]);
  distant = sqrt(distant);
  if ( fabs((qna.length + qnb.length) - 2*distant) < 1e-9  // 正对着
      || fabs((qna.length + qnb.length)*sqrt(2.0) - 2*distant) < 1e-9 ) // 斜对着
  {
    return 1;
  } 
  return 0;
}


// 计算网格点的一个正方形凸包
int cal_convex_hall(Mesh& mesh,Point& center,double& length)
{
	int i,j;
	double x_max[2],x_min[2];
	length = 0;
	for (j = 0; j < 2; j++)
	{
		x_max[j] = x_min[j]=0;
		//for (i = 0; i < mesh.TotleNodeNumber; i++)
		//{
			//x_max[j] = x_max[j] >= mesh.MeshPointList[i].x[j]? x_max[j] : mesh.MeshPointList[i].x[j];
			//x_min[j] = x_min[j] <= mesh.MeshPointList[i].x[j]? x_min[j] : mesh.MeshPointList[i].x[j];
		//}
		for (i = 0; i < mesh.nodenum; i++)
		{
			x_max[j] = x_max[j] >= mesh.node[i][j]? x_max[j] : mesh.node[i][j];
			x_min[j] = x_min[j] <= mesh.node[i][j]? x_min[j] : mesh.node[i][j];
		}
		center.x[j] = 0.5*(x_max[j] + x_min[j]);
		length = length >= (x_max[j] - x_min[j]) ? length: (x_max[j] - x_min[j]);
	}
	length *= 1.1;
	return 0;
}
// 创建四叉树
int quadtree_creat(Quadtree& qtree, Mesh& mesh)
{
	int i;
	int depth, nodenumber;
	int startIndex, endIndex;

	QuadtreeNode* root = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));

	cal_convex_hall(mesh,root->center,root->length);


	root -> number = 0;    // root 's global numbering
	root -> level = 0;     // root 's level
	root -> isLeaf = 0;    // is a leaf node
	root -> father = -1;   // root father's numbering is -1 
	root -> startIndex = 0; 
	root -> maxElem = mesh.elemnum;
	root -> coord = 0;
	
	// 初始化树
	quadtree_init(qtree,root->maxElem);
	// 插入根结点
	quadtree_insert(qtree,root);

	qtree.levelCell[0] = 0;
	qtree.numlevelCell[0] = 1;

	if (root -> maxElem <= maxEleminCell)
	{
		return 1;
	}

	// 逐层遍历生成树结构

	for (depth = 1; depth < maxTreeDepth; depth++)
	{
		nodenumber = qtree.numberTreenode;
		// 遍历depth-1层的树节点
		startIndex = qtree.levelCell[depth - 1] ; // depth-1层节点开始的下标
		endIndex = startIndex + qtree.numlevelCell[depth - 1] - 1; // depth-1层节点结束的下标

		for (i = startIndex; i <= endIndex; i++)
		{
			if (qtree.treeNodeList[i] -> maxElem > maxEleminCell)
			{
				qtree.treeNodeList[i] -> isLeaf = 0;
				quadtree_creat_childs(mesh,qtree,qtree.treeNodeList[i]); // 四分节点，并将非空树节点插入树节点列表
			}
			else
				qtree.treeNodeList[i] -> isLeaf = 1;
		}

		if (qtree.numberTreenode != nodenumber)	 // 如果有新节点产生,
			{
				qtree.levelCell[depth] = endIndex + 1;
				qtree.numlevelCell[depth] = qtree.numberTreenode - nodenumber;
			}
		else break;
	}
	qtree.treeDepth = depth;

  find_interact_list(qtree);

	//print_quadtree_info(qtree,"treeinfo1.txt");	
	plot_quadtree(qtree,"treeinfo.plt");

	return 0;
}

// 四分树节点
//
//      2----3
//      |    |
//      0----1
//
int quadtree_creat_childs(Mesh mesh,Quadtree& qtree,QuadtreeNode* ftnode)
{
	int i,j;
	int pos[4];
	int startIndex,endIndex,eindex,estartIndex;

	// 分配一个和父节点单元数一样的临时矩阵
	int** temparray = new int*[4];
	for (i = 0; i < 4; i++) {
		temparray[i] = new int[ftnode -> maxElem];
	}

	double coef_x[4]={-0.25,0.25,0.25,-0.25};
	double coef_y[4]={-0.25,-0.25,0.25,0.25};

	// 计算每一个象限的单元数
	for (i = 0; i < 4; i++) {
		pos[i] = 0;
	}

	// startIndex和endIndex表示在树的单元列表中开始和结束索引。
	startIndex = ftnode -> startIndex;
	endIndex = startIndex + ftnode -> maxElem;


	for(i = startIndex; i< endIndex; i++)  //  遍历父节点中的单元
	{
		// 计算单元中心的坐标
		eindex = qtree.elemList[i];

		double x1[2],x2[2];
/*
		x1[0] = mesh.MeshPointList[mesh.CellList[eindex].cellnumber[0]].x[0];
		x1[1] = mesh.MeshPointList[mesh.CellList[eindex].cellnumber[0]].x[1];

		x2[0] = mesh.MeshPointList[mesh.CellList[eindex].cellnumber[1]].x[0];
		x2[1] = mesh.MeshPointList[mesh.CellList[eindex].cellnumber[1]].x[1];
*/
		x1[0] = mesh.node[ mesh.elem[eindex][0] ][0];
		x1[1] = mesh.node[ mesh.elem[eindex][0] ][1];
		x2[0] = mesh.node[ mesh.elem[eindex][1] ][0];
		x2[1] = mesh.node[ mesh.elem[eindex][1] ][1];

		double x = 0.5*(x1[0] + x2[0]);
		double y = 0.5*(x1[1] + x2[1]);


		// 计算四个区域中的单元数
		if( y >= ftnode -> center.x[1]){ // N
			if (x >= ftnode -> center.x[0]){ // E
				temparray[2][pos[2]] = eindex; // 将该单元编号存入临时矩阵
				pos[2]++;
			}
			else{  // W
				temparray[3][pos[3]] = eindex;
				pos[3]++;
			}
		}
		else{	// S
			if (x >= ftnode -> center.x[0]){ // E
				temparray[1][pos[1]] = eindex; // 将该单元编号存入临时矩阵
				pos[1]++;
			}
			else{ // W
				temparray[0][pos[0]] = eindex; // 将该单元编号存入临时矩阵
				pos[0]++;
			}
		}
	}

	estartIndex = ftnode->startIndex;
	for (j = 0; j < 4; j++)
	{
		if (pos[j] != 0 ) // 创建非空节点
		{
			QuadtreeNode* childnode = new QuadtreeNode;
			// 赋层数；
			childnode -> level = ftnode -> level + 1;
			childnode -> coord = ftnode -> coord * 4 + j;
			// 父节点
			childnode -> father = ftnode -> number;
			// 节点的整体编号；
			main_number += 1;
			childnode -> number = main_number;
			// 单元开始下标；
			childnode -> startIndex = estartIndex;
			estartIndex += pos[j];
			// 节点中的单元个数；
			childnode -> maxElem = pos[j];
			// 单元中心点；
			childnode -> center.x[0] = ftnode -> center.x[0] + coef_x[j]* ftnode -> length;
			childnode -> center.x[1] = ftnode -> center.x[1] + coef_y[j]* ftnode -> length;

			// 单元边长；
			childnode -> length = 0.5 * ftnode -> length;

			for (int k = childnode->startIndex;k<childnode->startIndex + pos[j]; k++)
			{
				qtree.elemList[k] = temparray[j][k - childnode->startIndex];
			}

			//  初始化树节点的多极系数，局部系数的内存
			childnode -> mpcoeff.terms = mmExpantionTerms;
			childnode -> lccoeff.terms = mmExpantionTerms;

			childnode -> mpcoeff.mp_data = new _Complex double [mmExpantionTerms];
			childnode -> lccoeff.lc_data = new _Complex double [mmExpantionTerms];

			// 将节点插入节点列表；
			quadtree_insert(qtree,childnode);
		}
	}

	for (i = 0; i < 4; i++)
	{
		delete [] temparray[i];
	}
	delete [] temparray;
	return 0;
}

// 初始化一个空树
int quadtree_init(Quadtree& qtree, int numelem)
{
	qtree.numberTreenode = 0;
	qtree.numElem = numelem;
	qtree.treeDepth = 0;

	qtree.levelCell = new int[maxTreeDepth];
	qtree.numlevelCell = new int[maxTreeDepth];

	qtree.elemList = new int[numelem];

	for (int i = 0; i < numelem; i++) {
		qtree.elemList[i] = i;
	}

	qtree.treeNodeList = NULL;
	return 0;
}

// 将节点插入节点列表
int quadtree_insert(Quadtree& qtree, QuadtreeNode* tnode)
{
	int i;
	qtree.numberTreenode += 1;
	QuadtreeNode** templist;
	templist = new QuadtreeNode*[qtree.numberTreenode];

	for (i = 0; i < qtree.numberTreenode - 1; i++) {
		templist[i] = qtree.treeNodeList[i];
	}
	templist[qtree.numberTreenode - 1] = tnode;

	// 删除旧的树节点列表
	delete [] qtree.treeNodeList;

	qtree.treeNodeList = templist;
	return 0;
}

// 释放树占据的空间
int quadtree_destory(Quadtree& qtree)
{
	int i;
	delete [] qtree.levelCell;
	delete [] qtree.numlevelCell;
	delete [] qtree.elemList;

	for (i = 0; i < qtree.numberTreenode; i++) {
		//delete [] qtree.treeNodeList[i]->elementList;
		delete [] qtree.treeNodeList[i] -> lccoeff.lc_data;
		delete [] qtree.treeNodeList[i] -> mpcoeff.mp_data;
		delete qtree.treeNodeList[i];
	}
	delete [] qtree.treeNodeList;

	return 0;
}

int print_quadtree_info(Quadtree qtree,char filename[])
{
	FILE* fp;
	fp = fopen(filename,"w");
	// Export tree nodes;
	fprintf(fp,"number isleaf  coord    father   level    center_x\tcenter_y\tlength\tmaxnum\n");
	for(int i=0; i<qtree.numberTreenode; i++)
	{
		fprintf(fp,"%8d %8d %8d %8d %8d %12.8f %12.8f %12.8f %d\n",
				qtree.treeNodeList[i]->number,
				qtree.treeNodeList[i]->isLeaf,
				qtree.treeNodeList[i]->coord,
				qtree.treeNodeList[i]->father,
				qtree.treeNodeList[i]->level,
				qtree.treeNodeList[i]->center.x[0],
				qtree.treeNodeList[i]->center.x[1],
				qtree.treeNodeList[i]->length,
				qtree.treeNodeList[i]->maxElem);
	}
	fclose(fp);
	return 1;
}

int plot_quadtree(Quadtree qtree,char filename[])
{
	FILE* fp;
	double coef_x[4]={-0.5,0.5,0.5,-0.5};
	double coef_y[4]={-0.5,-0.5,0.5,0.5};

	fp =fopen(filename,"w");
	fprintf(fp,"zone N=%d,E=%d,T=\"quadtree\",F=FEPOINT,ET=QUADRILATERAL\n",4*qtree.numberTreenode,qtree.numberTreenode);
	for(int i=0; i<qtree.numberTreenode; i++)
	{
		for(int j =0; j<4;j++)
		{
			double x = qtree.treeNodeList[i]->center.x[0] + coef_x[j] * qtree.treeNodeList[i]->length;
			double y = qtree.treeNodeList[i]->center.x[1] + coef_y[j] * qtree.treeNodeList[i]->length;
			fprintf(fp,"%12.8f %12.8f\n",x,y);
		}
	}
	for(int i=0; i<qtree.numberTreenode; i++)
	{
		fprintf(fp, "%d %d %d %d\n",4*i+1,4*i+2,4*i+3,4*i+4);
	}
	return 0;
}
// 上行遍历树结构 计算每一个非空树节点的多极矩
// NOTE:   这里的参数data就是迭代的向量，或者说是矩阵乘向量的向量
//
int quadtree_upward(Quadtree& qtree, Mesh mesh, double* data)
{
	int maxLevel = qtree.treeDepth - 1;  // retrive max tree level.

	// 逐层向根方向遍历树结构 从树的最后一层到第二层
	for (int depth = maxLevel; depth >= 2; --depth)  
	{
		// 获取非空树节点的下标

		int startIndex = qtree.levelCell[ depth ];
		int endIndex = startIndex + qtree.numlevelCell[ depth ] - 1;

		// 遍历同层次的树节点
		for (int cellIndex = startIndex; cellIndex <= endIndex; ++cellIndex) 
		{
			// 这里需要查询父节点的编号  
			// 在传递过程中查询，并传递到父节点
			//{
				//int fatherIndex = qtree.treeNodeList[cellIndex] -> father;
			//}

			// 如果树节点不是叶子节点，则直接使用M2M将其累加到父节点中；
			// 如果该树节点是叶子节点，则计算多极矩，并使用M2M将其累加到其父节点中；
			if ( qtree.treeNodeList[cellIndex] -> isLeaf == 1)  // 是叶子节点
			{
				// 计算多极矩，（调用一个计算接口）， 
				// 输出是该节点的多极矩，输入是该节点的单元编号，
				// 计算编号为cellindex的树节点的多极矩
				// 计算多极矩的时候需要将迭代向量代入计算

				cal_multipole_moment(qtree,cellIndex,mesh,data); 
			}

			//if ( qtree.treeNodeList[cellIndex].isLeaf == 0)
			// 将多极矩传递到父节点中去， {不论当前节点是否是叶子节点}   
			{
				// 传递多极系数使用（M2M公式将当前节点的
				// 多极矩进行转换并累加到父节点中），输入是当前节点的多极矩，

				transfer_mm_to_its_father_by_m2m(qtree,cellIndex);
			}
		} // 遍历层结束

	} // 遍历树结束

	return 0;
}

// 下行遍历树结构，计算局部矩，这里的局部矩包括交互单元和远程单元
// 对整个积分的贡献。

int quadtree_downward(Quadtree& qtree)
{
	// 计算第二层的树节点的多极矩, 使用M2L将交互节点的多极矩转换为局部矩，
	// 累加到中心节点中, 初始计算参数
	{
		// 这一步可以统一到整个树的计算过程中。
	}	
	// 逐层向根方向遍历树结构
	int maxLevel = qtree.treeDepth - 1;  // retrive max tree level.
	for (int depth = 2; depth <= maxLevel; ++ depth)
	{
		int startIndex = qtree.levelCell[ depth ];
		int endIndex = startIndex + qtree.numlevelCell[ depth ] - 1;

		// 遍历同层次的树节点
		for (int cellIndex = startIndex; cellIndex <= endIndex; ++cellIndex) 
		{
			QuadtreeNode& qnode = *qtree.treeNodeList[cellIndex];
			// 遍历交互节点列表，使用M2L，将局部矩累加到节点上
			// 需要查询交互节点
			for (int i = 0; i < qnode.lenInterList; i++)
			{
				int iIndex = qnode.interList[i];
				transfer_lm_by_m2l(qtree, cellIndex, iIndex);	
			}
			// 将远程单元的贡献累加到当前单元中，
			// 只要使用L2L将父节点的局部系数转移到当前节点的中心上
			{
				if (depth == 2 ) continue; // 如果当前节点是第二层节点，则不必使用L2L
				transfer_lm_by_l2l(qtree, cellIndex);
			}
		} // 遍历层

	} // 遍历树
}

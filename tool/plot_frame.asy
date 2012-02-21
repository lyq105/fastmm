//
// This file is use to generate the animate of quad tree which is used 
// for fast multipole method.
//
// Author: Li Yiqiang (lyq105@163.com)

void plot_frame(string filename)
{
  file fin = input(filename);

  int number;
  real a,b;
  pair [] c;
  size(200);

  number = fin;
  for (int j = 0; j<number; ++j)
  {
    for (int i = 0; i<4; ++i)
    {
      a = fin;
      b = fin;
      c[i] = (a,b); 
      if(eof(fin))break;
    }
    if (j == 0)
      filldraw (c[0]--c[1]--c[2]--c[3]--cycle,fillpen=yellow,drawpen=red);
    else
      filldraw (c[0]--c[1]--c[2]--c[3]--cycle,fillpen=palered,drawpen=red);
  }

  number = fin;
  write(stdout,number);
  for (int j = 0; j<number; ++j)
  {
    for (int i = 0; i<4; ++i)
    {
      a = fin;
      b = fin;
      c[i] = (a,b); 
      if(eof(fin))break;
    }
    filldraw (c[0]--c[1]--c[2]--c[3]--cycle,fillpen=paleblue,drawpen=blue);
  }

  number = fin;
  for (int j = 0; j<number; ++j)
  {
    for (int i = 0; i<4; ++i)
    {
      a = fin;
      b = fin;
      c[i] = (a,b); 
      if(eof(fin))break;
    }
    draw (c[0]--c[1]--c[2]--c[3]--cycle);
  }
}


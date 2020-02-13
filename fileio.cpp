#include "headers.h"
#include <sstream>

//------------------------------------------------------------------------------------------  
std :: pair<int, int> dim(const string& name){

  std :: ifstream fin(name);
  if(!fin) ferror("can't open file\n");

  string s;
  std :: getline (fin, s);
  std :: istringstream ss{s};
  int col = 0, row = 1;
  
  while (std :: getline(ss, s, '>')) ++col;
  while (std :: getline(fin, s)) ++row;

  auto d = std :: make_pair(row, col);
  fin.close();
  
  return d;
}
//------------------------------------------------------------------------------------------  
void read_matrix(const string& name, Matrix2s& m)
{
  std :: ifstream fin(name);
  if(!fin) ferror("read_file can't open file\n");

  string s;
  for (int i=0; std :: getline(fin, s); ++i)
    {
      std :: istringstream ss{s};
      for(int j=0; std :: getline(ss, s, '>'); ++j)
	m(i,j) = s;
    }
  fin.close();

}
//------------------------------------------------------------------------------------------  
void init_matrix(const string& file, Matrix2s& m)
{
  auto d = dim(file);
  m = Matrix2s(d.first, d.second);
  read_matrix(file, m);
}
//------------------------------------------------------------------------------------------

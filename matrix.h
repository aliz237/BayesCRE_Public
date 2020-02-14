#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <iostream>
#include <string>
#include <stdio.h>

#define HERE fprintf(stderr, "HERE at %s:%d (%s)\n", __FILE__, __LINE__, __FUNCTION__)

struct error{
  std::string mess;
error(const std::string& s): mess(s) {}
};

template <class T> class matrix{
  T *vals;
  int rows;
  int cols;
public:
  T* getvals() { return vals;}
  int ncols() const { return cols;}
  int nrows() const { return rows;}
  int size() const { return rows * cols;}
  matrix(){ vals = 0; rows = cols = 0;}
  matrix(int r, int c){
    vals = new T[c*r];
    rows = r;
    cols = c;
  }
  matrix(const matrix& m){
    vals = new T[m.size()];
    std::copy(m.vals, m.vals+m.size(), vals);
    rows = m.rows;
    cols = m.cols;
  }
  matrix & operator = (const matrix& m){
    T *temp = new T[m.size()];
    std::copy(m.vals, m.vals+m.size(), temp);
    std::swap(temp, vals);
    rows = m.rows;
    cols = m.cols;
    return *this;
  }
    
  ~matrix(){
    delete[] vals;
  }
  
  void row_check(int r) const{
    if( r < 0 || r >= rows)
      throw(error("out of range row: " + std::to_string(r) + '\n'));
  }

  void col_check(int c) const{
    if( c < 0 || c >= cols)
      throw(error("out of range column: " + std::to_string(c) + '\n'));
  }
  
  T& operator () (int r, int c){
    row_check(r);
    col_check(c);
    return vals[r * cols + c];
  }

  T operator () (int r, int c) const{
    row_check(r);
    col_check(c);
    return vals[r * cols + c];
  }

  matrix<T>& operator *= (const T x){
    for (int i=0; i<size(); i++)
      vals[i] *= x;
    
    return *this;
  }
  
};  

template <class T>
matrix<T> transpose (const matrix<T>& m){
  matrix<T> tr(m.ncols(), m.nrows());
  
  for (int i=0; i<m.nrows(); i++){
    for (int j=0; j<m.ncols(); j++){
      tr(j,i) = m(i,j);
    }
  }
  // write a move constructor 
  return tr;
}
  
template <class T>
std::ostream& operator << (std::ostream& o, const matrix<T>& m)
{
  for (size_t i=0; i<m.nrows(); ++i)
    {
      for (size_t j=0; j<m.ncols(); ++j)
	o << m(i,j) << " ";

      o << "\n";
    }
  return o;
}

template <class T> class matrix_column{
  matrix<T>& data;
  int col;
 public:

 matrix_column(matrix<T>& m, int n):
  data(m)
    {
      m.col_check(n);
      col = n;
    }

  int size() const{ return data.nrows();}
  
  T& operator [] (int r){
    data.row_check(r);
    return data(r, col);
  }

  T operator [] (int r) const{
    data.row_check(r);
    return data(r, col);
  }

  //for range for
  T* begin() { return data.getvals();}
  //const T* begin() { return data.getval();}
  T* end() { return data.getvals() + data.size();}
  //const T* end() { return data.getvals() + data.size();}
  
};

template <class T>
std::ostream& operator << (std::ostream& o, const matrix_column<T>& m)
{
  for (size_t r=0; r<m.size(); r++)
    o << m[r] << " ";

  return o;
}

#endif

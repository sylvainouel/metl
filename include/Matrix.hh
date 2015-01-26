#ifndef MATRIX_HH
#define MATRIX_HH

/*
metl: A generic framework for sequential and parallel metaheuristics
Copyright (c) 2005-2015, Sylvain Ouellet


Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

*/

#include <assert.h>

namespace metl {


template <class T> class Matrix;


template<class T>
struct Matrix_iterator {
  Matrix_iterator(Matrix<T>& point_to)
    : m(point_to),
      i(0),
      j(0)
  {}

  Matrix_iterator(Matrix<T>& point_to, unsigned init_i, unsigned init_j)
    : m(point_to),
      i(init_i),
      j(init_j)
  {}


  inline T& operator*()
  {
    return m(i,j);
  }

  inline const T& operator*() const
  {
    return m(i,j);
  }

  inline Matrix_iterator& operator++()
  {
    if (++j >= m.cols) {
      ++i;
      j=0;
    }
    return *this;
  }

  inline bool operator==(const Matrix_iterator& rhs) const {
    return (i==rhs.i && j==rhs.j && &m==&rhs.m);
  }

  inline bool operator!=(const Matrix_iterator& rhs) const {
    return !(operator==(rhs));
  }


  inline unsigned get_i() const { return i; }
  inline unsigned get_j() const { return j; }

private:
  Matrix<T>& m;
  unsigned i, j;

  friend class Matrix<T>;
};



template<class T>
class Matrix {
public:
  typedef Matrix_iterator<T> iterator;

  Matrix();
  Matrix(unsigned nrows, unsigned ncols);
  ~Matrix();

  Matrix(const Matrix& m);
  Matrix& operator= (const Matrix& m);


  // Access methods to get the (i,j) element:
  inline T&       operator() (unsigned i, unsigned j) {  
    return data_[i][j];
  }
  inline const T& operator() (unsigned i, unsigned j) const {
    return data_[i][j];
  }

  unsigned get_rows() const { return rows; }
  unsigned get_cols() const { return cols; }

  iterator begin() { return iterator(*this); }
  const iterator end() {
    return iterator(*this, rows, 0);
  }

private:
  T** data_;
  unsigned rows, cols;

  friend class Matrix_iterator<T>;
};



template<class T>
Matrix<T>::Matrix()
  :data_(0), rows(0), cols(0)
{}


template<class T>
Matrix<T>::Matrix(const Matrix<T>& m)
  : data_(0),
    rows(m.rows),
    cols(m.cols)
{
  if (rows==0) return;

      // allocate memory
  data_ = new T*[rows];
  for (unsigned i=0; i<rows; ++i)
    data_[i] = new T[cols];

  for (unsigned i=0; i<rows; ++i)
    for (unsigned j=0; j<cols; ++j)
      data_[i][j] = m.data_[i][j];
}


template <class T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& m)
{
  if (this == &m) return *this;

  if (cols!=m.cols || rows!=m.rows) {
    for (unsigned i=0; i<rows; ++i)
      delete[] data_[i];

    if (rows!=m.rows) {
      if (rows!=0) delete[] data_;
      if (m.rows>0)
	data_ = new T*[m.rows];
      else
	data_ = 0;
    }

    rows = m.rows;
    cols = m.cols;

    for (unsigned i=0; i<m.rows; ++i) {
      data_[i] = new T[m.cols];
    }
  }

  // copy data from old matrix to new one
  for (unsigned i=0; i<rows; ++i)
    for (unsigned j=0; j<cols; ++j)
      data_[i][j] = m.data_[i][j];

  return *this;
}


template<class T>
Matrix<T>::Matrix(unsigned nrows, unsigned ncols)
  : data_(0),
    rows(nrows),
    cols(ncols)
{
  assert(rows!=0 || cols==0);

  if (rows==0) return;

  data_ = new T*[nrows];
  for (unsigned i=0; i<nrows; ++i)
    data_[i] = new T[ncols];

  for (unsigned i=0; i<nrows; ++i)
    for (unsigned j=0; j<ncols; ++j)
      data_[i][j] = 0;
}

template<class T>
Matrix<T>::~Matrix()
{
  if (rows>0) {
    for (unsigned i=0; i<rows; ++i)
      delete[] data_[i];

    delete[] data_;
  }
}

}

#endif

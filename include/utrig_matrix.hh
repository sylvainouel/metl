#ifndef UTRIG_MATRIX_HH
#define UTRIG_MATRIX_HH
/*
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

// that is an upper triangular matrx. i.e. only elements m(i,j) where
// j>i are accessibles


namespace metl {


template <class T> class utrig_matrix;


template<class T>
struct utrig_matrix_iterator {
  utrig_matrix_iterator(utrig_matrix<T>& point_to)
    : m(point_to),
      i(0),
      j(1)
  {}

  utrig_matrix_iterator(utrig_matrix<T>& point_to, unsigned init_i, unsigned init_j)
    : m(point_to),
      i(init_i),
      j(init_j)
  {
    assert(j>i);
  }


  inline T& operator*()
  {
    return m(i,j);
  }

  inline const T& operator*() const
  {
    return m(i,j);
  }

  inline utrig_matrix_iterator& operator++()
  {
    assert(j>i);
    assert(j<m.cols);

    if (++j == m.cols) {
      ++i;
      j=i+1;
    }
    return *this;
  }

  inline bool operator==(const utrig_matrix_iterator& rhs) const {
    return (&m==&rhs.m && i==rhs.i && j==rhs.j);
  }

  inline bool operator!=(const utrig_matrix_iterator& rhs) const {
    return !(operator==(rhs));
  }


  inline unsigned get_i() const { return i; }  // this is needed. we need to know were the iterator points.
  inline unsigned get_j() const { return j; }

private:
  utrig_matrix<T>& m;
  unsigned i, j;

  friend class utrig_matrix<T>;
};


template<class T>
class utrig_matrix {
public:
  typedef utrig_matrix_iterator<T> iterator;

  utrig_matrix() {};
  utrig_matrix(unsigned nrows, unsigned ncols);
  ~utrig_matrix();

  utrig_matrix(const utrig_matrix& m);
  utrig_matrix& operator= (const utrig_matrix& m);


  // Access methods to get the (i,j) element:
  inline T&       operator() (unsigned i, unsigned j) { 
    assert(j>i);
    return data_[i][j-i-1];
  };
  inline const T& operator() (unsigned i, unsigned j) const { 
    assert(j>i);
    return data_[i][j-i-1];
  };

  iterator begin() { return iterator(*this); }
  const iterator end() { return iterator(*this, rows-1, cols);  }

private:
  T** data_;
  unsigned rows, cols;
  friend class utrig_matrix_iterator<T>;
};


template<class T>
utrig_matrix<T>::utrig_matrix(const utrig_matrix<T>& m)
  : data_(0),
    rows(m.rows),
    cols(m.cols)
{
  if (rows==0) return;

  // allocate memory
  data_ = new T*[rows];
  for (unsigned i=0; i<rows; ++i)
    data_[i] = new T[cols-i-1];

  for (unsigned i=0; i<rows; ++i)
    for (unsigned j=i+1; j<cols; ++j)
      operator()(i,j) = m.operator()(i,j);
}


template <class T>
utrig_matrix<T>& utrig_matrix<T>::operator= (const utrig_matrix<T>& m)
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

    rows=m.rows;
    cols=m.cols;

    for (unsigned i=0; i<m.rows; ++i) {
      data_[i] = new T[m.cols];
    }
  }

  // copy data from old matrix to new one
  for (unsigned i=0; i<rows; ++i)
    for (unsigned j=i+1; j<cols; ++j)
      operator()(i,j) = m.operator()(i,j);

  return *this;
}


template<class T>
utrig_matrix<T>::utrig_matrix(unsigned nrows, unsigned ncols)
  : data_(0),
    rows(nrows),
    cols(ncols)
{
  if (rows==0) return;

  data_ = new T*[nrows];
  for (unsigned i=0; i<nrows; ++i)
    data_[i] = new T[ncols];

  for (unsigned i=0; i<nrows; ++i)
    for (unsigned j=i+1; j<ncols; ++j)
      operator()(i,j) = 0;
}

template<class T>
utrig_matrix<T>::~utrig_matrix()
{
  if (rows>0) {
    for (unsigned i=0; i<rows; ++i)
      delete[] data_[i];

    delete[] data_;
  }
}

}

#endif

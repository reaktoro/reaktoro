// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_sort
//! @{



template<typename eT>
class arma_ascend_sort_helper
  {
  public:
  
  arma_inline
  bool
  operator() (eT a, eT b) const
    {
    return (a < b);
    }
  };
  


template<typename eT>
class arma_descend_sort_helper
  {
  public:
  
  arma_inline
  bool
  operator() (eT a, eT b) const
    {
    return (a > b);
    }
  };
  


template<typename T>
class arma_ascend_sort_helper< std::complex<T> >
  {
  public:
  
  typedef typename std::complex<T> eT;
  
  inline
  bool
  operator() (const eT& a, const eT& b) const
    {
    return (std::abs(a) < std::abs(b));
    }
  };



template<typename T>
class arma_descend_sort_helper< std::complex<T> >
  {
  public:
  
  typedef typename std::complex<T> eT;
  
  inline
  bool
  operator() (const eT& a, const eT& b) const
    {
    return (std::abs(a) > std::abs(b));
    }
  };



template<typename eT>
inline 
void
op_sort::direct_sort(eT* X, const uword n_elem, const uword sort_type)
  {
  arma_extra_debug_sigprint();
  
  if(sort_type == 0)
    {
    arma_ascend_sort_helper<eT> comparator;
    
    std::sort(&X[0], &X[n_elem], comparator);
    }
  else
    {
    arma_descend_sort_helper<eT> comparator;
    
    std::sort(&X[0], &X[n_elem], comparator);
    }
  }



template<typename eT>
inline 
void
op_sort::direct_sort_ascending(eT* X, const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  arma_ascend_sort_helper<eT> comparator;
    
  std::sort(&X[0], &X[n_elem], comparator);
  }



template<typename eT>
inline 
void
op_sort::copy_row(eT* X, const Mat<eT>& A, const uword row)
  {
  const uword N = A.n_cols;
  
  uword i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    X[i] = A.at(row,i);
    X[j] = A.at(row,j);
    }
  
  if(i < N)
    {
    X[i] = A.at(row,i);
    }
  }



template<typename eT>
inline 
void
op_sort::copy_row(Mat<eT>& A, const eT* X, const uword row)
  {
  const uword N = A.n_cols;
  
  uword i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    A.at(row,i) = X[i]; 
    A.at(row,j) = X[j];
    }
  
  if(i < N)
    {
    A.at(row,i) = X[i];
    }
  }



template<typename T1>
inline
void
op_sort::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sort>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1>   tmp(in.m, out);
  const Mat<eT>&       X = tmp.M;
  
  const uword sort_type = in.aux_uword_a;
  const uword dim       = in.aux_uword_b;
  
  arma_debug_check( (sort_type > 1),          "sort(): incorrect usage. sort_type must be 0 or 1");
  arma_debug_check( (dim > 1),                "sort(): incorrect usage. dim must be 0 or 1"      );
  arma_debug_check( (X.is_finite() == false), "sort(): given object has non-finite elements"     );
  
  if( (X.n_rows * X.n_cols) <= 1 )
    {
    out = X;
    return;
    }
  
  
  if(dim == 0)  // sort the contents of each column
    {
    arma_extra_debug_print("op_sort::apply(), dim = 0");
    
    out = X;
    
    const uword n_rows = out.n_rows;
    const uword n_cols = out.n_cols;
      
    for(uword col=0; col < n_cols; ++col)
      {
      op_sort::direct_sort( out.colptr(col), n_rows, sort_type );
      }
    }
  else
  if(dim == 1)  // sort the contents of each row
    {
    if(X.n_rows == 1)  // a row vector
      {
      arma_extra_debug_print("op_sort::apply(), dim = 1, vector specific");
      
      out = X;
      op_sort::direct_sort(out.memptr(), out.n_elem, sort_type);
      }
    else  // not a row vector
      {
      arma_extra_debug_print("op_sort::apply(), dim = 1, generic");
      
      out.copy_size(X);
      
      const uword n_rows = out.n_rows;
      const uword n_cols = out.n_cols;
      
      podarray<eT> tmp_array(n_cols);
      
      for(uword row=0; row < n_rows; ++row)
        {
        op_sort::copy_row(tmp_array.memptr(), X, row);
        
        op_sort::direct_sort( tmp_array.memptr(), n_cols, sort_type );
        
        op_sort::copy_row(out, tmp_array.memptr(), row);
        }
      }
    }
  
  }


//! @}

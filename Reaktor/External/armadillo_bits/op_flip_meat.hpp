// Copyright (C) 2010 Conrad Sanderson
// Copyright (C) 2010 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_flip
//! @{



template<typename T1>
inline
void
op_flipud::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_flipud>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>  tmp(in.m);
  const Mat<eT> X = tmp.M;
  
  if(&out != &X)
    {
    out.copy_size(X);
    
    for(uword i=0; i<X.n_rows; ++i)
      {
      out.row(i) = X.row(X.n_rows-1 - i);
      }
    }
  else
    {
    const uword N = X.n_rows / 2;
    
    for(uword i=0; i<N; ++i)
      {
      out.swap_rows(i, X.n_rows-1 - i);
      }
    }
  }



template<typename T1>
inline
void
op_fliplr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_fliplr>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>  tmp(in.m);
  const Mat<eT> X = tmp.M;
  
  if(&out != &X)
    {
    out.copy_size(X);
    
    for(uword i=0; i<X.n_cols; ++i)
      {
      out.col(i) = X.col(X.n_cols-1 - i);
      }
    }
  else
    {
    const uword N = X.n_cols / 2;
    
    for(uword i=0; i<N; ++i)
      {
      out.swap_cols(i, X.n_cols-1 - i);
      }
    }
  }



//! @}

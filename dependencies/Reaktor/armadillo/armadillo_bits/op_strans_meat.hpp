// Copyright (C) 2008-2013 Conrad Sanderson
// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Ryan Curtin
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_strans
//! @{



//! for tiny square matrices (size <= 4x4)
template<typename eT, typename TA>
arma_hot
inline
void
op_strans::apply_mat_noalias_tinysq(Mat<eT>& out, const TA& A)
  {
  const eT*   Am =   A.memptr();
        eT* outm = out.memptr();
  
  switch(A.n_rows)
    {
    case 1:
      {
      outm[0] = Am[0];
      }
      break;
      
    case 2:
      {
      outm[pos<false,0,0>::n2] = Am[pos<true,0,0>::n2];
      outm[pos<false,1,0>::n2] = Am[pos<true,1,0>::n2];
      
      outm[pos<false,0,1>::n2] = Am[pos<true,0,1>::n2];
      outm[pos<false,1,1>::n2] = Am[pos<true,1,1>::n2];
      }
      break;
    
    case 3:
      {
      outm[pos<false,0,0>::n3] = Am[pos<true,0,0>::n3];
      outm[pos<false,1,0>::n3] = Am[pos<true,1,0>::n3];
      outm[pos<false,2,0>::n3] = Am[pos<true,2,0>::n3];
      
      outm[pos<false,0,1>::n3] = Am[pos<true,0,1>::n3];
      outm[pos<false,1,1>::n3] = Am[pos<true,1,1>::n3];
      outm[pos<false,2,1>::n3] = Am[pos<true,2,1>::n3];
      
      outm[pos<false,0,2>::n3] = Am[pos<true,0,2>::n3];
      outm[pos<false,1,2>::n3] = Am[pos<true,1,2>::n3];
      outm[pos<false,2,2>::n3] = Am[pos<true,2,2>::n3];
      }
      break;
    
    case 4:
      {
      outm[pos<false,0,0>::n4] = Am[pos<true,0,0>::n4];
      outm[pos<false,1,0>::n4] = Am[pos<true,1,0>::n4];
      outm[pos<false,2,0>::n4] = Am[pos<true,2,0>::n4];
      outm[pos<false,3,0>::n4] = Am[pos<true,3,0>::n4];
      
      outm[pos<false,0,1>::n4] = Am[pos<true,0,1>::n4];
      outm[pos<false,1,1>::n4] = Am[pos<true,1,1>::n4];
      outm[pos<false,2,1>::n4] = Am[pos<true,2,1>::n4];
      outm[pos<false,3,1>::n4] = Am[pos<true,3,1>::n4];
      
      outm[pos<false,0,2>::n4] = Am[pos<true,0,2>::n4];
      outm[pos<false,1,2>::n4] = Am[pos<true,1,2>::n4];
      outm[pos<false,2,2>::n4] = Am[pos<true,2,2>::n4];
      outm[pos<false,3,2>::n4] = Am[pos<true,3,2>::n4];
      
      outm[pos<false,0,3>::n4] = Am[pos<true,0,3>::n4];
      outm[pos<false,1,3>::n4] = Am[pos<true,1,3>::n4];
      outm[pos<false,2,3>::n4] = Am[pos<true,2,3>::n4];
      outm[pos<false,3,3>::n4] = Am[pos<true,3,3>::n4];
      }
      break;
    
    default:
      ;
    }
  
  }



//! Immediate transpose of a dense matrix
template<typename eT, typename TA>
arma_hot
inline
void
op_strans::apply_mat_noalias(Mat<eT>& out, const TA& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_cols = A.n_cols;
  const uword A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (TA::is_row) || (TA::is_col) || (A_n_cols == 1) || (A_n_rows == 1) )
    {
    arrayops::copy( out.memptr(), A.memptr(), A.n_elem );
    }
  else
    {
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      op_strans::apply_mat_noalias_tinysq(out, A);
      }
    else
      {
      for(uword k=0; k < A_n_cols; ++k)
        {
        uword i, j;
        
        const eT* colptr = A.colptr(k);
        
        for(i=0, j=1; j < A_n_rows; i+=2, j+=2)
          {
          const eT tmp_i = colptr[i];
          const eT tmp_j = colptr[j];
          
          out.at(k, i) = tmp_i;
          out.at(k, j) = tmp_j;
          }
        
        if(i < A_n_rows)
          {
          out.at(k, i) = colptr[i];
          }
        }
      }
    }
  }



template<typename eT>
arma_hot
inline
void
op_strans::apply_mat_inplace(Mat<eT>& out)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = out.n_rows;
  const uword n_cols = out.n_cols;
  
  if(n_rows == n_cols)
    {
    arma_extra_debug_print("op_strans::apply(): doing in-place transpose of a square matrix");
    
    const uword N = n_rows;
    
    for(uword k=0; k < N; ++k)
      {
      eT* colptr = out.colptr(k);
      
      uword i,j;
      
      for(i=(k+1), j=(k+2); j < N; i+=2, j+=2)
        {
        std::swap(out.at(k,i), colptr[i]);
        std::swap(out.at(k,j), colptr[j]);
        }
      
      if(i < N)
        {
        std::swap(out.at(k,i), colptr[i]);
        }
      }
    }
  else
    {
    Mat<eT> tmp;
    
    op_strans::apply_mat_noalias(tmp, out);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename TA>
arma_hot
inline
void
op_strans::apply_mat(Mat<eT>& out, const TA& A)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &A)
    {
    op_strans::apply_mat_noalias(out, A);
    }
  else
    {
    op_strans::apply_mat_inplace(out);
    }
  }



template<typename T1>
arma_hot
inline
void
op_strans::apply_proxy(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X);
  
  // allow detection of in-place transpose
  if( (is_Mat<typename Proxy<T1>::stored_type>::value == true) && (Proxy<T1>::fake_mat == false) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    op_strans::apply_mat(out, tmp.M);
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    const bool is_alias = P.is_alias(out);
    
    if( (resolves_to_vector<T1>::value == true) && (Proxy<T1>::prefer_at_accessor == false) )
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        eT* out_mem = out.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        uword i,j;
        for(i=0, j=1; j < n_elem; i+=2, j+=2)
          {
          const eT tmp_i = Pea[i];
          const eT tmp_j = Pea[j];
          
          out_mem[i] = tmp_i;
          out_mem[j] = tmp_j;
          }
        
        if(i < n_elem)
          {
          out_mem[i] = Pea[i];
          }
        }
      else  // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        eT* out_mem = out2.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        uword i,j;
        for(i=0, j=1; j < n_elem; i+=2, j+=2)
          {
          const eT tmp_i = Pea[i];
          const eT tmp_j = Pea[j];
          
          out_mem[i] = tmp_i;
          out_mem[j] = tmp_j;
          }
        
        if(i < n_elem)
          {
          out_mem[i] = Pea[i];
          }
        
        out.steal_mem(out2);
        }
      }
    else   // general matrix transpose
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
          {
          uword i, j;
          
          for(i=0, j=1; j < n_rows; i+=2, j+=2)
            {
            const eT tmp_i = P.at(i,k);
            const eT tmp_j = P.at(j,k);
            
            out.at(k,i) = tmp_i;
            out.at(k,j) = tmp_j;
            }
          
          if(i < n_rows)
            {
            out.at(k,i) = P.at(i,k);
            }
          }
        }
      else // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
          {
          uword i, j;
          
          for(i=0, j=1; j < n_rows; i+=2, j+=2)
            {
            const eT tmp_i = P.at(i,k);
            const eT tmp_j = P.at(j,k);
            
            out2.at(k,i) = tmp_i;
            out2.at(k,j) = tmp_j;
            }
          
          if(i < n_rows)
            {
            out2.at(k,i) = P.at(i,k);
            }
          }
        
        out.steal_mem(out2);
        }
      }
    }
  }



template<typename T1>
arma_hot
inline
void
op_strans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in)
  {
  arma_extra_debug_sigprint();
  
  op_strans::apply_proxy(out, in.m);
  }



//
// op_strans2



//! for tiny square matrices (size <= 4x4)
template<typename eT, typename TA>
arma_hot
inline
void
op_strans2::apply_noalias_tinysq(Mat<eT>& out, const TA& A, const eT val)
  {
  const eT* Am   =   A.memptr();
        eT* outm = out.memptr();
  
  switch(A.n_rows)
    {
    case 1:
      {
      outm[0] = val * Am[0];
      }
      break;
      
    case 2:
      {
      outm[pos<false,0,0>::n2] = val * Am[pos<true,0,0>::n2];
      outm[pos<false,1,0>::n2] = val * Am[pos<true,1,0>::n2];
      
      outm[pos<false,0,1>::n2] = val * Am[pos<true,0,1>::n2];
      outm[pos<false,1,1>::n2] = val * Am[pos<true,1,1>::n2];
      }
      break;
    
    case 3:
      {
      outm[pos<false,0,0>::n3] = val * Am[pos<true,0,0>::n3];
      outm[pos<false,1,0>::n3] = val * Am[pos<true,1,0>::n3];
      outm[pos<false,2,0>::n3] = val * Am[pos<true,2,0>::n3];
      
      outm[pos<false,0,1>::n3] = val * Am[pos<true,0,1>::n3];
      outm[pos<false,1,1>::n3] = val * Am[pos<true,1,1>::n3];
      outm[pos<false,2,1>::n3] = val * Am[pos<true,2,1>::n3];
      
      outm[pos<false,0,2>::n3] = val * Am[pos<true,0,2>::n3];
      outm[pos<false,1,2>::n3] = val * Am[pos<true,1,2>::n3];
      outm[pos<false,2,2>::n3] = val * Am[pos<true,2,2>::n3];
      }
      break;
    
    case 4:
      {
      outm[pos<false,0,0>::n4] = val * Am[pos<true,0,0>::n4];
      outm[pos<false,1,0>::n4] = val * Am[pos<true,1,0>::n4];
      outm[pos<false,2,0>::n4] = val * Am[pos<true,2,0>::n4];
      outm[pos<false,3,0>::n4] = val * Am[pos<true,3,0>::n4];
      
      outm[pos<false,0,1>::n4] = val * Am[pos<true,0,1>::n4];
      outm[pos<false,1,1>::n4] = val * Am[pos<true,1,1>::n4];
      outm[pos<false,2,1>::n4] = val * Am[pos<true,2,1>::n4];
      outm[pos<false,3,1>::n4] = val * Am[pos<true,3,1>::n4];
      
      outm[pos<false,0,2>::n4] = val * Am[pos<true,0,2>::n4];
      outm[pos<false,1,2>::n4] = val * Am[pos<true,1,2>::n4];
      outm[pos<false,2,2>::n4] = val * Am[pos<true,2,2>::n4];
      outm[pos<false,3,2>::n4] = val * Am[pos<true,3,2>::n4];
      
      outm[pos<false,0,3>::n4] = val * Am[pos<true,0,3>::n4];
      outm[pos<false,1,3>::n4] = val * Am[pos<true,1,3>::n4];
      outm[pos<false,2,3>::n4] = val * Am[pos<true,2,3>::n4];
      outm[pos<false,3,3>::n4] = val * Am[pos<true,3,3>::n4];
      }
      break;
    
    default:
      ;
    }
  
  }



template<typename eT, typename TA>
arma_hot
inline
void
op_strans2::apply_noalias(Mat<eT>& out, const TA& A, const eT val)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_cols = A.n_cols;
  const uword A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (TA::is_col) || (TA::is_row) || (A_n_cols == 1) || (A_n_rows == 1) )
    {
    const uword N = A.n_elem;
    
    const eT*   A_mem =   A.memptr();
          eT* out_mem = out.memptr();
    
    uword i,j;
    for(i=0, j=1; j < N; i+=2, j+=2)
      {
      const eT tmp_i = A_mem[i];
      const eT tmp_j = A_mem[j];
      
      out_mem[i] = val * tmp_i;
      out_mem[j] = val * tmp_j;
      }
    
    if(i < N)
      {
      out_mem[i] = val * A_mem[i];
      }
    }
  else
    {
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      op_strans2::apply_noalias_tinysq(out, A, val);
      }
    else
      {
      for(uword k=0; k < A_n_cols; ++k)
        {
        uword i, j;
        
        const eT* colptr = A.colptr(k);
        
        for(i=0, j=1; j < A_n_rows; i+=2, j+=2)
          {
          const eT tmp_i = colptr[i];
          const eT tmp_j = colptr[j];
          
          out.at(k, i) = val * tmp_i;
          out.at(k, j) = val * tmp_j;
          }
        
        if(i < A_n_rows)
          {
          out.at(k, i) = val * colptr[i];
          }
        }
      }
    }
  }



template<typename eT, typename TA>
arma_hot
inline
void
op_strans2::apply(Mat<eT>& out, const TA& A, const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &A)
    {
    op_strans2::apply_noalias(out, A, val);
    }
  else
    {
    const uword n_rows = out.n_rows;
    const uword n_cols = out.n_cols;
    
    if(n_rows == n_cols)
      {
      arma_extra_debug_print("op_strans2::apply(): doing in-place transpose of a square matrix");
      
      const uword N = n_rows;
      
      // TODO: do multiplication while swapping
      
      for(uword k=0; k < N; ++k)
        {
        eT* colptr = out.colptr(k);
        
        uword i,j;
        
        for(i=(k+1), j=(k+2); j < N; i+=2, j+=2)
          {
          std::swap(out.at(k,i), colptr[i]);
          std::swap(out.at(k,j), colptr[j]);
          }
        
        if(i < N)
          {
          std::swap(out.at(k,i), colptr[i]);
          }
        }
      
      arrayops::inplace_mul( out.memptr(), val, out.n_elem );
      }
    else
      {
      Mat<eT> tmp;
      op_strans2::apply_noalias(tmp, A, val);
      
      out.steal_mem(tmp);
      }
    }
  }



template<typename T1>
arma_hot
inline
void
op_strans2::apply_proxy(Mat<typename T1::elem_type>& out, const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X);
  
  // allow detection of in-place transpose
  if( (is_Mat<typename Proxy<T1>::stored_type>::value == true) && (Proxy<T1>::fake_mat == false) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    op_strans2::apply(out, tmp.M, val);
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    const bool is_alias = P.is_alias(out);
    
    if( (resolves_to_vector<T1>::value == true) && (Proxy<T1>::prefer_at_accessor == false) )
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        eT* out_mem = out.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        uword i,j;
        for(i=0, j=1; j < n_elem; i+=2, j+=2)
          {
          const eT tmp_i = Pea[i];
          const eT tmp_j = Pea[j];
          
          out_mem[i] = val * tmp_i;
          out_mem[j] = val * tmp_j;
          }
        
        if(i < n_elem)
          {
          out_mem[i] = val * Pea[i];
          }
        }
      else  // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        eT* out_mem = out2.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        uword i,j;
        for(i=0, j=1; j < n_elem; i+=2, j+=2)
          {
          const eT tmp_i = Pea[i];
          const eT tmp_j = Pea[j];
          
          out_mem[i] = val * tmp_i;
          out_mem[j] = val * tmp_j;
          }
        
        if(i < n_elem)
          {
          out_mem[i] = val * Pea[i];
          }
        
        out.steal_mem(out2);
        }
      }
    else   // general matrix transpose
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
          {
          uword i, j;
          
          for(i=0, j=1; j < n_rows; i+=2, j+=2)
            {
            const eT tmp_i = P.at(i,k);
            const eT tmp_j = P.at(j,k);
            
            out.at(k,i) = val * tmp_i;
            out.at(k,j) = val * tmp_j;
            }
          
          if(i < n_rows)
            {
            out.at(k,i) = val * P.at(i,k);
            }
          }
        }
      else // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
          {
          uword i, j;
          
          for(i=0, j=1; j < n_rows; i+=2, j+=2)
            {
            const eT tmp_i = P.at(i,k);
            const eT tmp_j = P.at(j,k);
            
            out2.at(k,i) = val * tmp_i;
            out2.at(k,j) = val * tmp_j;
            }
          
          if(i < n_rows)
            {
            out2.at(k,i) = val * P.at(i,k);
            }
          }
        
        out.steal_mem(out2);
        }
      }
    }
  }



//! @}

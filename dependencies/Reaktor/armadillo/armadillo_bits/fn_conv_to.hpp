// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_conv_to
//! @{



//! conversion from Armadillo Base and BaseCube objects to scalars
//! (kept only for compatibility with old code; use as_scalar() instead for Base objects like Mat)
template<typename out_eT>
class conv_to
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static out_eT from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);

  template<typename in_eT, typename T1>
  inline static out_eT from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static out_eT from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static out_eT from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem != 1), "conv_to(): given object doesn't have exactly one element" );
  
  return out_eT(X.mem[0]);
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem != 1), "conv_to(): given object doesn't have exactly one element" );
  
  out_eT out;
  
  arrayops::convert_cx_scalar(out, X.mem[0]);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const unwrap_cube<T1>  tmp(in.get_ref());
  const Cube<in_eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem != 1), "conv_to(): given object doesn't have exactly one element" );
  
  return out_eT(X.mem[0]);
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const unwrap_cube<T1>  tmp(in.get_ref());
  const Cube<in_eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem != 1), "conv_to(): given object doesn't have exactly one element" );
  
  out_eT out;
  
  arrayops::convert_cx_scalar(out, X.mem[0]);
  
  return out;
  }



//! conversion to Armadillo matrices from Armadillo Base objects, as well as from std::vector
template<typename out_eT>
class conv_to< Mat<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  
  
  
  template<typename in_eT>
  inline static Mat<out_eT> from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT>
  inline static Mat<out_eT> from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  Mat<out_eT> out(X.n_rows, X.n_cols);
  
  arrayops::convert( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  Mat<out_eT> out(X.n_rows, X.n_cols);
  
  arrayops::convert_cx( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Mat<out_eT> out(N, 1);
  
  if(N > 0)
    {
          out_eT* out_mem = out.memptr();
    const  in_eT*  in_mem = &(in[0]);
    
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = out_eT( in_mem[i] );
      }
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  Mat<out_eT> out(in.size(), 1);
  
  typename std::vector<in_eT>::const_iterator in_begin = in.begin();
  typename std::vector<in_eT>::const_iterator in_end   = in.end();
  
  typename Mat<out_eT>::iterator out_begin = out.begin();
  typename Mat<out_eT>::iterator out_end   = out.end();
  
  typename std::vector<in_eT>::const_iterator in_it;
  typename Mat<out_eT>::iterator              out_it;
  
  for(in_it = in_begin, out_it = out_begin; (in_it != in_end) && (out_it != out_end); ++in_it, ++out_it)
    {
          out_eT& out_elem = (*out_it);
    const in_eT&  in_elem  = (*in_it);
    
    arrayops::convert_cx_scalar(out_elem, in_elem);
    }
  
  return out;
  }



//! conversion to Armadillo row vectors from Armadillo Base objects, as well as from std::vector
template<typename out_eT>
class conv_to< Row<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Row<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static Row<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  
  
  
  template<typename in_eT>
  inline static Row<out_eT> from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT>
  inline static Row<out_eT> from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object can't be interpreted as a vector" );
  
  Row<out_eT> out(X.n_elem);
  
  arrayops::convert( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object can't be interpreted as a vector" );
  
  Row<out_eT> out(X.n_rows, X.n_cols);
  
  arrayops::convert_cx( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Row<out_eT> out(N);
  
  if(N > 0)
    {
          out_eT* out_mem = out.memptr();
    const  in_eT*  in_mem = &(in[0]);
    
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = out_eT( in_mem[i] );
      }
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  Row<out_eT> out( in.size() );
  
  typename std::vector<in_eT>::const_iterator in_begin = in.begin();
  typename std::vector<in_eT>::const_iterator in_end   = in.end();
  
  typename Row<out_eT>::iterator out_begin = out.begin();
  typename Row<out_eT>::iterator out_end   = out.end();
  
  typename std::vector<in_eT>::const_iterator in_it;
  typename Row<out_eT>::iterator              out_it;
  
  for(in_it = in_begin, out_it = out_begin; (in_it != in_end) && (out_it != out_end); ++in_it, ++out_it)
    {
          out_eT& out_elem = (*out_it);
    const in_eT&  in_elem  = (*in_it);
    
    arrayops::convert_cx_scalar(out_elem, in_elem);
    }
  
  return out;
  }



//! conversion to Armadillo column vectors from Armadillo Base objects, as well as from std::vector
template<typename out_eT>
class conv_to< Col<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Col<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static Col<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  
  
  
  template<typename in_eT>
  inline static Col<out_eT> from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT>
  inline static Col<out_eT> from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object can't be interpreted as a vector" );
  
  Col<out_eT> out(X.n_elem);
  
  arrayops::convert( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object can't be interpreted as a vector" );
  
  Col<out_eT> out(X.n_rows, X.n_cols);
  
  arrayops::convert_cx( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Col<out_eT> out(N);
  
  if(N > 0)
    {
          out_eT* out_mem = out.memptr();
    const  in_eT*  in_mem = &(in[0]);
    
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = out_eT( in_mem[i] );
      }
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  Col<out_eT> out( in.size() );
  
  typename std::vector<in_eT>::const_iterator in_begin = in.begin();
  typename std::vector<in_eT>::const_iterator in_end   = in.end();
  
  typename Col<out_eT>::iterator out_begin = out.begin();
  typename Col<out_eT>::iterator out_end   = out.end();
  
  typename std::vector<in_eT>::const_iterator in_it;
  typename Col<out_eT>::iterator              out_it;
  
  for(in_it = in_begin, out_it = out_begin; (in_it != in_end) && (out_it != out_end); ++in_it, ++out_it)
    {
          out_eT& out_elem = (*out_it);
    const in_eT&  in_elem  = (*in_it);
    
    arrayops::convert_cx_scalar(out_elem, in_elem);
    }
  
  return out;
  }



//! conversion to Armadillo cubes from Armadillo BaseCube objects
template<typename out_eT>
class conv_to< Cube<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Cube<out_eT> from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static Cube<out_eT> from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_cube<T1>  tmp( in.get_ref() );
  const Cube<in_eT>& X = tmp.M;
  
  Cube<out_eT> out(X.n_rows, X.n_cols, X.n_slices);
  
  arrayops::convert( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_cube<T1>  tmp( in.get_ref() );
  const Cube<in_eT>& X = tmp.M;
  
  Cube<out_eT> out(X.n_rows, X.n_cols, X.n_slices);
  
  arrayops::convert_cx( out.memptr(), X.memptr(), out.n_elem );
  
  return out;
  }



//! conversion to std::vector from Armadillo Base objects
template<typename out_eT>
class conv_to< std::vector<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static std::vector<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = 0);
  
  template<typename in_eT, typename T1>
  inline static std::vector<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = 0);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
inline
std::vector<out_eT>
conv_to< std::vector<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object can't be interpreted as a vector" );
  
  const uword N = X.n_elem;
  
  std::vector<out_eT> out(N);
  
  if(N > 0)
    {
          out_eT* out_mem = &(out[0]);
    const  in_eT*   X_mem = X.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = out_eT( X_mem[i] );
      }
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
std::vector<out_eT>
conv_to< std::vector<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap<T1>      tmp(in.get_ref());
  const Mat<in_eT>& X = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object can't be interpreted as a vector" );
  
  std::vector<out_eT> out(X.n_elem);
  
  typename Mat<in_eT>::const_iterator X_begin = X.begin();
  typename Mat<in_eT>::const_iterator X_end   = X.end();

  typename std::vector<out_eT>::iterator out_begin = out.begin();
  typename std::vector<out_eT>::iterator out_end   = out.end();
  
  typename Mat<in_eT>::const_iterator    X_it;
  typename std::vector<out_eT>::iterator out_it;
  
  for(X_it = X_begin, out_it = out_begin; (X_it != X_end) && (out_it != out_end); ++X_it, ++out_it)
    {
          out_eT& out_elem = (*out_it);
    const in_eT&  X_elem   = (*X_it);
    
    arrayops::convert_cx_scalar(out_elem, X_elem);
    }
  
  return out;
  }



//! @}

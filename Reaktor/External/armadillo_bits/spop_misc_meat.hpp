// Copyright (C) 2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup spop_misc
//! @{



namespace priv
  {
  template<typename eT>
  struct functor_scalar_times
    {
    const eT k;
    
    functor_scalar_times(const eT in_k) : k(in_k) {}
    
    arma_inline eT operator()(const eT val) const { return val * k; }
    };
  }



template<typename T1>
inline
void
spop_scalar_times::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(in.aux != eT(0))
    {
    out.init_xform(in.m, priv::functor_scalar_times<eT>(in.aux));
    }
  else
    {
    out.reset();
    
    const SpProxy<T1> P(in.m);
    
    out.set_size( P.get_n_rows(), P.get_n_cols() );
    }
  }



namespace priv
  {
  struct functor_square
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return val*val; }
    };
  }



template<typename T1>
inline
void
spop_square::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_square>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_square());
  }



namespace priv
  {
  struct functor_sqrt
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return eop_aux::sqrt(val); }
    };
  }



template<typename T1>
inline
void
spop_sqrt::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_sqrt>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_sqrt());
  }



namespace priv
  {
  struct functor_abs
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return eop_aux::arma_abs(val); }
    };
  }



template<typename T1>
inline
void
spop_abs::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_abs>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_abs());
  }



namespace priv
  {
  struct functor_cx_abs
    {
    template<typename T>
    arma_inline T operator()(const std::complex<T>& val) const { return std::abs(val); }
    };
  }



template<typename T1>
inline
void
spop_cx_abs::apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_abs>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform_mt(in.m, priv::functor_cx_abs());
  }



//! @}

// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_vectorise
//! @{



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_vectorise_col>
  >::result
vectorise(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_vectorise_col>(X);
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_vectorise_all>
  >::result
vectorise(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (dim > 1), "vectorise(): dim must be 0 or 1");
  
  return Op<T1, op_vectorise_all>(X, dim, 0);
  }



//! @}

// Copyright (C) 2013-2014 Ryan Curtin
// Copyright (C) 2013-2014 Conrad Sanderson
// Copyright (C) 2013-2014 NICTA
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup sp_auxlib
//! @{


inline
sp_auxlib::form_type
sp_auxlib::interpret_form_str(const char* form_str)
  {
  arma_extra_debug_sigprint();
  
  // the order of the 3 if statements below is important
  if( form_str    == NULL    )  { return form_none; }
  if( form_str[0] == char(0) )  { return form_none; }
  if( form_str[1] == char(0) )  { return form_none; }
  
  const char c1 = form_str[0];
  const char c2 = form_str[1];
  
  if(c1 == 'l')
    {
    if(c2 == 'm')  { return form_lm; }
    if(c2 == 'r')  { return form_lr; }
    if(c2 == 'i')  { return form_li; }
    }
  else
  if(c1 == 's')
    {
    if(c2 == 'm')  { return form_sm; }
    if(c2 == 'r')  { return form_sr; }
    if(c2 == 'i')  { return form_si; }
    }
  
  return form_none;
  }



//! immediate eigendecomposition of symmetric real sparse object
template<typename eT, typename T1>
inline
bool
sp_auxlib::eigs_sym(Col<eT>& eigval, Mat<eT>& eigvec, const SpBase<eT, T1>& X, const uword n_eigvals, const char* form_str, const eT default_tol)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ARPACK)
    {
    const form_type form_val = sp_auxlib::interpret_form_str(form_str);
    
    arma_debug_check( (form_val != form_lm) && (form_val != form_sm), "eigs_sym(): unknown form specified" );
    
    char  which_sm[3] = "SM";
    char  which_lm[3] = "LM";
    char* which      = (form_val == form_sm) ? which_sm : which_lm;  // select which eigenvalues we want: smallest magnitude or largest magnitude
    
    // Make a sparse proxy object.
    SpProxy<T1> p(X.get_ref());
    
    // Make sure it's square.
    arma_debug_check( (p.get_n_rows() != p.get_n_cols()), "eigs_sym(): given sparse matrix is not square");
    
    // Make sure we aren't asking for every eigenvalue.
    arma_debug_check( (n_eigvals + 1 >= p.get_n_rows()), "eigs_sym(): n_eigvals + 1 must be less than the number of rows in the matrix");
    
    // If the matrix is empty, the case is trivial.
    if(p.get_n_cols() == 0) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    // Set up variables that get used for neupd().
    blas_int n, ncv, ldv, lworkl, info;
    eT tol = default_tol;
    podarray<eT> resid, v, workd, workl;
    podarray<blas_int> iparam, ipntr;
    podarray<eT> rwork; // Not used in this case.
    
    run_aupd(n_eigvals, which, p, true /* sym, not gen */, n, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
    
    if(info != 0)
      {
      return false;
      }
    
    // The process has converged, and now we need to recover the actual eigenvectors using seupd()
    blas_int rvec = 1; // .TRUE
    blas_int nev  = n_eigvals;
    
    char howmny = 'A';
    char bmat   = 'I'; // We are considering the standard eigenvalue problem.
    
    podarray<blas_int> select(ncv); // Logical array of dimension NCV.
    blas_int ldz = n;
    
    // seupd() will output directly into the eigval and eigvec objects.
    eigval.set_size(n_eigvals);
    eigvec.set_size(n, n_eigvals);
    
    arpack::seupd(&rvec, &howmny, select.memptr(), eigval.memptr(), eigvec.memptr(), &ldz, (eT*) NULL, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, &info);
    
    // Check for errors.
    if(info != 0)
      {
      std::stringstream tmp;
      tmp << "eigs_sym(): ARPACK error " << info << " in seupd()";
      arma_debug_warn(true, tmp.str());
      return false;
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_stop("eigs_sym(): use of ARPACK needs to be enabled");
    return false;
    }
  #endif
  }



//! immediate eigendecomposition of non-symmetric real sparse object
template<typename T, typename T1>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase<T, T1>& X, const uword n_eigvals, const char* form_str, const T default_tol)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ARPACK)
    {
    const form_type form_val = sp_auxlib::interpret_form_str(form_str);
    
    arma_debug_check( (form_val == form_none), "eigs_gen(): unknown form specified" );
    
    char which_lm[3] = "LM";
    char which_sm[3] = "SM";
    char which_lr[3] = "LR";
    char which_sr[3] = "SR";
    char which_li[3] = "LI";
    char which_si[3] = "SI";
    
    char* which;
    
    switch(form_val)
      {
      case form_lm:  which = which_lm;  break;
      case form_sm:  which = which_sm;  break;
      case form_lr:  which = which_lr;  break;
      case form_sr:  which = which_sr;  break;
      case form_li:  which = which_li;  break;
      case form_si:  which = which_si;  break;
      
      default:       which = which_lm;
      }
    
    
    // Make a sparse proxy object.
    SpProxy<T1> p(X.get_ref());
    
    // Make sure it's square.
    arma_debug_check( (p.get_n_rows() != p.get_n_cols()), "eigs_gen(): given sparse matrix is not square");
    
    // Make sure we aren't asking for every eigenvalue.
    arma_debug_check( (n_eigvals + 1 >= p.get_n_rows()), "eigs_gen(): n_eigvals + 1 must be less than the number of rows in the matrix");
    
    // If the matrix is empty, the case is trivial.
    if(p.get_n_cols() == 0) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    // Set up variables that get used for neupd().
    blas_int n, ncv, ldv, lworkl, info;
    T tol = default_tol;
    podarray<T> resid, v, workd, workl;
    podarray<blas_int> iparam, ipntr;
    podarray<T> rwork; // Not used in the real case.
    
    run_aupd(n_eigvals, which, p, false /* gen, not sym */, n, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
    
    if(info != 0)
      {
      return false;
      }

    // The process has converged, and now we need to recover the actual eigenvectors using neupd().
    blas_int rvec = 1; // .TRUE
    blas_int nev = n_eigvals;
    
    char howmny = 'A';
    char bmat   = 'I'; // We are considering the standard eigenvalue problem.
    
    podarray<blas_int> select(ncv); // Logical array of dimension NCV.
    podarray<T> dr(nev + 1); // Real array of dimension NEV + 1.
    podarray<T> di(nev + 1); // Real array of dimension NEV + 1.
    podarray<T> z(n * (nev + 1)); // Real N by NEV array if HOWMNY = 'A'.
    blas_int ldz = n;
    podarray<T> workev(3 * ncv);
    
    arpack::neupd(&rvec, &howmny, select.memptr(), dr.memptr(), di.memptr(), z.memptr(), &ldz, (T*) NULL, (T*) NULL, workev.memptr(), &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
    
    // Check for errors.
    if(info != 0)
      {
      std::stringstream tmp;
      tmp << "eigs_gen(): ARPACK error " << info << " in neupd()";
      arma_debug_warn(true, tmp.str());
      return false;
      }
    
    // Put it into the outputs.
    eigval.set_size(n_eigvals);
    eigvec.set_size(n, n_eigvals);
    
    for (uword i = 0; i < n_eigvals; ++i)
      {
      eigval[i] = std::complex<T>(dr[i], di[i]);
      }
    
    // Now recover the eigenvectors.
    for (uword i = 0; i < n_eigvals; ++i)
      {
      // ARPACK ?neupd lays things out kinda odd in memory; so does LAPACK
      // ?geev (see auxlib::eig_gen()).
      if((i < n_eigvals - 1) && (eigval[i] == std::conj(eigval[i + 1])))
        {
        for (uword j = 0; j < n; ++j)
          {
          eigvec.at(j, i)     = std::complex<T>(z[n * i + j], z[n * (i + 1) + j]);
          eigvec.at(j, i + 1) = std::complex<T>(z[n * i + j], -z[n * (i + 1) + j]);
          }
        ++i; // Skip the next one.
        }
      else
      if((i == n_eigvals - 1) && (std::complex<T>(eigval[i]).imag() != 0.0))
        {
        // We don't have the matched conjugate eigenvalue.
        for (uword j = 0; j < n; ++j)
          {
          eigvec.at(j, i) = std::complex<T>(z[n * i + j], z[n * (i + 1) + j]);
          }
        }
      else
        {
        // The eigenvector is entirely real.
        for (uword j = 0; j < n; ++j)
          {
          eigvec.at(j, i) = std::complex<T>(z[n * i + j], T(0));
          }
        }
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_stop("eigs_gen(): use of ARPACK needs to be enabled");
    return false;
    }
  #endif
  }




//! immediate eigendecomposition of non-symmetric complex sparse object
template<typename T, typename T1>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase< std::complex<T>, T1>& X, const uword n_eigvals, const char* form_str, const T default_tol)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ARPACK)
    {
    const form_type form_val = sp_auxlib::interpret_form_str(form_str);
    
    arma_debug_check( (form_val == form_none), "eigs_gen(): unknown form specified" );
    
    char which_lm[3] = "LM";
    char which_sm[3] = "SM";
    char which_lr[3] = "LR";
    char which_sr[3] = "SR";
    char which_li[3] = "LI";
    char which_si[3] = "SI";
    
    char* which;
    
    switch(form_val)
      {
      case form_lm:  which = which_lm;  break;
      case form_sm:  which = which_sm;  break;
      case form_lr:  which = which_lr;  break;
      case form_sr:  which = which_sr;  break;
      case form_li:  which = which_li;  break;
      case form_si:  which = which_si;  break;
      
      default:       which = which_lm;
      }
    
    
    // Make a sparse proxy object.
    SpProxy<T1> p(X.get_ref());
    
    // Make sure it's square.
    arma_debug_check( (p.get_n_rows() != p.get_n_cols()), "eigs_gen(): given sparse matrix is not square");
    
    // Make sure we aren't asking for every eigenvalue.
    arma_debug_check( (n_eigvals + 1 >= p.get_n_rows()), "eigs_gen(): n_eigvals + 1 must be less than the number of rows in the matrix");
    
    // If the matrix is empty, the case is trivial.
    if(p.get_n_cols() == 0) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    // Set up variables that get used for neupd().
    blas_int n, ncv, ldv, lworkl, info;
    T tol = default_tol;
    podarray< std::complex<T> > resid, v, workd, workl;
    podarray<blas_int> iparam, ipntr;
    podarray<T> rwork;
    
    run_aupd(n_eigvals, which, p, false /* gen, not sym */, n, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
    
    if(info != 0)
      {
      return false;
      }
    
    // The process has converged, and now we need to recover the actual eigenvectors using neupd().
    blas_int rvec = 1; // .TRUE
    blas_int nev  = n_eigvals;
    
    char howmny = 'A';
    char bmat   = 'I'; // We are considering the standard eigenvalue problem.
    
    podarray<blas_int> select(ncv); // Logical array of dimension NCV.
    podarray<std::complex<T> > d(nev + 1); // Real array of dimension NEV + 1.
    podarray<std::complex<T> > z(n * nev); // Real N by NEV array if HOWMNY = 'A'.
    blas_int ldz = n;
    podarray<std::complex<T> > workev(2 * ncv);
    
    // Prepare the outputs; neupd() will write directly to them.
    eigval.set_size(n_eigvals);
    eigvec.set_size(n, n_eigvals);
    std::complex<T> sigma;
    
    arpack::neupd(&rvec, &howmny, select.memptr(), eigval.memptr(),
(std::complex<T>*) NULL, eigvec.memptr(), &ldz, (std::complex<T>*) &sigma, (std::complex<T>*) NULL, workev.memptr(), &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
    
    // Check for errors.
    if(info != 0)
      {
      std::stringstream tmp;
      tmp << "eigs_gen(): ARPACK error " << info << " in neupd()";
      arma_debug_warn(true, tmp.str());
      return false;
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_stop("eigs_gen(): use of ARPACK needs to be enabled");
    return false;
    }
  #endif
  }



template<typename eT, typename T, typename T1>
inline
void
sp_auxlib::run_aupd
  (
  const uword n_eigvals, char* which, const SpProxy<T1>& p, const bool sym,
  blas_int& n, eT& tol,
  podarray<T>& resid, blas_int& ncv, podarray<T>& v, blas_int& ldv,
  podarray<blas_int>& iparam, podarray<blas_int>& ipntr,
  podarray<T>& workd, podarray<T>& workl, blas_int& lworkl, podarray<eT>& rwork,
  blas_int& info
  )
  {
  #if defined(ARMA_USE_ARPACK)
    {
    // ARPACK provides a "reverse communication interface" which is an
    // entertainingly archaic FORTRAN software engineering technique that
    // basically means that we call saupd()/naupd() and it tells us with some
    // return code what we need to do next (usually a matrix-vector product) and
    // then call it again.  So this results in some type of iterative process
    // where we call saupd()/naupd() many times.
    blas_int ido = 0; // This must be 0 for the first call.
    char bmat = 'I'; // We are considering the standard eigenvalue problem.
    n = p.get_n_rows(); // The size of the matrix.
    blas_int nev = n_eigvals;
    
    resid.set_size(n);
    
    // "NCV must satisfy the two inequalities 2 <= NCV-NEV and NCV <= N".
    // "It is recommended that NCV >= 2 * NEV".
    ncv = 2 + nev;
    if (ncv < 2 * nev) { ncv = 2 * nev; }
    if (ncv > n)       { ncv = n; }
    v.set_size(n * ncv); // Array N by NCV (output).
    rwork.set_size(ncv); // Work array of size NCV for complex calls.
    ldv = n; // "Leading dimension of V exactly as declared in the calling program."
    
    // IPARAM: integer array of length 11.
    iparam.zeros(11);
    iparam(0) = 1; // Exact shifts (not provided by us).
    iparam(2) = 1000; // Maximum iterations; all the examples use 300, but they were written in the ancient times.
    iparam(6) = 1; // Mode 1: A * x = lambda * x.
    
    // IPNTR: integer array of length 14 (output).
    ipntr.set_size(14);
    
    // Real work array used in the basic Arnoldi iteration for reverse communication.
    workd.set_size(3 * n);
    
    // lworkl must be at least 3 * NCV^2 + 6 * NCV.
    lworkl = 3 * (ncv * ncv) + 6 * ncv;
    
    // Real work array of length lworkl.
    workl.set_size(lworkl);
    
    info = 0; // Set to 0 initially to use random initial vector.
    
    // All the parameters have been set or created.  Time to loop a lot.
    while (ido != 99)
      {
      // Call saupd() or naupd() with the current parameters.
      if(sym)
        arpack::saupd(&ido, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, &info);
      else
        arpack::naupd(&ido, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
      
      // What do we do now?
      switch (ido)
        {
        case -1:
        case 1:
          {
          // We need to calculate the matrix-vector multiplication y = OP * x
          // where x is of length n and starts at workd(ipntr(0)), and y is of
          // length n and starts at workd(ipntr(1)).
          
          // operator*(sp_mat, vec) doesn't properly put the result into the
          // right place so we'll just reimplement it here for now...
          
          // Set the output to point at the right memory.  We have to subtract
          // one from FORTRAN pointers...
          Col<T> out(workd.memptr() + ipntr(1) - 1, n, false /* don't copy */);
          // Set the input to point at the right memory.
          Col<T> in(workd.memptr() + ipntr(0) - 1, n, false /* don't copy */);
          
          out.zeros();
          typename SpProxy<T1>::const_iterator_type x_it     = p.begin();
          typename SpProxy<T1>::const_iterator_type x_it_end = p.end();
          
          while(x_it != x_it_end)
            {
            out[x_it.row()] += (*x_it) * in[x_it.col()];
            ++x_it;
            }
          
          // No need to modify memory further since it was all done in-place.
          
          break;
          }
        case 99:
          // Nothing to do here, things have converged.
          break;
        default:
          {
          return; // Parent frame can look at the value of info.
          }
        }
      }
    
    // The process has ended; check the return code.
    if( (info != 0) && (info != 1) )
      {
      // Print warnings if there was a failure.
      std::stringstream tmp;
      
      if(sym)
        {
        tmp << "eigs_sym(): ARPACK error " << info << " in saupd()";
        }
      else
        {
        tmp << "eigs_gen(): ARPACK error " << info << " in naupd()";
        }
      
      arma_debug_warn(true, tmp.str());
      
      return; // Parent frame can look at the value of info.
      }
    }
  #endif
  }

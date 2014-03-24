#ifndef __TTNS_DAVIDSON_H
#define __TTNS_DAVIDSON_H 1

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <functional>

#include <btas/QSPARSE/QSDArray.h>

#include "Input.h"

namespace ttns
{

namespace davidson
{

template<size_t N, class Q>
using Functor = std::function<void(const btas::QSDArray<N, Q>&, btas::QSDArray<N, Q>&)>;

//
// Davidson's precondition
//

template<size_t N, class Q>
void precondition (const double& eval, const btas::QSDArray<N, Q>& diag, btas::QSDArray<N, Q>& errv)
{
   for(auto ir = errv.begin(); ir != errv.end(); ++ir)
   {
      auto id = diag.find(ir->first);
      if(id != diag.end())
      {
         auto irx = ir->second->begin();
         auto idx = id->second->begin();
         for(; irx != ir->second->end(); ++irx, ++idx)
         {
            double denm = eval - *idx;
            if(fabs(denm) < 1.0e-12) denm = 1.0e-12;
            *irx /= denm;
         }
      }
      else
      {
         btas::Dscal(1.0/eval, (*ir->second));
      }
   }
}

//
// Davidson eigen solver
//

template<size_t N, class Q>
double diagonalize (
      const Functor<N, Q>& f_contract,
      const btas::QSDArray<N, Q>& diag,
            btas::QSDArray<N, Q>& wfnc,
      const SweepParameters& param)
{
   int max_ritz = param.MaxDimDavidson();

   double eval = 0.0;

   // reserve working space
   std::vector<btas::QSDArray<N, Q>> trial(max_ritz);
   std::vector<btas::QSDArray<N, Q>> sigma(max_ritz);

   btas::QSDcopy(wfnc, trial[0]);

   int niter = 0;
   int iconv = 0;
   std::cout << "\t\t\t----------------------------------------------------------------" << std::endl;
   while(iconv < 1 && niter < param.MaxIterDavidson())
   {
      std::cout << "\t\t\t\tmacro iteration [ " << std::setw(2) << niter << " ] " << std::endl;
      std::cout << "\t\t\t----------------------------------------------------------------" << std::endl;
      // to keep numerical stability
      btas::Normalize(trial[0]);
      sigma[0].clear();
      f_contract(trial[0], sigma[0]);
      for(int m = 1; m <= max_ritz; ++m)
      {
         // compute small Hamiltonian matrix
         btas::DArray<2> heff(m, m);
         btas::DArray<2> ovlp(m, m);
         for(int i = 0; i < m; ++i)
         {
            heff(i, i) = btas::QSDdotc(trial[i], sigma[i]);
            ovlp(i, i) = btas::QSDdotc(trial[i], trial[i]);
            for(int j = 0; j < i; ++j)
            {
               double hij = btas::QSDdotc(trial[i], sigma[j]);
               heff(i, j) = hij;
               heff(j, i) = hij;
               double sij = btas::QSDdotc(trial[i], trial[j]);
               ovlp(i, j) = sij;
               ovlp(j, i) = sij;
            }
         }
         // solve eigenvalue problem to obtain Ritz value & vector
         btas::DArray<2> rvec;
         btas::DArray<1> rval;
         Dsyev('V', 'U', heff, rval, rvec);
         eval = rval(0);
         std::cout << "\t\t\t\tmicro iteration [ " << std::setw(2) << m << " ] :: " << std::setprecision(16) << std::setw(22) << std::fixed << eval << std::endl;
         // rotate trial & sigma vectors by Ritz vector
         std::vector<btas::QSDArray<N, Q>> trial_save(m);
         std::vector<btas::QSDArray<N, Q>> sigma_save(m);
         for(int i = 0; i < m; ++i)
         {
            btas::QSDcopy(trial[i], trial_save[i]);
            btas::QSDcopy(sigma[i], sigma_save[i]);
            btas::QSDscal(rvec(i, i), trial[i]);
            btas::QSDscal(rvec(i, i), sigma[i]);
         }
         for(int i = 0; i < m; ++i)
         {
            for(int j = 0; j < m; ++j)
            {
               if(i != j)
               {
                  btas::QSDaxpy(rvec(i, j), trial_save[i], trial[j]);
                  btas::QSDaxpy(rvec(i, j), sigma_save[i], sigma[j]);
               }
            }
         }
         // compute error vector
         btas::QSDArray<N, Q> evec;
         btas::QSDArray<N, Q> errv;
         btas::QSDcopy( trial[0], evec);
         btas::QSDcopy( sigma[0], errv);
         btas::QSDaxpy(-eval, evec, errv);
         double rnorm = btas::QSDdotc(errv, errv);
         if(rnorm < param.ToleDavidson()) { ++iconv; break; }
         // solve correction equation
         if(m < max_ritz)
         {
            precondition(eval, diag, errv);
            for(int i = 0; i < m; ++i)
            {
               btas::Normalize(errv);
               btas::Orthogonalize(trial[i], errv);
            }
            btas::Normalize(errv);
            btas::QSDcopy(errv, trial[m]);
            sigma[m].clear();
            f_contract(trial[m], sigma[m]);
         }
      }
      ++niter;
      std::cout << "\t\t\t----------------------------------------------------------------" << std::endl;
   }
   btas::QSDcopy(trial[0], wfnc);

   return eval;
}

} // namespace davidson

} // namespace ttns

#endif // __TTNS_DAVIDSON_H

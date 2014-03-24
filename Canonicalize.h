#ifndef __TTNS_CANONICALIZE_H
#define __TTNS_CANONICALIZE_H

#include <btas/QSPARSE/QSDArray.h>
#include "ttns_assert.h"

namespace ttns
{

/// tree graph direction
enum Direction { Forward, Backward };

/// von-Neumann entropy
inline double vonNeumannEntropy (const btas::SDArray<1>& sval)
{
   double entropy = 0.0;

   for(auto it = sval.begin(); it != sval.end(); ++it)
   {
      for(double& si : *it->second) entropy += (si*si*log(si));
   }

   return -2.0*entropy/log(2.0);
// return -2.0*entropy*log(2.0); // ?
}

/// Canonicalize for onedot
template<size_t N, class Q>
double CanonicalizeOneDot (btas::QSDArray<N, Q>& tns, size_t iSvd, size_t D, btas::QSDArray<2, Q>& res)
{
   TTNS_DEBUG("CanonicalizeOneDot : Started.");

   btas::QSDArray<N, Q> tp;

   btas::IVector<N> reorder;
   for(size_t i = 0; i < iSvd; ++i) reorder[i] = i;
   for(size_t i = iSvd; i < N-1; ++i) reorder[i] = i+1;
   reorder[N-1] = iSvd;

   btas::Permute(tns, reorder, tp);

   btas::SDArray<1> sval;
   btas::QSDArray<N, Q> up;

   btas::Gesvd(tp, sval, up, res, D);

   // FIXME: this should be done before SVD on tp?
   // parity operation for total quanta
   if(tns.q().parity() && iSvd < N-1)
   {
      std::vector<int> pi;
      for(size_t i = 0; i < iSvd; ++i) pi.push_back(i);
      up.parity(pi);
   }

   if(iSvd == N-1)
   {
      tns = up;
   }
   else
   {
      btas::IVector<N> back_reorder;
      for(size_t i = 0; i < iSvd; ++i) back_reorder[i] = i;
      for(size_t i = iSvd+1; i < N; ++i) back_reorder[i] = i-1;
      back_reorder[iSvd] = N-1;
      btas::Permute(up, back_reorder, tns);
   }

   btas::Dimm(sval, res);

   TTNS_DEBUG("CanonicalizeOneDot : Finished.");

   return vonNeumannEntropy(sval);
}

/// Canonicalize for twodot
template<size_t M, size_t N, class Q>
double CanonicalizeTwoDot (const btas::QSDArray<M+N-2, Q>& psi, btas::QSDArray<M, Q>& tnsA, btas::QSDArray<N, Q>& tnsB, size_t indexB, size_t D, Direction dir)
{
   TTNS_DEBUG("CanonicalizeTwoDot : Started.");

   btas::SDArray<1> sval;

   if(dir == Forward)
   {
      btas::QSDArray<N, Q> tnsBp;

      btas::Gesvd(psi, sval, tnsA, tnsBp, D);

      btas::Dimm(sval, tnsBp);

      btas::IVector<N> back_reorderB;
      back_reorderB[indexB] = 0;
      for(size_t i = 0; i < indexB; ++i) back_reorderB[i] = i+1;
      for(size_t i = indexB+1; i < N; ++i) back_reorderB[i] = i;

      btas::Permute(tnsBp, back_reorderB, tnsB);
   }
   else
   {
      btas::QSDArray<M, Q> tnsAp;
      btas::QSDArray<N, Q> tnsBp;

      {
         btas::QSDArray<M+N-2, Q> pp;
         btas::IVector<M+N-2> reorder;
         for(size_t i = 0; i < N-1; ++i) reorder[i] = i+M-1;
         for(size_t i = 0; i < M-1; ++i) reorder[i+N-1] = i;
         btas::Permute(psi, reorder, pp);

         btas::Gesvd(pp, sval, tnsBp, tnsAp, D);
      }

      btas::Dimm(sval, tnsAp);

      btas::IVector<M> back_reorderA;
      back_reorderA[M-1] = 0;
      for(size_t i = 0; i < M-1; ++i) back_reorderA[i] = i+1;
      btas::Permute(tnsAp, back_reorderA, tnsA);

      btas::IVector<N> back_reorderB;
      back_reorderB[indexB] = N-1;
      for(size_t i = 0; i < indexB; ++i) back_reorderB[i] = i;
      for(size_t i = indexB+1; i < N; ++i) back_reorderB[i] = i-1;
      btas::Permute(tnsBp, back_reorderB, tnsB);
   }

   TTNS_DEBUG("CanonicalizeTwoDot : Finished.");

   return vonNeumannEntropy(sval);
}

} // namespace ttns

#endif // __TTNS_CANONICALIZE_H

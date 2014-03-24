#ifndef __TTNS_OP_COMPONENTS_UTIL_H
#define __TTNS_OP_COMPONENTS_UTIL_H

#include <btas/common/TVector.h>
#include <btas/SPARSE/SDArray.h>

#include "ttns_assert.h"

namespace ttns
{

template<size_t N>
inline void ScaledCopy (
      const double& alpha,
      const btas::SDArray<N>& x,
            btas::SDArray<N>& y)
{
// TTNS_DEBUG("ScaledCopy : Started.");

   btas::Copy(x, y);
   btas::Scal(alpha, y);

// TTNS_DEBUG("ScaledCopy : Finished.");
}

template<size_t N>
inline void PermutedAxpy (
      const double& alpha,
      const btas::SDArray<N>& x,
      const btas::IVector<N>& reorder,
            btas::SDArray<N>& y)
{
// TTNS_DEBUG("PermutedAxpy : Started.");

   if(x.size() > 0)
   {
      btas::SDArray<N> x_scr;
      btas::Permute(x, reorder, x_scr);

      btas::Axpy(alpha, x_scr, y);
   }

// TTNS_DEBUG("PermutedAxpy : Finished.");
}

//
// (1) computing sigma vector: N = K
// (2) renormalization       : N = K + 1
//
//     N : rank of wavefunction
//     K : # blocks to be contracted
//

template<size_t N>
inline void IndexPermute (
            btas::IVector<N>& reorder,
      const size_t& op1)
{
   for(size_t i = 0;     i < op1; ++i) reorder[i] = i+1;
   reorder[op1] = 0;
   for(size_t i = op1+1; i < N;   ++i) reorder[i] = i;
}

template<size_t N>
inline void IndexPermute (
            btas::IVector<N>& reorder,
      const size_t& op1,
      const size_t& op2)
{
   for(size_t i = 0; i < N;   ++i) reorder[i] = i;
   for(size_t i = 0; i < op1; ++i) ++reorder[i];
   for(size_t i = 0; i < op2; ++i) ++reorder[i];
   reorder[op1] = 1;
   reorder[op2] = 0;
}

template<size_t N>
inline void IndexPermute (
            btas::IVector<N>& reorder,
      const size_t& op1,
      const size_t& op2,
      const size_t& op3)
{
   for(size_t i = 0; i < N;   ++i) reorder[i] = i;
   for(size_t i = 0; i < op1; ++i) ++reorder[i];
   for(size_t i = 0; i < op2; ++i) ++reorder[i];
   for(size_t i = 0; i < op3; ++i) ++reorder[i];
   reorder[op1] = 2;
   reorder[op2] = 1;
   reorder[op3] = 0;
}

template<size_t N>
inline void IndexPermute (
            btas::IVector<N>& reorder,
      const size_t& op1,
      const size_t& op2,
      const size_t& op3,
      const size_t& op4)
{
   for(size_t i = 0; i < N;   ++i) reorder[i] = i;
   for(size_t i = 0; i < op1; ++i) ++reorder[i];
   for(size_t i = 0; i < op2; ++i) ++reorder[i];
   for(size_t i = 0; i < op3; ++i) ++reorder[i];
   for(size_t i = 0; i < op4; ++i) ++reorder[i];
   reorder[op1] = 3;
   reorder[op2] = 2;
   reorder[op3] = 1;
   reorder[op4] = 0;
}

//
// Create parity indices
//
template<size_t N>
inline void IndexParity (
            double& opSign,
            std::vector<int>& parIndex,
      const btas::IVector<N>& order,
      const size_t& op1,
      const size_t& op2)
{
   parIndex.clear();
   size_t lb = order[op1]; // lower bound
   size_t ub = order[op2]; // upper bound

   if(lb > ub)
   {
      opSign = -opSign;
      std::swap(lb, ub);
   }

   for(size_t i = 0; i < N; ++i)
      if(lb < order[i] && order[i] <= ub) parIndex.push_back(i);
}

} // namespace ttns

#endif // __TTNS_OP_COMPONENTS_UTIL_H

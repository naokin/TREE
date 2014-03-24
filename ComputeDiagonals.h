#ifndef __TTNS_COMPUTE_DIAGONALS_H
#define __TTNS_COMPUTE_DIAGONALS_H

#include "Block.h"
#include "ttns_assert.h"

namespace ttns
{

inline void DiagonalElement (const btas::SDArray<2>& m, btas::SDArray<1>& d)
{
   btas::Tie(m, btas::shape(0, 1), d);
}

template<size_t N, class Q>
void Get_Ham_diag (
      const btas::SDArray<1>& op_Ham_diag, size_t op1,
            btas::QSDArray<N, Q>& diag)
{
   TTNS_DEBUG("Get_Ham_diag : Started.");
// TTNS_DEBUG("Get_Ham_diag : diag(input) : " << diag);

   size_t str_00f = diag.stride(op1);
   size_t str_0ff = str_00f*diag.shape(op1);

   size_t ext_i = diag.size()/str_0ff;
   size_t ext_k = str_00f;

   for(auto it1 = op_Ham_diag.begin(); it1 != op_Ham_diag.end(); ++it1)
   {
      size_t idx_0j0 = it1->first*str_00f;

      for(size_t idx_i00 = 0; idx_i00 < ext_i; ++idx_i00)
      {
         size_t idx_ij0 = idx_i00*str_0ff+idx_0j0;

         for(size_t idx_00k = 0; idx_00k < str_00f; ++idx_00k)
         {
            size_t idx_ijk = idx_ij0+idx_00k;

            auto it2 = diag.reserve(idx_ijk);

            if(it2 == diag.end()) continue;

            btas::DArray<1>& x = *it1->second;
            btas::DArray<N>& y = *it2->second;

            size_t d_ext_k = y.stride(op1);
            size_t d_ext_j = x.size();
            size_t d_ext_i = y.size()/d_ext_j/d_ext_k;

            double* py = y.data();
            for(size_t d_i = 0; d_i < d_ext_i; ++d_i)
            {
               const double* px = x.data();
               for(size_t d_j = 0; d_j < d_ext_j; ++d_j, ++px)
               {
                  for(size_t d_k = 0; d_k < d_ext_k; ++d_k, ++py)
                  {
                     (*py) = (*px);
                  }
               }
            }
         }
      }
   }

   TTNS_DEBUG("Get_Ham_diag : Finished.");
}

template<size_t N, class Q>
void Get_BijQij_diag (
      const btas::SDArray<1>& op_Bij_diag, size_t op1,
      const btas::SDArray<1>& op_Qij_diag, size_t op2,
            btas::QSDArray<N, Q>& diag)
{
   TTNS_DEBUG("Get_BijQij_diag : Started.");
   TTNS_ASSERT(op1 <= op2, "Get_BijQij_diag: assumed op1 <= op2.");

   size_t str_0000f = diag.stride(op2);
   size_t str_000ff = str_0000f*diag.shape(op2);
   size_t str_00fff = diag.stride(op1);
   size_t str_0ffff = str_00fff*diag.shape(op1);

   size_t ext_i = diag.size()/str_0ffff;
   size_t ext_k = str_00fff/str_000ff;
   size_t ext_m = str_0000f;

   for(auto it1 = op_Bij_diag.begin(); it1 != op_Bij_diag.end(); ++it1)
   {
      size_t idx_0j000 = it1->first*str_00fff;

      for(auto it2 = op_Qij_diag.begin(); it2 != op_Qij_diag.end(); ++it2)
      {
         size_t idx_0j0l0 = idx_0j000 + it2->first*str_0000f;

         for(size_t i = 0; i < ext_i; ++i)
         {
            size_t idx_ij0l0 = i*str_0ffff + idx_0j0l0;

            for(size_t k = 0; k < ext_k; ++k)
            {
               size_t idx_ijkl0 = idx_ij0l0 + k*str_000ff;

               for(size_t m = 0; m < ext_m; ++m)
               {
                  size_t idx_ijklm = idx_ijkl0 + m;

                  auto it3 = diag.reserve(idx_ijklm);

                  if(it3 == diag.end()) continue;

                  btas::DArray<1>& x = *it1->second;
                  btas::DArray<1>& y = *it2->second;
                  btas::DArray<N>& z = *it3->second;

                  size_t d_ext_m = z.stride(op2);
                  size_t d_ext_l = y.size();
                  size_t d_ext_k = z.stride(op1)/d_ext_l/d_ext_m;
                  size_t d_ext_j = x.size();
                  size_t d_ext_i = z.size()/d_ext_j/z.stride(op1);

                  double* pz = z.data();
                  for(size_t d_i = 0; d_i < d_ext_i; ++d_i)
                  {
                     const double* px = x.data();
                     for(size_t d_j = 0; d_j < d_ext_j; ++d_j, ++px)
                     {
                        for(size_t d_k = 0; d_k < d_ext_k; ++d_k)
                        {
                           const double* py = y.data();
                           for(size_t d_l = 0; d_l < d_ext_l; ++d_l, ++py)
                           {
                              double value = (*px) * (*py);
                              for(size_t d_m = 0; d_m < d_ext_m; ++d_m, ++pz)
                              {
                                 (*pz) = value;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   TTNS_DEBUG("Get_BijQij_diag : Finished.");
}

template<size_t N, class Q>
void ContractDiagonals (const btas::TVector<Block*, N>& bs, btas::QSDArray<N, Q>& diag)
{
   TTNS_DEBUG("ContractDiagonals : Started.");

   for(size_t op1 = 0; op1 < N; ++op1)
   {
      if(!bs[op1]) continue;

      btas::SDArray<1> op_Ham_diag;
      DiagonalElement(bs[op1]->Ham(), op_Ham_diag);
      Get_Ham_diag(op_Ham_diag, op1, diag); // diag += op_Ham_diag * 1...1

      auto loop_index = bs[op1]->inside();

      for(size_t op2 = op1+1; op2 < N; ++op2)
      {
         if(!bs[op2]) continue;

         for(size_t i = 0; i < loop_index.size(); ++i)
         {
            size_t ix = loop_index[i];

            if(bs[op2]->CreDesComp(ix, ix).size() == 0) continue;

            btas::SDArray<1> op_Bii_diag;
            DiagonalElement(bs[op1]->CreDes(ix, ix), op_Bii_diag);

            btas::SDArray<1> op_Qii_diag;
            DiagonalElement(bs[op2]->CreDesComp(ix, ix), op_Qii_diag);

            Get_BijQij_diag(op_Bii_diag, op1, op_Qii_diag, op2, diag); // diag += op_Bii_diag * op_Qii_diag * 1...1

            for(size_t j = i+1; j < loop_index.size(); ++j)
            {
               size_t jx = loop_index[j];

               if(bs[op2]->CreDesComp(ix, jx).size() == 0) continue;

               btas::SDArray<1> op_Bij_diag;
               DiagonalElement(bs[op1]->CreDes(ix, jx), op_Bij_diag);

               btas::SDArray<1> op_Qij_diag;
               DiagonalElement(bs[op2]->CreDesComp(ix, jx), op_Qij_diag);

               btas::Scal(2.0, op_Bij_diag); // think of ( Bij * Qij + Bji * Qji )
               Get_BijQij_diag(op_Bij_diag, op1, op_Qij_diag, op2, diag); // diag += op_Bij_diag * op_Qij_diag * 1...1
            }
         }
      }
   }

   TTNS_DEBUG("ContractDiagonals : Finished.");
}

/// compute diagonal Ham. diag is assumed to be resized already
template<size_t N, class Q>
void ComputeDiagonals (const btas::TVector<Block*, N>& bs, btas::QSDArray<N, Q>& diag)
{
   TTNS_DEBUG("ComputeDiagonals(OneDot) : Started.");

   ContractDiagonals(bs, diag);

   TTNS_DEBUG("ComputeDiagonals(OneDot) : Finished.");
}

/// compute diagonal Ham. diag is assumed to be resized already
template<size_t M, size_t N, class Q>
void ComputeDiagonals (const btas::TVector<Block*, M>& bsA, const btas::TVector<Block*, N>& bsB, btas::QSDArray<M+N, Q>& diag)
{
   TTNS_DEBUG("ComputeDiagonals(TwoDot) : Started.");

   btas::TVector<Block*, M+N> bs;
   for(size_t i = 0; i < M; ++i) bs[i]   = bsA[i];
   for(size_t i = 0; i < N; ++i) bs[i+M] = bsB[i];
   ContractDiagonals(bs, diag);

   TTNS_DEBUG("ComputeDiagonals(TwoDot) : Finished.");
}

} // namespace ttns

#endif // __TTNS_COMPUTE_DIAGONALS_H

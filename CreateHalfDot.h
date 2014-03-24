#ifndef __TTNS_CREATE_HALF_DOT_H
#define __TTNS_CREATE_HALF_DOT_H

#include <algorithm>

#include "Renormalize.h"
#include "Randomize.h"

namespace ttns
{

template<size_t N, class Q>
void CreateHalfDot (
      const btas::QSDArray<N, Q>& tns,
      const btas::TVector<Block*, N-2>& bs,
      const size_t& depth,
      const size_t& index,
            btas::QSDArray<3, Q>& mps,
            btas::QSDArray<N-1, Q>& res,
            Block* ops)
{
   TTNS_DEBUG("CreateHalfDot : Started.");

   if(depth == 0) // trunk
   {
      btas::QSDArray<N, Q> tns_scr;

      {
         btas::QSDArray<N, Q> tns_cpy(tns);

         std::vector<int> p1;
         std::vector<int> p2;
         for(size_t i = 0; i < index; ++i)
         {
            p1.push_back(i);
            p2.push_back(index);
         }
         tns_cpy.parity(p1, p2);

         btas::IVector<N> reorder;
         for(size_t i = 0; i < index; ++i) reorder[i] = i;
         for(size_t i = index; i < N-2; ++i) reorder[i] = i+1;
         reorder[N-2] = index;
         reorder[N-1] = N-1;

         btas::Permute(tns_cpy, reorder, tns_scr);
      }

      btas::SDArray<1> sval;
      btas::QSDArray<3, Q> mps_scr;
      btas::QSDArray<N-1, Q> res_scr;

      btas::Gesvd(tns_scr, sval, res_scr, mps_scr, 0);
      btas::Dimm(sval, mps_scr);

      btas::Permute(mps_scr, btas::shape(1,0,2), mps);

      btas::IVector<N-1> reorder_back;
      reorder_back[0] = N-2;
      for(size_t i = 0; i < N-2; ++i) reorder_back[i+1] = i;
      btas::Permute(res_scr, reorder_back, res);

      Blocking(res, res, bs, 0, ops);
   }
   else // branch
   {
      if(index == N-1) // forward
      {
         btas::SDArray<1> sval;

         btas::Gesvd(tns, sval, res, mps, 0);
         btas::Dimm(sval, mps);

         Blocking(res, res, bs, N-2, ops);
      }
      else // backward
      {
         btas::QSDArray<N, Q> tns_scr;

         {
            btas::QSDArray<N, Q> tns_cpy(tns);

            std::vector<int> p1;
            std::vector<int> p2;
            for(int i = 0; i < index; ++i)
            {
               p1.push_back(i);
               p2.push_back(index);
               p1.push_back(i);
               p2.push_back(N-2);
            }
            for(int i = index+1; i < N-2; ++i)
            {
               p1.push_back(i);
               p2.push_back(N-2);
            }
            tns_cpy.parity(p1, p2);

            btas::IVector<N> reorder;
            for(size_t i = 0; i < index; ++i) reorder[i] = i;
            for(size_t i = index; i < N-3; ++i) reorder[i] = i+1;
            reorder[N-3] = N-1;
            reorder[N-2] = index;
            reorder[N-1] = N-2;

            btas::Permute(tns_cpy, reorder, tns_scr);
         }

         btas::SDArray<1> sval;
         btas::QSDArray<3, Q> mps_scr;
         btas::QSDArray<N-1, Q> res_scr;

         btas::Gesvd(tns_scr, sval, res_scr, mps_scr, 0);
         btas::Dimm(sval, mps_scr);

         btas::Permute(mps_scr, btas::shape(1,2,0), mps);

         btas::IVector<N-1> reorder_back;
         reorder_back[0] = N-2;
         for(size_t i = 0; i < N-2; ++i) reorder_back[i+1] = i;
         btas::Permute(res_scr, reorder_back, res);

         Blocking(res, res, bs, 0, ops);
      }
   }

   TTNS_DEBUG("CreateHalfDot : Finished.");
}

template<size_t N, class Q>
void CreateRandHalfDot (
      const btas::QSDArray<N, Q>& tns,
      const btas::TVector<Block*, N-2>& bs,
      const size_t& depth,
      const size_t& index,
      const btas::QSDArray<3, Q>& mps,
      const SweepParameters& param,
            btas::QSDArray<4, Q>& psi,
            btas::QSDArray<N-1, Q>& res,
            Block* ops)
{
   TTNS_DEBUG("CreateRandHalfDot : Started.");

   btas::Dshapes bonds = tns.dshape(index);
   size_t D = std::accumulate(bonds.begin(), bonds.end(), 0);
   D = std::max(D, param.M());

   std::cout << "\t\t\tRandomized Half-Block: Keeping 4 x " << D << " states " << std::endl;
   D *= 4;

   if(depth == 0) // trunk
   {
      btas::QSDArray<N+1, Q> psi_full_scr;

      {
         btas::QSDArray<N, Q> tns_cpy(tns);

         std::vector<int> p1;
         std::vector<int> p2;
         for(size_t i = 0; i < index; ++i)
         {
            p1.push_back(i);
            p2.push_back(index);
         }
         tns_cpy.parity(p1, p2);

         btas::QSDArray<N+1, Q> psi_full;
         btas::Contract(1.0, mps, btas::shape(2), tns_cpy, btas::shape(index), 1.0, psi_full);

         btas::TVector<btas::Dshapes, N+1> psi_dshape = psi_full.dshape();
         psi_dshape[1] = btas::Dshapes(4, 1);
         psi_dshape[N] = btas::Dshapes(4, 1);
         Randomize(param.Noise(), psi_full, psi_dshape);

         btas::IVector<N+1> reorder;
         for(size_t i = 0; i < N-2; ++i) reorder[i] = i+2;
         reorder[N-2] = 0;
         reorder[N-1] = 1;
         reorder[N] = N;

         btas::Permute(psi_full, reorder, psi_full_scr);
      }

      btas::SDArray<1> sval;
      btas::QSDArray<4, Q> psi_scr;
      btas::QSDArray<N-1, Q> res_scr;

      btas::Gesvd(psi_full_scr, sval, res_scr, psi_scr, D);
      btas::Dimm(sval, psi_scr);

      btas::Permute(psi_scr, btas::shape(1,2,0,3), psi);

      btas::IVector<N-1> reorder_back;
      reorder_back[0] = N-2;
      for(size_t i = 0; i < N-2; ++i) reorder_back[i+1] = i;
      btas::Permute(res_scr, reorder_back, res);

      Blocking(res, res, bs, 0, ops);
   }
   else // branch
   {
      if(index == N-1) // backward
      {
         btas::QSDArray<N+1, Q> psi_full;
         btas::Contract(1.0, tns, btas::shape(index), mps, btas::shape(0), 1.0, psi_full);

         btas::TVector<btas::Dshapes, N+1> psi_dshape = psi_full.dshape();
         psi_dshape[N-2] = btas::Dshapes(4, 1);
         if(depth == 1)
            psi_dshape[N] = btas::Dshapes(4, 1);
         else
            psi_dshape[N-1] = btas::Dshapes(4, 1);
         Randomize(param.Noise(), psi_full, psi_dshape);

         btas::SDArray<1> sval;

         btas::Gesvd(psi_full, sval, res, psi, D);
         btas::Dimm(sval, psi);

         Blocking(res, res, bs, N-2, ops);
      }
      else // forward
      {
         btas::QSDArray<N+1, Q> psi_full_scr;

         {
            btas::QSDArray<N, Q> tns_cpy(tns);

            std::vector< int > p1;
            std::vector< int > p2;
            for(int i = 0; i < index; ++i)
            {
               p1.push_back(i);
               p2.push_back(index);
               p1.push_back(i);
               p2.push_back(N-2);
            }
            for(int i = index + 1; i < N - 2; ++i)
            {
               p1.push_back(i);
               p2.push_back(N-2);
            }
            tns_cpy.parity(p1, p2);

            btas::QSDArray<N+1, Q> psi_full;
            btas::Contract(1.0, mps, btas::shape(2), tns_cpy, btas::shape(index), 1.0, psi_full);

            btas::TVector<btas::Dshapes, N+1> psi_dshape = psi_full.dshape();
            psi_dshape[1] = btas::Dshapes(4, 1);
            psi_dshape[N-1] = btas::Dshapes(4, 1);
            Randomize(param.Noise(), psi_full, psi_dshape);

            btas::IVector<N+1> reorder;
            for(size_t i = 0; i < N-3; ++i) reorder[i] = i+2;
            reorder[N-3] = N;
            reorder[N-2] = 0;
            reorder[N-1] = 1;
            reorder[N] = N-1;

            btas::Permute(psi_full, reorder, psi_full_scr);
         }

         btas::SDArray<1> sval;
         btas::QSDArray<4, Q> psi_scr;
         btas::QSDArray<N-1, Q> res_scr;

         btas::Gesvd(psi_full_scr, sval, res_scr, psi_scr, D);
         btas::Dimm(sval, psi_scr);

         btas::Permute(psi_scr, btas::shape(1,2,3,0), psi);

         btas::IVector<N-1> reorder_back;
         reorder_back[0] = N-2;
         for(size_t i = 0; i < N-2; ++i) reorder_back[i+1] = i;
         btas::Permute(res_scr, reorder_back, res);

         Blocking(res, res, bs, 0, ops);
      }
   }

   TTNS_DEBUG("CreateRandHalfDot : Finished.");
}

template<size_t N, class Q>
void BackTransform (
      const btas::QSDArray<3, Q>& mps,
      const btas::QSDArray<N-1, Q>& res,
      const size_t& depth,
      const size_t& index,
            btas::QSDArray<N, Q>& tns)
{
   TTNS_DEBUG("BackTransform : Started.");

   if(depth == 0)
   {
      btas::QSDArray<N, Q> tns_scr;
      btas::Contract(1.0, mps, btas::shape(1), res, btas::shape(0), 1.0, tns_scr);

      btas::IVector<N> reorder;
      for(size_t i = 0; i < index; ++i) reorder[i] = i+2;
      for(size_t i = index+1; i < N-1; ++i) reorder[i] = i+1;
      reorder[index] = 0;
      reorder[N-1] = 1;
      btas::Permute(tns_scr, reorder, tns);

      std::vector<int> p1;
      std::vector<int> p2;
      for(size_t i = 0; i < index; ++i)
      {
         p1.push_back(i);
         p2.push_back(index);
      }
      tns.parity(p1, p2);
   }
   else
   {
      if(index == N-1)
      {
         tns.clear();
         btas::Contract(1.0, res, btas::shape(N-2), mps, btas::shape(0), 1.0, tns);
      }
      else
      {
         btas::QSDArray<N, Q> tns_scr;
         btas::Contract(1.0, mps, btas::shape(2), res, btas::shape(0), 1.0, tns_scr);

         btas::IVector<N> reorder;
         for(size_t i = 0; i < index; ++i) reorder[i] = i+2;
         for(size_t i = index+1; i < N-2; ++i) reorder[i] = i+1;
         reorder[index] = 0;
         reorder[N-2] = 1;
         reorder[N-1] = N-1;
         btas::Permute(tns_scr, reorder, tns);

         std::vector<int> p1;
         std::vector<int> p2;
         for(size_t i = 0; i < index; ++i)
         {
            p1.push_back(i);
            p2.push_back(index);
            p1.push_back(i);
            p2.push_back(N-2);
         }
         for(size_t i = index+1; i < N-2; ++i)
         {
            p1.push_back(i);
            p2.push_back(N-2);
         }
         tns.parity(p1, p2);
      }
   }

   TTNS_DEBUG("BackTransform : Finished.");
}

} // namespace ttns

#endif // __TTNS_CREATE_HALF_DOT_H

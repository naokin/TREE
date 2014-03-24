#ifndef __TTNS_COMPUTE_SIGMA_VECTOR_H
#define __TTNS_COMPUTE_SIGMA_VECTOR_H

#include "Block.h"
#include "ComputeHam.h"

namespace ttns
{

/// Compute sigma vector: V = H * C
template<size_t N, class Q>
void ComputeSigmaVector (
      const btas::TVector<Block*, N>& bs,
      const btas::IVector<N>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::QSDArray<N, Q>& sgv)
{
   // sorting blocks by block-size
   btas::IVector<N> st_order;
   btas::IVector<N> reorder;
   btas::IVector<N> reorder_back;

   {
      std::multimap<size_t, size_t> block_size_map;

      for(size_t i = 0; i < N; ++i)
      {
         if(bs[i])
            block_size_map.insert(std::make_pair(bs[i]->size(), i));
         else
            block_size_map.insert(std::make_pair(0, i));
      }

      auto it = block_size_map.begin();

      for(size_t i = 0; i < N; ++i, ++it)
      {
         reorder[i] = it->second;
         reorder_back[it->second] = i;
      }

      for(size_t i = 0; i < N; ++i)
      {
         st_order[i] = order[reorder[i]];
      }
   }

   // permute psi with reorder
   btas::QSDArray<N, Q> psi_scr;
   btas::Permute(psi, reorder, psi_scr);

   // scratch storage of sigma vector
   btas::QSDArray<N, Q> sgv_scr(psi_scr.q(), psi_scr.qshape());

   // permute blocks
   btas::TVector<Block*, N> st_bs = btas::permute(bs, reorder);

   // compute Ci x TPS stored as intermediate
// btas::TArray<btas::QSDArray<N, Q>, 1> CreStore;
// btas::TArray<btas::QSDArray<N, Q>, 1> DesStore;
// ComputeCreStorage(st_bs, psi_scr, CreStore, DesStore);

// std::cout << "psi.size() = " << psi_scr.size() << std::endl;

   // computing sigma vector
   ComputeHamHam(st_bs, st_order, psi_scr, sgv_scr);

   ComputeHamCiRi(st_bs, st_order, psi_scr, sgv_scr);

   ComputeHamCijRij(st_bs, st_order, psi_scr, sgv_scr);

   ComputeHamCiCjRij(st_bs, st_order, psi_scr, sgv_scr);

   ComputeHamCiCjDkDl(st_bs, st_order, psi_scr, sgv_scr);

   // reorder back to original order
   btas::Permute(sgv_scr, reorder_back, sgv);
}

/// Compute Cre * TNS as intermediate
template<size_t N, class Q>
void ComputeCreStorage (
      const btas::TVector<Block*, N>& bs,
      const btas::IVector<N>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::TArray<btas::QSDArray<N, Q>, 1>& CreStore,
            btas::TArray<btas::QSDArray<N, Q>, 1>& DesStore)
{
   CreStore.resize(Block::orbitals());
   DesStore.resize(Block::orbitals());

   for(size_t op1 = 0; op1 < N; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(size_t i = 0; i < op1_loop_index.size(); ++i)
      {
         size_t ix = op1_loop_index[i];

         btas::Contract(1.0, bs[op1]->Cre(ix), btas::shape(0), psi, btas::shape(op1), 1.0, CreStore(ix));
         btas::Contract(1.0, bs[op1]->Cre(ix), btas::shape(1), psi, btas::shape(op1), 1.0, DesStore(ix));
      }
   }
}

} // namespace ttns

#endif // __TTNS_COMPUTE_SIGMA_VECTOR_H

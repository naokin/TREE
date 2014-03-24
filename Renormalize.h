#ifndef __TTNS_RENORMALIZE_H
#define __TTNS_RENORMALIZE_H

#include "Block.h"
#include "ComputeHam.h"
#include "OpComponentsNormal.h"
#include "OpComponentsComp.h"

namespace ttns
{

// renormalize superblock
template<class Q>
inline void DoRenormalize (
      const DmrgOperator& superBlockOps,
      const btas::QSDArray<2, Q>& rotBra,
      const btas::QSDArray<2, Q>& rotKet,
            DmrgOperator& newBlockOps)
{
   DmrgOperator opsScr;
   btas::Contract( 1.0, rotBra, btas::shape(1), superBlockOps, btas::shape(0), 1.0, opsScr);
   btas::Contract( 1.0, opsScr, btas::shape(1), rotKet, btas::shape(1), 1.0, newBlockOps);
}

template<class Q>
void Renormalize (
      const Block* superBlock,
      const btas::QSDArray<2, Q>& rotBra,
      const btas::QSDArray<2, Q>& rotKet,
            Block* newBlock)
{
   std::vector<size_t> loop_index(superBlock->inside());
   std::vector<size_t> comp_index(superBlock->outside());
   newBlock->resize(loop_index);

   DoRenormalize(superBlock->Ham(), rotBra, rotKet, newBlock->Ham());

   for(size_t i = 0; i < loop_index.size(); ++i)
   {
      size_t ix = loop_index[i];

      DoRenormalize(superBlock->Cre(ix), rotBra, rotKet, newBlock->Cre(ix));
      DoRenormalize(superBlock->CreDes(ix, ix), rotBra, rotKet, newBlock->CreDes(ix, ix));

      for(size_t j = i+1; j < loop_index.size(); ++j)
      {
         size_t jx = loop_index[j];

         DoRenormalize(superBlock->CreCre(ix, jx), rotBra, rotKet, newBlock->CreCre(ix, jx));
         ScaledCopy(-1.0, newBlock->CreCre(ix, jx), newBlock->CreCre(jx, ix));

         DoRenormalize(superBlock->CreDes(ix, jx), rotBra, rotKet, newBlock->CreDes(ix, jx));
         btas::Permute(newBlock->CreDes(ix, jx), btas::shape(1, 0), newBlock->CreDes(jx, ix));
      }
   }

   for(size_t i = 0; i < comp_index.size(); ++i)
   {
      size_t ix = comp_index[i];

      DoRenormalize(superBlock->CreComp(ix), rotBra, rotKet, newBlock->CreComp(ix));
      DoRenormalize(superBlock->CreDesComp(ix, ix), rotBra, rotKet, newBlock->CreDesComp(ix, ix));

      for(size_t j = i + 1; j < comp_index.size(); ++j)
      {
         size_t jx = comp_index[j];

         DoRenormalize(superBlock->CreCreComp(ix, jx), rotBra, rotKet, newBlock->CreCreComp(ix, jx));
         ScaledCopy(-1.0, newBlock->CreCreComp(ix, jx), newBlock->CreCreComp(jx, ix));

         DoRenormalize(superBlock->CreDesComp(ix, jx), rotBra, rotKet, newBlock->CreDesComp(ix, jx));
         btas::Permute(newBlock->CreDesComp(ix, jx), btas::shape(1, 0), newBlock->CreDesComp(jx, ix));
      }
   }
}

// blocking and decimation 
template<size_t K, class Q>
void Blocking (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const size_t& indexCurrent,
            Block* newBlock)
{
   TTNS_DEBUG("Blocking : Started.");

   // resizing newBlock
   std::vector<size_t> new_loop_index;
   for(size_t i = 0; i < K; ++i)
   {
      if(!bs[i]) continue;

      std::vector<size_t> opi_loop_index(bs[i]->inside());
      new_loop_index.insert(new_loop_index.end(), opi_loop_index.begin(), opi_loop_index.end());
   }
   newBlock->resize(new_loop_index);

   // preconditioning
   btas::QSDArray<K+1, Q> braCopy(bra);
   btas::QSDArray<K+1, Q> ketCopy(ket);
   if(indexCurrent < K)
   {
      std::vector<int> p1;
      std::vector<int> p2;
      for(size_t i = 0; i < indexCurrent; ++i)
      {
         p1.push_back(i);
         p2.push_back(indexCurrent);
      }
      braCopy.parity(p1, p2);
      ketCopy.parity(p1, p2);
   }

   btas::IVector<K> reorderBlocks;
   btas::IVector<K+1> reorderTns;

   reorderTns[K] = indexCurrent;
   {
      std::multimap<size_t, size_t> sortingMap;
      for(size_t i = 0; i < K; ++i)
      {
         if(bs[i])
            sortingMap.insert(std::make_pair(bs[i]->size(), i));
         else
            sortingMap.insert(std::make_pair(0, i));
      }

      auto it = sortingMap.begin();
      for(size_t i = 0; i < K; ++i, ++it)
      {
         size_t ib = it->second;
         reorderBlocks[i] = ib;

         if(ib < indexCurrent)
            reorderTns[i] = ib;
         else
            reorderTns[i] = ib+1;
      }
   }

   // sorting tps indices
   btas::QSDArray<K+1, Q> braSorted;
   btas::Permute(braCopy, reorderTns, braSorted);
   btas::QSDArray<K+1, Q> ketSorted;
   btas::Permute(ketCopy, reorderTns, ketSorted);

   // sorting bs
   btas::TVector<Block*, K> bsSorted;
   for(size_t i = 0; i < K; ++i) bsSorted[i] = bs[reorderBlocks[i]];

   // renormalization
   DoBlocking(braSorted, ketSorted, bsSorted, reorderTns, newBlock);

   TTNS_DEBUG("Blocking : Finished.");
}

//
// R: renormalized, T: trunk-side, B: branch-side
//
// branch site;
//              RB RB RB   T
// ordering: [( 0, 1, 2 ), 3 ] = renormalized to trunk (forward)
//              RT RT RT   B
// ordering: [( 1, 2, 3 ), 0 ] = renormalized to branch 0 (backward)
// ordering: [( 0, 2, 3 ), 1 ] = renormalized to branch 1 (backward)
//
// center site;
//              RT RT RT   B
// ordering: [( 1, 2, 3 ), 0 ] = renormalized to branch 0 (backward)
// ordering: [( 0, 2, 3 ), 1 ] = renormalized to branch 0 (backward)
// ordering: [( 0, 1, 3 ), 2 ] = renormalized to branch 0 (backward)

// normal renormalization
template<size_t K, class Q>
void DoBlocking (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K+1>& ordering,
            Block* newBlock)
{
   TTNS_DEBUG("DoBlocking : Started.");

   ComputeHamBlock      (bra, ket, bs, ordering, newBlock);

   ComputeNormalBlockCi (bra, ket, bs, ordering, newBlock);
   ComputeNormalBlockCij(bra, ket, bs, ordering, newBlock);

   ComputeCompBlockRi   (bra, ket, bs, ordering, newBlock);
   ComputeCompBlockRij  (bra, ket, bs, ordering, newBlock);

   TTNS_DEBUG("DoBlocking : Finished.");
}

// compute local ham of newBlock
template<size_t K, class Q>
void ComputeHamBlock (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K+1>& ordering,
            Block* newBlock)
{
   TTNS_DEBUG("ComputeHamBlock : Started.");

   btas::IVector<K> blockOrdering;
   for(size_t i = 0; i < K; ++i) blockOrdering[i] = ordering[i];

   btas::QSDArray<K+1, Q> sgv(bra.q(), bra.qshape());

   TTNS_DEBUG("ComputeHamBlock : Calling ComputeHamHam.");
   ComputeHamHam     (bs, blockOrdering, bra, sgv);
   TTNS_DEBUG("ComputeHamBlock : Calling ComputeHamCiRi.");
   ComputeHamCiRi    (bs, blockOrdering, bra, sgv);
   TTNS_DEBUG("ComputeHamBlock : Calling ComputeHamCijRij.");
   ComputeHamCijRij  (bs, blockOrdering, bra, sgv);
   TTNS_DEBUG("ComputeHamBlock : Calling ComputeHamCiCjRij.");
   ComputeHamCiCjRij (bs, blockOrdering, bra, sgv);
   TTNS_DEBUG("ComputeHamBlock : Calling ComputeHamCiCjDkDl.");
   ComputeHamCiCjDkDl(bs, blockOrdering, bra, sgv);

   btas::IVector<K> contracts;
   for(size_t i = 0; i < K; ++i) contracts[i] = i;

   btas::Contract( 1.0, sgv, contracts, ket, contracts, 1.0, newBlock->Ham());

   TTNS_DEBUG("ComputeHamBlock : Finished.");
}

template<size_t K, class Q>
void ComputeNormalBlockCi (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K+1>& ordering,
            Block* newBlock)
{
   for(long op1 = 0; op1 < K; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      CiBlockCi(bra, ket, bs[op1], op1, op1_loop_index, ordering, newBlock);
   }
}

template<size_t K, class Q>
void ComputeNormalBlockCij (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K+1>& ordering,
            Block* newBlock)
{
   long N_small_blocks = 0;
   for(size_t i = 0; i < K - 1; ++i)
   {
      if(bs[i] && bs[i]->size() > 2) break; // theoretical
      ++N_small_blocks;
   }

   for(long op1 = 0; op1 < N_small_blocks; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      CijBlockCij(bra, ket, bs[op1], op1, op1_loop_index, ordering, newBlock);

      for(long op2 = op1 + 1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

//       CijBlockCiCjS(bra, ket, bs[op1], op1, op1_loop_index,
//       CijBlockCiCjL(bra, ket, bs[op1], op1, op1_loop_index,
         CijBlockCiCjM(bra, ket, bs[op1], op1, op1_loop_index,
                                 bs[op2], op2, op2_loop_index, ordering, newBlock);
      }
   }

   for(long op1 = N_small_blocks; op1 < K; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      CijBlockCij(bra, ket, bs[op1], op1, op1_loop_index, ordering, newBlock);

      for(long op2 = op1 + 1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

//       CijBlockCiCjL(bra, ket, bs[op1], op1, op1_loop_index,
         CijBlockCiCjM(bra, ket, bs[op1], op1, op1_loop_index,
                                 bs[op2], op2, op2_loop_index, ordering, newBlock);
      }
   }
}

template<size_t K, class Q>
void ComputeCompBlockRi (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K+1>& ordering,
            Block* newBlock)
{
   TTNS_DEBUG("ComputeCompBlockRi : Started.");

   long N_small_blocks = 0;
   for(size_t i = 0; i < K - 1; ++i)
   {
      if(bs[i] && bs[i]->size() > 2) break; // theoretical
      ++N_small_blocks;
   }

   std::vector<size_t> new_comp_index(newBlock->outside());
   btas::TArray<btas::SDArray<K+1>, 1> CreCompRi;

   TTNS_DEBUG("ComputeCompBlockRi : Small(1) * Large(2)");

   for(long op1 = 0; op1 < N_small_blocks; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      TTNS_DEBUG("ComputeCompBlockRi : Ri -> Ri [" << op1 << "]");

      // Ri x 1
      RiBlockRi(bra, bs[op1], op1, ordering, CreCompRi, new_comp_index);

      for(long op2 = op1 + 1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         TTNS_DEBUG("ComputeCompBlockRi : CjRij -> Ri [" << op1 << "," << op2 << "]");

         RiBlockCjRijS(bra, bs[op1], op1, op1_loop_index,
                            bs[op2], op2, ordering, CreCompRi, new_comp_index);

         RiBlockCjRijS(bra, bs[op2], op2, op2_loop_index,
                            bs[op1], op1, ordering, CreCompRi, new_comp_index);

         for(long op3 = op2 + 1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;
            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            TTNS_DEBUG("ComputeCompBlockRi : CjCkDl -> Ri [" << op1 << "," << op2 << "," << op3 << "]");

            // Cj x Ck x Dl
            RiBlockCjCkDlL(bra, bs[op1], op1, op1_loop_index,
//          RiBlockCjCkDlS(bra, bs[op1], op1, op1_loop_index,
                                bs[op2], op2, op2_loop_index,
                                bs[op3], op3, op3_loop_index, ordering, CreCompRi, new_comp_index);
         }
      }
   }

   TTNS_DEBUG("ComputeCompBlockRi : Large(1) * Large(2)");

   for(long op1 = N_small_blocks; op1 < K; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      TTNS_DEBUG("ComputeCompBlockRi : Ri -> Ri [" << op1 << "]");

      // Ri x 1
      RiBlockRi(bra, bs[op1], op1, ordering, CreCompRi, new_comp_index);

      for(long op2 = op1 + 1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         TTNS_DEBUG("ComputeCompBlockRi : CjRij -> Ri [" << op1 << "," << op2 << "]");

         RiBlockCjRijL(bra, bs[op1], op1, op1_loop_index,
                            bs[op2], op2, ordering, CreCompRi, new_comp_index);

         RiBlockCjRijL(bra, bs[op2], op2, op2_loop_index,
                            bs[op1], op1, ordering, CreCompRi, new_comp_index);

         for(long op3 = op2 + 1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;
            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            TTNS_DEBUG("ComputeCompBlockRi : CjCkDl -> Ri [" << op1 << "," << op2 << "," << op3 << "]");

            // Cj x Ck x Dl
            RiBlockCjCkDlL(bra, bs[op1], op1, op1_loop_index,
                                bs[op2], op2, op2_loop_index,
                                bs[op3], op3, op3_loop_index,
                                ordering, CreCompRi, new_comp_index);
         }
      }
   }

   TTNS_DEBUG("ComputeCompBlockRi : Intermediate Ri * Ket");

   if(CreCompRi.size() > 0)
   {
      btas::TVector< size_t, K > contracts;
      for(size_t i = 0; i < K; ++i) contracts[i] = i;

      for(size_t i = 0; i < new_comp_index.size(); ++i)
      {
         size_t ix = new_comp_index[i];
         if(CreCompRi(i).size() > 0)
            btas::Contract( 1.0, CreCompRi(i), contracts, ket, contracts,
                            1.0, newBlock->CreComp(ix));
      }
   }

   TTNS_DEBUG("ComputeCompBlockRi : Finished.");
}

template<size_t K, class Q>
void ComputeCompBlockRij (
      const btas::QSDArray<K+1, Q>& bra,
      const btas::QSDArray<K+1, Q>& ket,
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K+1>& ordering,
            Block* newBlock)
{
    std::vector<size_t> new_comp_index(newBlock->outside());
    for(long op1 = 0; op1 < K; ++op1)
    {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      // Rij x 1
      RijBlockRij(bra, ket, bs[op1], op1, ordering, newBlock, new_comp_index);
//    RijBlockCkl(op1_loop_index, newBlock, new_comp_index);

      for(long op2 = op1 + 1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         // Ck(a) x Cl(b)
         RijBlockCkCl(op1_loop_index, op2_loop_index, newBlock, new_comp_index);
      }
   }
}

} // namespace ttns

#endif // __TTNS_RENORMALIZE_H

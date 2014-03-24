#ifndef __TTNS_OP_COMPONENTS_COMP_H
#define __TTNS_OP_COMPONENTS_COMP_H

#include "OpComponentsUtil.h"
#include "ttns_assert.h"

namespace ttns
{

#ifndef THRE_INT
#define THRE_INT 1.0e-16
#endif

//==========================================================================================//
// Compute Complementary Blocks for renormalization.                                        //
//==========================================================================================//

//
// CreComp
//

//
// Ri := Ri x 1
//
template<size_t N, class Q>
void RiBlockRi (
      const btas::QSDArray<N, Q>& bra,
            Block* blockRi, const size_t& indexRi,
      const btas::IVector<N>& order,
            btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
      const std::vector<size_t>& i_loop_index)
{
   TTNS_DEBUG("RiBlockRi : Started.");

   const size_t K = N-1;

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra); 

   std::vector<int> p1;
   for(size_t i = 0; i < K; ++i)
      if(order[indexRi] < order[i]) p1.push_back(i);
   braCopy.parity(p1);

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexRi);

   size_t ni = i_loop_index.size();
   if(CreCompRi.size() == 0) CreCompRi.resize(ni);

   for(size_t i = 0; i < ni; ++i)
   {
      size_t ix = i_loop_index[i];
      if(blockRi->CreComp(ix).size() == 0) continue;

//    TTNS_DEBUG("RiBlockRi : Ri * Bra -> Ri [" << ix << "]");
//    TTNS_DEBUG("RiBlockRi : blockRi->CreComp(ix) : " << blockRi->CreComp(ix));
//    TTNS_DEBUG("RiBlockRi : braCopy : " << braCopy);

      // Ri
      btas::SDArray<N> CreCompScr;
      btas::Contract(1.0, blockRi->CreComp(ix), btas::shape(0), braCopy, btas::shape(indexRi), 1.0, CreCompScr);
      PermutedAxpy(1.0, CreCompScr, reorder, CreCompRi(i));
   }

   TTNS_DEBUG("RiBlockRi : Finished.");
}

//
// Ri := Cj x Rij
//
template<size_t N, class Q>
void RiBlockCjRijS (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<N>& order,
            btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
      const std::vector<size_t>& i_loop_index)
{
   const size_t K = N-1;

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra); 
   std::vector<int> p1;
   for(size_t i = 0; i < K; ++i)
      if(order[indexCj] < order[i]) p1.push_back(i);
   braCopy.parity(p1);
  
   btas::IVector<N> reorder;
   IndexPermute(reorder, indexRij, indexCj);

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();

   if(CreCompRi.size() == 0) CreCompRi.resize(ni);

   for(size_t i = 0; i < ni; ++i)
   {
      size_t ix = i_loop_index[i];

      btas::SDArray<4> SuperBlockPij;
      btas::SDArray<4> SuperBlockQij;

      for(size_t j = 0; j < nj; ++j)
      {
         size_t jx = j_loop_index[j];
         // R(i) = Cj x Pij + Dj x Qij
         if(blockRij->CreCreComp(ix, jx).size() > 0)
            btas::Ger( 1.0, blockCj->Cre(jx), blockRij->CreCreComp(ix, jx), SuperBlockPij);

         if(blockRij->CreDesComp(ix, jx).size() > 0)
            btas::Ger( 1.0, blockCj->Cre(jx), blockRij->CreDesComp(ix, jx), SuperBlockQij);
      }

      btas::SDArray<4> CreCompSuperBlock;

      if(SuperBlockPij.size() > 0)
         PermutedAxpy(-1.0, SuperBlockPij, btas::shape(1, 0, 2, 3), CreCompSuperBlock);

      if(SuperBlockQij.size() > 0)
         PermutedAxpy( 1.0, SuperBlockQij, btas::shape(0, 1, 3, 2), CreCompSuperBlock);

      if(CreCompSuperBlock.size() == 0) continue;

      btas::SDArray<N> CreCompScr;
//    btas::Contract(-1.0, CreCompSuperBlock, btas::shape(0, 2),  /*   Di x CreCompRi */
      btas::Contract( 1.0, CreCompSuperBlock, btas::shape(0, 2),  /* - Di x CreCompRi */
                           braCopy, btas::shape(indexCj, indexRij), 1.0, CreCompScr);
      PermutedAxpy( 1.0, CreCompScr, reorder, CreCompRi(i));
   }
}

// computing CreComp only for renormalization
template<size_t N, class Q>
void RiBlockCjRijL (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<N>& order,
            btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
      const std::vector<size_t>& i_loop_index)
{
   const size_t K = N-1;

   btas::QSDArray<N, Q> braCopy(bra); 

   // parity operation
   std::vector<int> p1;
   for(size_t i = 0; i < K; ++i)
      if(order[indexCj] < order[i]) p1.push_back(i);
   braCopy.parity(p1);
  
   size_t iOffRij = 0; if(indexCj > indexRij) iOffRij++;

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCj, indexRij);

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();
   if(CreCompRi.size() == 0) CreCompRi.resize(ni);

   btas::TArray<btas::SDArray<N>, 1> CreCompScr(ni);

   // R(i)'= Pij'x [ Dj x bra ] + Qij'x [ Cj x bra ]
   for(size_t j = 0; j < nj; ++j)
   {
      size_t jx = j_loop_index[j];
      // Cj x bra
      btas::SDArray<N> CreScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0), braCopy, btas::shape(indexCj), 1.0, CreScr);
      // Dj x bra
      btas::SDArray<N> DesScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(1), braCopy, btas::shape(indexCj), 1.0, DesScr);

      for(size_t i = 0; i < ni; ++i)
      {
         size_t ix = i_loop_index[i];
         if(blockRij->CreCreComp(ix, jx).size() > 0)
            // Ri':= Pij'x Dj
            btas::Contract(-1.0, blockRij->CreCreComp(ix, jx), btas::shape(0), DesScr, btas::shape(indexRij+iOffRij), 1.0, CreCompScr(i));

         if(blockRij->CreDesComp(ix, jx).size() > 0)
            // Ri':= Qij'x Cj
            btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(1), CreScr, btas::shape(indexRij+iOffRij), 1.0, CreCompScr(i));
      }
   }

   for(size_t i = 0; i < ni; ++i)
   {
      if(CreCompScr(i).size() > 0)
//       PermutedAxpy(-1.0, CreCompScr(i), reorder, CreCompRi(i)); // + Di x CreCompRi
         PermutedAxpy( 1.0, CreCompScr(i), reorder, CreCompRi(i)); // - Di x CreCompRi
   }
}

template<size_t N, class Q>
void RiBlockCjCkDlL (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockCk,
      const size_t& indexCk,
      const std::vector<size_t>& k_loop_index,
            Block* blockDl,
      const size_t& indexDl,
      const std::vector<size_t>& l_loop_index,
      const btas::IVector<N>& order,
            btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
      const std::vector<size_t>& i_loop_index)
{
   const size_t K = N - 1;

   // must be indexCj < indexDk < indexDl
   {
      std::vector<size_t> originalOrder;
      originalOrder.push_back(indexCj);
      originalOrder.push_back(indexCk);
      originalOrder.push_back(indexDl);

      std::vector<size_t> sortedOrder(originalOrder);
      std::sort(sortedOrder.begin(), sortedOrder.end());

      if(!std::equal(sortedOrder.begin(), sortedOrder.end(), originalOrder.begin()))
        TTNS_ASSERT(false, "ttns::CreCompCjDkDl; contraction order must be sorted before calling");
   }

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   std::vector<int> p1;
   {
      std::vector<size_t> bounds;
      bounds.push_back(order[indexCj]);
      bounds.push_back(order[indexCk]);
      bounds.push_back(order[indexDl]);
      std::sort(bounds.begin(), bounds.end());

      for(size_t i = 0; i < K; ++i)
      {
         if(bounds[0] < order[i] && order[i] <= bounds[1])
            p1.push_back(i);
         if(bounds[2] < order[i])
            p1.push_back(i);
      }
   }
   braCopy.parity(p1);

   double opSign = 1.0;
   if(order[indexCj] > order[indexCk]) opSign = -opSign;
   if(order[indexCj] > order[indexDl]) opSign = -opSign;
   if(order[indexCk] > order[indexDl]) opSign = -opSign;
  
   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCk, indexDl, indexCj);

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();
   size_t kBlockSize = k_loop_index.size();
   size_t lBlockSize = l_loop_index.size();
   if(CreCompRi.size() == 0) CreCompRi.resize(ni);

   btas::TArray<btas::SDArray<N>, 2> ijCreCreComp(ni, nj);
   btas::TArray<btas::SDArray<N>, 2> ijDesCreComp(ni, nj);

   for(size_t k = 0; k < kBlockSize; ++k)
   {
      size_t kx = k_loop_index[k];
      int kSpin = 1; if(kx % 2 != 0) kSpin = -1;
      btas::SDArray<N> kCreScr;
      btas::Contract( 1.0, blockCk->Cre(kx), btas::shape(0), braCopy, btas::shape(indexCk), 1.0, kCreScr);

      btas::SDArray<N> kDesScr;
      btas::Contract( 1.0, blockCk->Cre(kx), btas::shape(1), braCopy, btas::shape(indexCk), 1.0, kDesScr);

      for(size_t l = 0; l < lBlockSize; ++l)
      {
         size_t lx = l_loop_index[l];
         int lSpin = 1; if(lx % 2 != 0) lSpin = -1;
         // NOTE: check sign problems
         btas::SDArray<N> klCreCreScr;
         btas::Contract(-1.0, blockDl->Cre(lx), btas::shape(0), kCreScr, btas::shape(indexDl), 1.0, klCreCreScr);

         btas::SDArray<N> klCreDesScr;
         btas::Contract(-1.0, blockDl->Cre(lx), btas::shape(1), kCreScr, btas::shape(indexDl), 1.0, klCreDesScr);

         btas::SDArray<N> klDesCreScr;
         btas::Contract(-1.0, blockDl->Cre(lx), btas::shape(0), kDesScr, btas::shape(indexDl), 1.0, klDesCreScr);

         for(size_t i = 0; i < ni; ++i)
         {
            size_t ix = i_loop_index[i];
            int iSpin = 1; if(ix % 2 != 0) iSpin = -1;
            for(size_t j = 0; j < nj; ++j)
            {
               size_t jx = j_loop_index[j];
               int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

               double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
               if(fabs(vikjl) >= THRE_INT && (iSpin == jSpin) && (kSpin == lSpin))
               {
                  btas::Axpy(+opSign*vikjl, klDesCreScr, ijDesCreComp(i, j));
                  btas::Axpy(-opSign*vikjl, klCreDesScr, ijDesCreComp(i, j));
               }

               double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
               if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
               {
                  btas::Axpy(-opSign*vijkl, klCreCreScr, ijCreCreComp(i, j));
                  btas::Axpy(+opSign*vijkl, klCreDesScr, ijDesCreComp(i, j));
               }

               double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
               if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
               {
                  btas::Axpy(+opSign*vijlk, klCreCreScr, ijCreCreComp(i, j));
                  btas::Axpy(-opSign*vijlk, klDesCreScr, ijDesCreComp(i, j));
               }
            }
         }
      }
   }

   for(size_t i = 0; i < ni; ++i)
   {
      size_t ix = i_loop_index[i];
      btas::SDArray<N> iCreCompScr;
      for(size_t j = 0; j < nj; ++j)
      {
         size_t jx = j_loop_index[j];
         // CiComp = Cj x DiCjComp + Dj x CiCjComp
         if(ijDesCreComp(i, j).size() > 0)
            btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0), ijDesCreComp(i, j), btas::shape(indexCj+2), 1.0, iCreCompScr);

         if(ijCreCreComp(i, j).size() > 0)
            btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(1), ijCreCreComp(i, j), btas::shape(indexCj+2), 1.0, iCreCompScr);
      }

      if(iCreCompScr.size() > 0)
         PermutedAxpy( 1.0, iCreCompScr, reorder, CreCompRi(i));
   }
}

// CreCreComp Block
template<size_t N, class Q>
void RijBlockRij (
      const btas::QSDArray<N, Q>& bra,
      const btas::QSDArray<N, Q>& ket,
            Block* blockRij, const size_t& indexRij,
      const btas::IVector<N>& order,
            Block* newBlock,
      const std::vector<size_t>& ij_loop_index)
{
   const size_t K = N - 1;

   btas::IVector<K> braContracts;
   braContracts[indexRij] = 0;
   for(size_t i = 0;            i < indexRij; ++i) braContracts[i] = i + 1;
   for(size_t i = indexRij + 1; i < K;   ++i) braContracts[i] = i;

   btas::IVector<K> ketContracts;
   for(size_t i = 0; i < K; ++i) ketContracts[i] = i;

   for(size_t i = 0; i < ij_loop_index.size(); ++i)
   {
      size_t ix = ij_loop_index[i];
      // Qii
      if(blockRij->CreDesComp(ix, ix).size() > 0)
      {
         btas::SDArray<N> iiCreDesCompScr;
         btas::Contract( 1.0, blockRij->CreDesComp(ix, ix), btas::shape(0), bra, btas::shape(indexRij), 1.0, iiCreDesCompScr);
         btas::Contract( 1.0, iiCreDesCompScr, braContracts, ket, ketContracts, 1.0, newBlock->CreDesComp(ix, ix));
      }

      for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
      {
         size_t jx = ij_loop_index[j];
         // Pij
         if(blockRij->CreCreComp(ix, jx).size() > 0)
         {
            btas::SDArray<N> ijCreCreCompScr;
            btas::Contract( 1.0, blockRij->CreCreComp(ix, jx), btas::shape(0), bra, btas::shape(indexRij), 1.0, ijCreCreCompScr);
            btas::Contract( 1.0, ijCreCreCompScr, braContracts, ket, ketContracts, 1.0, newBlock->CreCreComp(ix, jx));
            ScaledCopy(-1.0, newBlock->CreCreComp(ix, jx), newBlock->CreCreComp(jx, ix));
         }
         // Qij
         if(blockRij->CreDesComp(ix, jx).size() > 0)
         {
            btas::SDArray<N> ijCreDesCompScr;
            btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(0), bra, btas::shape(indexRij), 1.0, ijCreDesCompScr);
            btas::Contract( 1.0, ijCreDesCompScr, braContracts, ket, ketContracts, 1.0, newBlock->CreDesComp(ix, jx));
            btas::Permute(newBlock->CreDesComp(ix, jx), btas::shape(1, 0), newBlock->CreDesComp(jx, ix));
         }
      }
   }
}

void RijBlockCkl (
      const std::vector<size_t>& kl_loop_index,
            Block* newBlock,
      const std::vector<size_t>& ij_loop_index)
{
   for(size_t k = 0; k < kl_loop_index.size(); ++k)
   {
      size_t kx = kl_loop_index[k];
      int kSpin = 1; if(kx % 2 != 0) kSpin = -1;

      DmrgOperator BkkScr = newBlock->CreDes(kx, kx);

      for(size_t i = 0; i < ij_loop_index.size(); ++i)
      {
         size_t ix = ij_loop_index[i];
         int iSpin = 1; if(ix % 2 != 0) iSpin = -1;

         // Qii
         DmrgOperator& QiiRef = newBlock->CreDesComp(ix, ix);

         double vikik = integrals::moeri::TwoInt(ix, kx, ix, kx);
         if(fabs(vikik) >= THRE_INT)
            btas::Axpy( vikik, BkkScr, QiiRef);

         double viikk = integrals::moeri::TwoInt(ix, ix, kx, kx);
         if(fabs(viikk) >= THRE_INT && (iSpin == kSpin))
            btas::Axpy(-viikk, BkkScr, QiiRef);

         for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
         {
            size_t jx = ij_loop_index[j];
            int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

            // Qij
            DmrgOperator& QijRef = newBlock->CreDesComp(ix, jx);

            double vikjk = integrals::moeri::TwoInt(ix, kx, jx, kx);
            if(fabs(vikjk) >= THRE_INT && (iSpin == jSpin))
               btas::Axpy( vikjk, BkkScr, QijRef);

            double vijkk = integrals::moeri::TwoInt(ix, jx, kx, kx);
            if(fabs(vijkk) >= THRE_INT && (iSpin == kSpin) && (jSpin == kSpin))
               btas::Axpy(-vijkk, BkkScr, QijRef);
         }
      }
   }

   for(size_t k = 0; k < kl_loop_index.size(); ++k)
   {
      size_t kx = kl_loop_index[k];
      int kSpin = 1; if(kx % 2 != 0) kSpin = -1;

      for(size_t l = k + 1; l < kl_loop_index.size(); ++l)
      {
         size_t lx = kl_loop_index[l];
         int lSpin = 1; if(lx % 2 != 0) lSpin = -1;

         DmrgOperator AklScr = newBlock->CreCre(kx, lx);
         DmrgOperator BklScr = newBlock->CreDes(kx, lx);
         DmrgOperator BlkScr;
         btas::Permute(BklScr, btas::shape(1, 0), BlkScr);

         for(size_t i = 0; i < ij_loop_index.size(); ++i)
         {
            size_t ix = ij_loop_index[i];
            int iSpin = 1; if(ix % 2 != 0) iSpin = -1;

            // Qii
            DmrgOperator& QiiRef = newBlock->CreDesComp(ix, ix);

            double vikil = integrals::moeri::TwoInt(ix, kx, ix, lx);
            if(fabs(vikil) >= THRE_INT && (kSpin == lSpin)) // iSpin == iSpin is trivial
            {
               btas::Axpy( vikil, BklScr, QiiRef);
               btas::Axpy( vikil, BlkScr, QiiRef);
            }

            double viikl = integrals::moeri::TwoInt(ix, ix, kx, lx);
            if(fabs(viikl) >= THRE_INT && (iSpin == kSpin) && (iSpin == lSpin))
            {
               btas::Axpy(-viikl, BklScr, QiiRef);
               btas::Axpy(-viikl, BlkScr, QiiRef);
            }

            for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
            {
               size_t jx = ij_loop_index[j];
               int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

               // Pij & Qij
               DmrgOperator& QijRef = newBlock->CreDesComp(ix, jx);
               DmrgOperator& PijRef = newBlock->CreCreComp(ix, jx);

               double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
               if(fabs(vikjl) >= THRE_INT && (iSpin == jSpin) && (kSpin == lSpin))
               {
                  btas::Axpy( vikjl, BklScr, QijRef);
                  btas::Axpy( vikjl, BlkScr, QijRef);
               }

               double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
               if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
               {
                  btas::Axpy(-vijkl, AklScr, PijRef);
                  btas::Axpy(-vijkl, BlkScr, QijRef);
               }

               double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
               if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
               {
                  btas::Axpy( vijlk, AklScr, PijRef);
                  btas::Axpy(-vijlk, BklScr, QijRef);
               }
            }
         }
      }
   }

   for(size_t i = 0; i < ij_loop_index.size(); ++i)
   {
      size_t ix = ij_loop_index[i];

      for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
      {
         size_t jx = ij_loop_index[j];

         // Pji = -Pij
         if(newBlock->CreCreComp(ix, jx).size() > 0)
            ScaledCopy(-1.0, newBlock->CreCreComp(ix, jx),
                                   newBlock->CreCreComp(jx, ix));
         // Qji = Qij'
         if(newBlock->CreDesComp(ix, jx).size() > 0)
            btas::Permute(newBlock->CreDesComp(ix, jx), btas::shape(1, 0), newBlock->CreDesComp(jx, ix));
      }
   }
}

void RijBlockCkCl (
      const std::vector<size_t>& k_loop_index,
      const std::vector<size_t>& l_loop_index,
            Block* newBlock,
      const std::vector<size_t>& ij_loop_index)
{
   for(size_t k = 0; k < k_loop_index.size(); ++k)
   {
      size_t kx = k_loop_index[k];
      int kSpin = 1; if(kx % 2 != 0) kSpin = -1;

      for(size_t l = 0; l < l_loop_index.size(); ++l)
      {
         size_t lx = l_loop_index[l];
         int lSpin = 1; if(lx % 2 != 0) lSpin = -1;

         DmrgOperator AklScr = newBlock->CreCre(kx, lx);
         DmrgOperator BklScr = newBlock->CreDes(kx, lx);
         DmrgOperator BlkScr;
         btas::Permute(BklScr, btas::shape(1, 0), BlkScr);

         for(size_t i = 0; i < ij_loop_index.size(); ++i)
         {
            size_t ix = ij_loop_index[i];
            int iSpin = 1; if(ix % 2 != 0) iSpin = -1;

            // Qii
            DmrgOperator& QiiRef = newBlock->CreDesComp(ix, ix);

            double vikil = integrals::moeri::TwoInt(ix, kx, ix, lx);
            if(fabs(vikil) >= THRE_INT && (kSpin == lSpin)) // iSpin == iSpin is trivial
            {
               btas::Axpy( vikil, BklScr, QiiRef);
               btas::Axpy( vikil, BlkScr, QiiRef);
            }

            double viikl = integrals::moeri::TwoInt(ix, ix, kx, lx);
            if(fabs(viikl) >= THRE_INT && (iSpin == kSpin) && (iSpin == lSpin))
            {
               btas::Axpy(-viikl, BklScr, QiiRef);
               btas::Axpy(-viikl, BlkScr, QiiRef);
            }

            for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
            {
               size_t jx = ij_loop_index[j];
               int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

               // Pij & Qij
               DmrgOperator& QijRef = newBlock->CreDesComp(ix, jx);
               DmrgOperator& PijRef = newBlock->CreCreComp(ix, jx);

               double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
               if(fabs(vikjl) >= THRE_INT && (iSpin == jSpin) && (kSpin == lSpin))
               {
                  btas::Axpy( vikjl, BklScr, QijRef);
                  btas::Axpy( vikjl, BlkScr, QijRef);
               }

               double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
               if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
               {
                  btas::Axpy(-vijkl, AklScr, PijRef);
                  btas::Axpy(-vijkl, BlkScr, QijRef);
               }

               double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
               if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
               {
                  btas::Axpy( vijlk, AklScr, PijRef);
                  btas::Axpy(-vijlk, BklScr, QijRef);
               }
            }
         }
      }
   }

   for(size_t i = 0; i < ij_loop_index.size(); ++i)
   {
      size_t ix = ij_loop_index[i];

      for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
      {
         size_t jx = ij_loop_index[j];

         // Pji = -Pij
         if(newBlock->CreCreComp(ix, jx).size() > 0)
            ScaledCopy(-1.0, newBlock->CreCreComp(ix, jx), newBlock->CreCreComp(jx, ix));

         // Qji = Qij'
         if(newBlock->CreDesComp(ix, jx).size() > 0)
            btas::Permute(newBlock->CreDesComp(ix, jx), btas::shape(1, 0), newBlock->CreDesComp(jx, ix));
      }
   }
}

} // namespace ttns

#endif // __TTNS_OP_COMPONENTS_COMP_H

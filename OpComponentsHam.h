#ifndef __TTNS_OP_COMPONENTS_HAM_H
#define __TTNS_OP_COMPONENTS_HAM_H

#include "OpComponentsUtil.h"
#include "ttns_assert.h"

namespace ttns
{

#ifndef THRE_INT
#define THRE_INT 1.0e-16
#endif

//==========================================================================================//
// Compute Ham, Cre, CreCre, and CreDes blocks for sigma vector compt. and renormalization. //
//==========================================================================================//

//
// Ham Block from Local H
//
template<size_t N, class Q>
void HamHam (
      const btas::QSDArray<N, Q>& bra,
            Block* blockHam,
      const size_t& indexHam,
            btas::QSDArray<N, Q>& hamHam)
{
   TTNS_DEBUG("HamHam : Started.");

   btas::SDArray<N> hamScr;
   btas::Contract( 1.0, blockHam->Ham(), btas::shape(0), bra, btas::shape(indexHam), 1.0, hamScr);

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexHam);
   PermutedAxpy( 1.0, hamScr, reorder, hamHam);

   TTNS_DEBUG("HamHam : Finished.");
}

//
// H := bra x [ Ci x Ri + Rj x Cj ] ( algorithm for small block )
//
template<size_t N, size_t K, class Q>
void HamCiCjS (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiCj)
{
   TTNS_DEBUG("HamCiCjS : Started.");

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   double opSign = 1.0; std::vector<int> p1;
   IndexParity(opSign, p1, order, indexCi, indexCj);
   braCopy.parity(p1);

   btas::SDArray<N> hamScr;

   btas::SDArray<4> hamSuperBlock;
   // Ci x Ri
   btas::SDArray<4> hamSuperBlockScr;
   for(size_t i = 0; i < i_loop_index.size(); ++i)
   {
      size_t ix = i_loop_index[i];
      if(blockCj->CreComp(ix).size() == 0) continue;
      btas::Ger(+opSign, blockCi->Cre(ix), blockCj->CreComp(ix), hamSuperBlockScr);
   }
   // Rj x Cj
   for(size_t j = 0; j < j_loop_index.size(); ++j)
   {
      size_t jx = j_loop_index[j];
      if(blockCi->CreComp(jx).size() == 0) continue;
      btas::Ger(+opSign, blockCi->CreComp(jx), blockCj->Cre(jx), hamSuperBlockScr);
   }

   if(hamSuperBlockScr.size() > 0)
   {
      // CreComp x Des
      PermutedAxpy( 1.0, hamSuperBlockScr, btas::shape(0, 1, 3, 2), hamSuperBlock);
      // DesComp x Cre
      PermutedAxpy(-1.0, hamSuperBlockScr, btas::shape(1, 0, 2, 3), hamSuperBlock);

      btas::Contract( 1.0, hamSuperBlock, btas::shape(0, 2), braCopy, btas::shape(indexCi, indexCj), 1.0, hamScr);

   TTNS_DEBUG("HamCiCjS : indexCi = " << indexCi << ", indexCj = " << indexCj);
      btas::IVector<N> reorder;
      IndexPermute(reorder, indexCj, indexCi);
      PermutedAxpy( 1.0, hamScr, reorder, hamCiCj);
   }

   TTNS_DEBUG("HamCiCjS : Finished.");
}

//
// H := Ci x [ Cj x bra ] ( algorithm for large block )
//
template<size_t N, size_t K, class Q>
void HamCiCjL (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiCj)
{
   TTNS_DEBUG("HamCiCjL : Started.");
   TTNS_ASSERT(indexCi <= indexCj, "ttns::HamCiCjL; indices must be indexCi <= indexCj");

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   double opSign = 1.0; std::vector<int> p1;
   IndexParity(opSign, p1, order, indexCi, indexCj);
   braCopy.parity(p1);

   btas::SDArray<N> hamScr;
   // Ci x Ri
   for(size_t i = 0; i < i_loop_index.size(); ++i)
   {
      size_t ix = i_loop_index[i];
      if(blockCj->CreComp(ix).size() == 0) continue;
      // Cre x DesComp
      {
         btas::SDArray<N> CreScr;
         btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0), braCopy, btas::shape(indexCi), 1.0, CreScr);
         btas::Contract(+opSign, blockCj->CreComp(ix), btas::shape(1), CreScr, btas::shape(indexCj), 1.0, hamScr);
      }
      // Des x CreComp
      {
         btas::SDArray<N> DesScr;
         btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(1), braCopy, btas::shape(indexCi), 1.0, DesScr);
         btas::Contract(-opSign, blockCj->CreComp(ix), btas::shape(0), DesScr, btas::shape(indexCj), 1.0, hamScr);
      }
   }
   // Rj x Cj
   for(size_t j = 0; j < j_loop_index.size(); ++j)
   {
      size_t jx = j_loop_index[j];
      if(blockCi->CreComp(jx).size() == 0) continue;
      // CreComp x Des
      {
         btas::SDArray<N> CreCompScr;
         btas::Contract( 1.0, blockCi->CreComp(jx), btas::shape(0), braCopy, btas::shape(indexCi), 1.0, CreCompScr);
         btas::Contract(+opSign, blockCj->Cre(jx), btas::shape(1), CreCompScr, btas::shape(indexCj), 1.0, hamScr);
      }
      // DesComp x Cre
      {
         btas::SDArray<N> DesCompScr;
         btas::Contract( 1.0, blockCi->CreComp(jx), btas::shape(1), braCopy, btas::shape(indexCi), 1.0, DesCompScr);
         btas::Contract(-opSign, blockCj->Cre(jx), btas::shape(0), DesCompScr, btas::shape(indexCj), 1.0, hamScr);
      }
   }

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCi, indexCj);
   PermutedAxpy( 1.0, hamScr, reorder, hamCiCj);

   TTNS_DEBUG("HamCiCjL : Finished.");
}

//
// H := Ci x DesCompRi - Di x CreCompRi ( using pre-contracted [ Ri x bra ] )
//
template<size_t N, size_t K, class Q>
void HamCiRi (
      const btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
      const btas::TArray<btas::SDArray<N>, 1>& DesCompRi,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiRi)
{
   TTNS_DEBUG("HamCiRi : Started.");

   btas::SDArray<N> hamScr;
   // Cre x DesComp
   if(CreCompRi.size() > 0)
   {
      for(size_t i = 0; i < i_loop_index.size(); ++i)
      {
         size_t ix = i_loop_index[i];
         if(DesCompRi(i).size() > 0)
            btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0), DesCompRi(i), btas::shape(indexCi), 1.0, hamScr);
      }
   }
   // Des x CreComp
   if(DesCompRi.size() > 0)
   {
      for(size_t i = 0; i < i_loop_index.size(); ++i)
      {
         size_t ix = i_loop_index[i];
         if(CreCompRi(i).size() > 0)
            btas::Contract(-1.0, blockCi->Cre(ix), btas::shape(1), CreCompRi(i), btas::shape(indexCi), 1.0, hamScr);
      }
   }

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCi);
   PermutedAxpy( 1.0, hamScr, reorder, hamCiRi);

   TTNS_DEBUG("HamCiRi : Finished.");
}

//
// H := bra x [ Aij x Pij + Bij x Qij ] ( algorithm for small block )
//
template<size_t N, size_t K, class Q>
void HamCijRijS (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCij,
      const size_t& indexCij,
      const std::vector<size_t>& ij_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCijRij)
{
   TTNS_DEBUG("HamCijRijS : Started.");

   btas::SDArray<4> hamSuperBlock;

   btas::SDArray<4> hamSuperBlockBiiQii;
   btas::SDArray<4> hamSuperBlockAijPij;
   btas::SDArray<4> hamSuperBlockBijQij;

   for(size_t i = 0; i < ij_loop_index.size(); ++i)
   {
      size_t ix = ij_loop_index[i];
      // Bii x Qii
      if(blockRij->CreDesComp(ix, ix).size() > 0)
      {
         btas::Ger( 1.0, blockCij->CreDes(ix, ix), blockRij->CreDesComp(ix, ix), hamSuperBlockBiiQii);
      }

      for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
      {
         size_t jx = ij_loop_index[j];
         // Aij x Pij
         if(blockRij->CreCreComp(ix, jx).size() > 0)
         {
            btas::Ger(-1.0, blockCij->CreCre(ix, jx), blockRij->CreCreComp(ix, jx), hamSuperBlockAijPij);
         }
         // Bij x Qij
         if(blockRij->CreDesComp(ix, jx).size() > 0)
         {
            btas::Ger( 1.0, blockCij->CreDes(ix, jx), blockRij->CreDesComp(ix, jx), hamSuperBlockBijQij);
         }
      }
   }
   if(hamSuperBlockBiiQii.size() > 0)
   {
      // H := [ Bii x Qii ]
      btas::Axpy ( 1.0, hamSuperBlockBiiQii, hamSuperBlock);
   }
   if(hamSuperBlockAijPij.size() > 0)
   {
      // H := [ Aij x Pij' ] + [ Aij' x Pij ]
      PermutedAxpy( 1.0, hamSuperBlockAijPij, btas::shape(0, 1, 3, 2), hamSuperBlock);
      PermutedAxpy( 1.0, hamSuperBlockAijPij, btas::shape(1, 0, 2, 3), hamSuperBlock);
   }
   if(hamSuperBlockBijQij.size() > 0)
   {
      // H := [ Bij x Qij ] + [ Bij' x Qij' ]
      btas::Axpy ( 1.0, hamSuperBlockBijQij, hamSuperBlock);
      PermutedAxpy( 1.0, hamSuperBlockBijQij, btas::shape(1, 0, 3, 2), hamSuperBlock);
   }

   if(hamSuperBlock.size() > 0)
   {
      btas::SDArray<N> hamScr;
      btas::Contract( 1.0, hamSuperBlock, btas::shape(0, 2), bra, btas::shape(indexCij, indexRij), 1.0, hamScr);

      btas::IVector<N> reorder;
      IndexPermute(reorder, indexRij, indexCij);
      PermutedAxpy( 1.0, hamScr, reorder, hamCijRij);
   }

   TTNS_DEBUG("HamCijRijS : Finished.");
}

//
// H := Aij x [ Pij x bra ] + Bij x [ Qij x bra ] ( algorithm for large block )
//
template<size_t N, size_t K, class Q>
void HamCijRijL(const btas::QSDArray<N, Q>& bra,
                        Block* blockCij, const size_t& indexCij, const std::vector<size_t>& ij_loop_index,
                        Block* blockRij, const size_t& indexRij,
                  const btas::IVector<K>& order,
                        btas::QSDArray<N, Q>& hamCijRij)
{
   TTNS_DEBUG("HamCijRijL : Started.");
   TTNS_ASSERT(indexCij <= indexRij, "HamCijRijL; indices must be indexCij <= indexRij");

   btas::SDArray<N> hamScr;

   for(size_t i = 0; i < ij_loop_index.size(); ++i)
   {
      size_t ix = ij_loop_index[i];
      if(blockRij->CreDesComp(ix, ix).size() > 0)
      {
         // Bii x Qii
         btas::SDArray<N> CreDesScr;
         btas::Contract( 1.0, blockCij->CreDes(ix, ix), btas::shape(0), bra, btas::shape(indexCij), 1.0, CreDesScr);
         btas::Contract( 1.0, blockRij->CreDesComp(ix, ix), btas::shape(0), CreDesScr, btas::shape(indexRij), 1.0, hamScr);
      }

      for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
      {
         size_t jx = ij_loop_index[j];
         // Aij x Pij
         if(blockRij->CreCreComp(ix, jx).size() > 0)
         {
            btas::SDArray<N> CreCreScr;
            btas::Contract( 1.0, blockCij->CreCre(ix, jx), btas::shape(0), bra, btas::shape(indexCij), 1.0, CreCreScr);
            btas::Contract(-1.0, blockRij->CreCreComp(ix, jx), btas::shape(1), CreCreScr, btas::shape(indexRij), 1.0, hamScr);

            btas::SDArray<N> DesDesScr;
            btas::Contract(-1.0, blockCij->CreCre(ix, jx), btas::shape(1), bra, btas::shape(indexCij), 1.0, DesDesScr);
            btas::Contract( 1.0, blockRij->CreCreComp(ix, jx), btas::shape(0), DesDesScr, btas::shape(indexRij), 1.0, hamScr);
         }

         // Bij x Qij
         if(blockRij->CreDesComp(ix, jx).size() > 0)
         {
            btas::SDArray<N> CreDesScr;
            btas::Contract( 1.0, blockCij->CreDes(ix, jx), btas::shape(0), bra, btas::shape(indexCij), 1.0, CreDesScr);
            btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(0), CreDesScr, btas::shape(indexRij), 1.0, hamScr);

            btas::SDArray<N> DesCreScr;
            btas::Contract( 1.0, blockCij->CreDes(ix, jx), btas::shape(1), bra, btas::shape(indexCij), 1.0, DesCreScr);
            btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(1), DesCreScr, btas::shape(indexRij), 1.0, hamScr);
         }
      }
   }

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCij, indexRij);
   PermutedAxpy( 1.0, hamScr, reorder, hamCijRij);

   TTNS_DEBUG("HamCijRijL : Finished.");
}

//
// H := Ci x [ bra x [ Cj x Pij + Dj x Qij ] ] ( algorithm for small block )
//
template<size_t N, size_t K, class Q>
void HamCiCjRijS (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiCjRij)
{
   TTNS_DEBUG("HamCiCjRijS : Started.");

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   double opSign = 1.0; std::vector<int> p1;
   IndexParity(opSign, p1, order, indexCi, indexCj);
   braCopy.parity(p1);

   size_t iOffCi  = 0;
   if(indexCj  > indexCi) iOffCi++;
   if(indexRij > indexCi) iOffCi++;

   btas::SDArray<N> hamScr;

   for(size_t i = 0; i < i_loop_index.size(); ++i)
   {
      size_t ix = i_loop_index[i];

      btas::SDArray<4> SuperBlockPij;
      btas::SDArray<4> SuperBlockQij;

      for(size_t j = 0; j < j_loop_index.size(); ++j)
      {
         size_t jx = j_loop_index[j];
         // R(i) = Cj x Pij + Dj x Qij
         if(blockRij->CreCreComp(ix, jx).size() > 0)
         {
            btas::Ger( 1.0, blockCj->Cre(jx), blockRij->CreCreComp(ix, jx), SuperBlockPij);
         }

         if(blockRij->CreDesComp(ix, jx).size() > 0)
         {
            btas::Ger( 1.0, blockCj->Cre(jx), blockRij->CreDesComp(ix, jx), SuperBlockQij);
         }
      }

      btas::SDArray<4> CreCompSuperBlock;

      if(SuperBlockPij.size() > 0)
      {
         PermutedAxpy(-1.0, SuperBlockPij, btas::shape(1, 0, 2, 3), CreCompSuperBlock);
      }
      if(SuperBlockQij.size() > 0)
      {
         PermutedAxpy( 1.0, SuperBlockQij, btas::shape(0, 1, 3, 2), CreCompSuperBlock);
      }

      if(CreCompSuperBlock.size() == 0) continue;

      btas::SDArray<N> CreCompScr;
      btas::Contract( 1.0, CreCompSuperBlock, btas::shape(0, 2), braCopy, btas::shape(indexCj, indexRij), 1.0, CreCompScr);
      btas::Contract(-opSign, blockCi->Cre(ix), btas::shape(1), CreCompScr, btas::shape(indexCi+iOffCi), 1.0, hamScr);

      btas::SDArray<N> DesCompScr;
      btas::Contract( 1.0, CreCompSuperBlock, btas::shape(1, 3), braCopy, btas::shape(indexCj, indexRij), 1.0, DesCompScr);
      btas::Contract(+opSign, blockCi->Cre(ix), btas::shape(0), DesCompScr, btas::shape(indexCi+iOffCi), 1.0, hamScr);
   }

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexRij, indexCj, indexCi);
   PermutedAxpy( 1.0, hamScr, reorder, hamCiCjRij);

   TTNS_DEBUG("HamCiCjRijS : Finished.");
}

//
// H := Ci x [ Pij x [ Cj x bra ] + Qij x [ Dj x bra ] ] ( algorithm for large block )
//
template<size_t N, size_t K, class Q>
void HamCiCjRijL (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiCjRij)
{
   TTNS_DEBUG("HamCiCjRijL : Started.");

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   double opSign = 1.0; std::vector<int> p1;
   IndexParity(opSign, p1, order, indexCi, indexCj);
   braCopy.parity(p1);

   size_t iOffRij = 0;
   if(indexCj  > indexRij) iOffRij++;

   size_t iOffCi  = 0;
   if(indexCj  > indexCi) iOffCi++;
   if(indexRij > indexCi) iOffCi++;

   btas::SDArray<N> hamScr;

   btas::TArray<btas::SDArray<N>, 1 > CreCompScr(i_loop_index.size());
   btas::TArray<btas::SDArray<N>, 1 > DesCompScr(i_loop_index.size());

   // R(i) = Pij x [ Cj x bra ] + Qij x [ Dj x bra ]
   // R(i)'= Pij'x [ Dj x bra ] + Qij'x [ Cj x bra ]
   for(size_t j = 0; j < j_loop_index.size(); ++j)
   {
      size_t jx = j_loop_index[j];
      // Cj x bra
      btas::SDArray<N> CreScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0), braCopy, btas::shape(indexCj), 1.0, CreScr);
      // Dj x bra
      btas::SDArray<N> DesScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(1), braCopy, btas::shape(indexCj), 1.0, DesScr);

      for(size_t i = 0; i < i_loop_index.size(); ++i)
      {
        size_t ix = i_loop_index[i];
        if(blockRij->CreCreComp(ix, jx).size() > 0)
        {
           // Ri := Pij x Cj
           btas::Contract(-1.0, blockRij->CreCreComp(ix, jx), btas::shape(1), CreScr, btas::shape(indexRij+iOffRij), 1.0, DesCompScr(i));
           // Ri':= Pij'x Dj
           btas::Contract(-1.0, blockRij->CreCreComp(ix, jx), btas::shape(0), DesScr, btas::shape(indexRij+iOffRij), 1.0, CreCompScr(i));
        }
        if(blockRij->CreDesComp(ix, jx).size() > 0)
        {
           // Ri := Qij x Dj
           btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(0), DesScr, btas::shape(indexRij+iOffRij), 1.0, DesCompScr(i));
           // Ri':= Qij'x Cj
           btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(1), CreScr, btas::shape(indexRij+iOffRij), 1.0, CreCompScr(i));
        }
      }
   }
   // H := Ci x R(i) - Di x Ri(i)'
   for(size_t i = 0; i < i_loop_index.size(); ++i)
   {
      size_t ix = i_loop_index[i];
      if(DesCompScr(i).size() > 0)
      {
         btas::Contract(+opSign, blockCi->Cre(ix), btas::shape(0), DesCompScr(i), btas::shape(indexCi+iOffCi), 1.0, hamScr);
      }
      if(CreCompScr(i).size() > 0)
      {
         btas::Contract(-opSign, blockCi->Cre(ix), btas::shape(1), CreCompScr(i), btas::shape(indexCi+iOffCi), 1.0, hamScr);
      }
   }
   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCj, indexRij, indexCi);
   PermutedAxpy( 1.0, hamScr, reorder, hamCiCjRij);

   TTNS_DEBUG("HamCiCjRijL : Finished.");
}

//
// Ri := bra x [ Cj x Pij + Dj x Qij ] ( algorithm for small block )
//
template<size_t N, size_t K, class Q>
void CreCompCjRijS (
      const btas::QSDArray<N, Q>& bra,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<K>& order,
            btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
            btas::TArray<btas::SDArray<N>, 1>& DesCompRi)
{
   TTNS_DEBUG("CreCompCjRijS : Started.");

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   double opSign = 1.0; std::vector<int> p1;
   IndexParity(opSign, p1, order, indexCi, indexCj);
   braCopy.parity(p1);

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexRij, indexCj);

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();

   if(CreCompRi.size() == 0) CreCompRi.resize(ni);
   if(DesCompRi.size() == 0) DesCompRi.resize(ni);

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
         {
            btas::Ger( 1.0, blockCj->Cre(jx), blockRij->CreCreComp(ix, jx), SuperBlockPij);
         }

         if(blockRij->CreDesComp(ix, jx).size() > 0)
         {
            btas::Ger( 1.0, blockCj->Cre(jx), blockRij->CreDesComp(ix, jx), SuperBlockQij);
         }
      }

      btas::SDArray<4> CreCompSuperBlock;
      if(SuperBlockPij.size() > 0)
      {
         PermutedAxpy(-1.0, SuperBlockPij, btas::shape(1, 0, 2, 3), CreCompSuperBlock);
      }
      if(SuperBlockQij.size() > 0)
      {
         PermutedAxpy( 1.0, SuperBlockQij, btas::shape(0, 1, 3, 2), CreCompSuperBlock);
      }

      if(CreCompSuperBlock.size() == 0) continue;

      btas::SDArray<N> CreCompScr;
//    btas::Contract(-opSign, CreCompSuperBlock, btas::shape(0, 2),  /*   Di x CreCompRi */
      btas::Contract(+opSign, CreCompSuperBlock, btas::shape(0, 2),  /* - Di x CreCompRi */
                              braCopy, btas::shape(indexCj, indexRij), 1.0, CreCompScr);
      PermutedAxpy( 1.0, CreCompScr, reorder, CreCompRi(i));

      btas::SDArray<N> DesCompScr;
      btas::Contract(+opSign, CreCompSuperBlock, btas::shape(1, 3), braCopy, btas::shape(indexCj, indexRij), 1.0, DesCompScr);
      PermutedAxpy( 1.0, DesCompScr, reorder, DesCompRi(i));
   }

   TTNS_DEBUG("CreCompCjRijS : Finished.");
}

//
// Ri := Pij x [ Cj x bra ] + Qij x [ Dj x bra ] ( algorithm for large block )
//
template<size_t N, size_t K, class Q>
void CreCompCjRijL (
      const btas::QSDArray<N, Q>& bra,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockRij,
      const size_t& indexRij,
      const btas::IVector<K>& order,
            btas::TArray<btas::SDArray<N>, 1>& CreCompRi,
            btas::TArray<btas::SDArray<N>, 1>& DesCompRi)
{
   TTNS_DEBUG("CreCompCjRijL : Started.");

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   double opSign = 1.0; std::vector<int> p1;
   IndexParity(opSign, p1, order, indexCi, indexCj);
   braCopy.parity(p1);

   size_t iOffRij = 0;
   if(indexCj  > indexRij) iOffRij++;

   btas::IVector<N> reorder;
   IndexPermute(reorder, indexCj, indexRij);

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();

   if(CreCompRi.size() == 0) CreCompRi.resize(ni);
   if(DesCompRi.size() == 0) DesCompRi.resize(ni);

   btas::TArray<btas::SDArray<N>, 1> CreCompScr(ni);
   btas::TArray<btas::SDArray<N>, 1> DesCompScr(ni);

   // R(i) = Pij x [ Cj x bra ] + Qij x [ Dj x bra ]
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
         {
            // Ri := Pij x Cj
            btas::Contract(-1.0, blockRij->CreCreComp(ix, jx), btas::shape(1), CreScr, btas::shape(indexRij+iOffRij), 1.0, DesCompScr(i));
            // Ri':= Pij'x Dj
            btas::Contract(-1.0, blockRij->CreCreComp(ix, jx), btas::shape(0), DesScr, btas::shape(indexRij+iOffRij), 1.0, CreCompScr(i));
         }
         if(blockRij->CreDesComp(ix, jx).size() > 0)
         {
            // Ri := Qij x Dj
            btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(0), DesScr, btas::shape(indexRij+iOffRij), 1.0, DesCompScr(i));
            // Ri':= Qij'x Cj
            btas::Contract( 1.0, blockRij->CreDesComp(ix, jx), btas::shape(1), CreScr, btas::shape(indexRij+iOffRij), 1.0, CreCompScr(i));
         }
      }
   }

   for(size_t i = 0; i < ni; ++i)
   {
      if(CreCompScr(i).size() > 0)
//       PermutedAxpy(-opSign, CreCompScr(i), reorder, CreCompRi(i)); // + Di x CreCompRi
         PermutedAxpy(+opSign, CreCompScr(i), reorder, CreCompRi(i)); // - Di x CreCompRi
      if(DesCompScr(i).size() > 0)
         PermutedAxpy(+opSign, DesCompScr(i), reorder, DesCompRi(i)); // + Ci x DesCompRi
   }

   TTNS_DEBUG("CreCompCjRijL : Finished.");
}

//
// H := Cl x [ Ci x [ vijkl [ Ck x [ Cj x bra ] ] ] ] ( Ci, Cj are in small block )
//
template<size_t N, size_t K, class Q>
void HamCiCjDkDlS (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockDk,
      const size_t& indexDk,
      const std::vector<size_t>& k_loop_index,
            Block* blockDl,
      const size_t& indexDl,
      const std::vector<size_t>& l_loop_index,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiCjDkDl)
{
   TTNS_DEBUG("HamCiCjDkDlS : Started.");

   // must be indexCi < indexCj < indexDk < indexDl
   {
      std::vector<size_t> originalOrder;
      originalOrder.push_back(indexCi);
      originalOrder.push_back(indexCj);
      originalOrder.push_back(indexDk);
      originalOrder.push_back(indexDl);

      std::vector<size_t> sortedOrder(originalOrder);
      std::sort(sortedOrder.begin(), sortedOrder.end());

      TTNS_ASSERT(std::equal(sortedOrder.begin(), sortedOrder.end(), originalOrder.begin()),
      "ttns::HamCiCjDkDlL; contraction order must be sorted before calling");
   }

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   std::vector<int> p1;
   {
      std::vector<size_t> bounds;
      bounds.push_back(order[indexCi]);
      bounds.push_back(order[indexCj]);
      bounds.push_back(order[indexDk]);
      bounds.push_back(order[indexDl]);
      std::sort(bounds.begin(), bounds.end());

      for(size_t iBlock = 0; iBlock < K; ++iBlock)
      {
         if(bounds[0] < order[iBlock] && order[iBlock] <= bounds[1])
            p1.push_back(iBlock);
         if(bounds[2] < order[iBlock] && order[iBlock] <= bounds[3])
            p1.push_back(iBlock);
      }
   }
   braCopy.parity(p1);

   btas::SDArray<N> hamScr;

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();
   size_t nk = k_loop_index.size();
   size_t nl = l_loop_index.size();

   btas::TArray<btas::SDArray<N>, 2> ilCreCreComp(ni, nl);
   btas::TArray<btas::SDArray<N>, 2> ilCreDesComp(ni, nl);
   btas::TArray<btas::SDArray<N>, 2> ilDesCreComp(ni, nl);
   btas::TArray<btas::SDArray<N>, 2> ilDesDesComp(ni, nl);

   double opSign = 1.0;
   if(order[indexCi] > order[indexCj]) opSign = -opSign;
   if(order[indexCi] > order[indexDk]) opSign = -opSign;
   if(order[indexCi] > order[indexDl]) opSign = -opSign;
   if(order[indexCj] > order[indexDk]) opSign = -opSign;
   if(order[indexCj] > order[indexDl]) opSign = -opSign;
   if(order[indexDk] > order[indexDl]) opSign = -opSign;

   for(size_t k = 0; k < nk; ++k)
   {
      size_t kx = k_loop_index[k];
      int kSpin = 1.0; if(kx % 2 != 0) kSpin = -1.0;
      btas::SDArray<N> kCreScr;
      btas::Contract( 1.0, blockDk->Cre(kx), btas::shape(0), braCopy, btas::shape(indexDk), 1.0, kCreScr);

      btas::SDArray<N> kDesScr;
      btas::Contract( 1.0, blockDk->Cre(kx), btas::shape(1), braCopy, btas::shape(indexDk), 1.0, kDesScr);

      for(size_t j = 0; j < nj; ++j)
      {
         size_t jx = j_loop_index[j];
         int jSpin = 1.0; if(jx % 2 != 0) jSpin = -1.0;
         btas::SDArray<N> jkCreCreScr;
         btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0), kCreScr, btas::shape(indexCj+1), 1.0, jkCreCreScr);

         btas::SDArray<N> jkCreDesScr;
         btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0), kDesScr, btas::shape(indexCj+1), 1.0, jkCreDesScr);

         btas::SDArray<N> jkDesCreScr;
         btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(1), kCreScr, btas::shape(indexCj+1), 1.0, jkDesCreScr);

         btas::SDArray<N> jkDesDesScr;
         btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(1), kDesScr, btas::shape(indexCj+1), 1.0, jkDesDesScr);

         for(size_t i = 0; i < ni; ++i)
         {
            size_t ix = i_loop_index[i];
            int iSpin = 1.0; if(ix % 2 != 0) iSpin = -1.0;
            for(size_t l = 0; l < nl; ++l)
            {
               size_t lx = l_loop_index[l];
               int lSpin = 1.0; if(lx % 2 != 0) lSpin = -1.0;

               double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
               if(fabs(vikjl) >= THRE_INT && (iSpin == jSpin) && (kSpin == lSpin))
               {
                  // Pil
                  btas::Axpy(-opSign*vikjl, jkCreCreScr, ilCreCreComp(i, l));
                  btas::Axpy(-opSign*vikjl, jkDesDesScr, ilDesDesComp(i, l));
                  // Qil
                  btas::Axpy(+opSign*vikjl, jkDesCreScr, ilCreDesComp(i, l));
                  btas::Axpy(+opSign*vikjl, jkCreDesScr, ilDesCreComp(i, l));
               }

               double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
               if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
               {
                  // Pil
                  btas::Axpy(+opSign*vijkl, jkCreCreScr, ilCreCreComp(i, l));
                  btas::Axpy(+opSign*vijkl, jkDesDesScr, ilDesDesComp(i, l));
                  // Qil
                  btas::Axpy(-opSign*vijkl, jkCreDesScr, ilCreDesComp(i, l));
                  btas::Axpy(-opSign*vijkl, jkDesCreScr, ilDesCreComp(i, l));
               }

               double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
               if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
               {
                  // Qil
                  btas::Axpy(+opSign*vijlk, jkCreDesScr, ilCreDesComp(i, l));
                  btas::Axpy(-opSign*vijlk, jkDesCreScr, ilCreDesComp(i, l));
                  btas::Axpy(+opSign*vijlk, jkDesCreScr, ilDesCreComp(i, l));
                  btas::Axpy(-opSign*vijlk, jkCreDesScr, ilDesCreComp(i, l));
               }
            }
         }
      }
   }

   for(size_t l = 0; l < nl; ++l)
   {
      size_t lx = l_loop_index[l];
      btas::SDArray<N> lCreCompScr;
      btas::SDArray<N> lDesCompScr;

      for(size_t i = 0; i < ni; ++i)
      {
         size_t ix = i_loop_index[i];
         // ClComp
         if(ilCreDesComp(i, l).size() > 0)
            btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0), ilCreDesComp(i, l), btas::shape(indexCi+2), 1.0, lCreCompScr);
         if(ilCreCreComp(i, l).size() > 0)
            btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(1), ilCreCreComp(i, l), btas::shape(indexCi+2), 1.0, lCreCompScr);

         // DlComp
         if(ilDesCreComp(i, l).size() > 0)
            btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(1), ilDesCreComp(i, l), btas::shape(indexCi+2), 1.0, lDesCompScr);
         if(ilDesDesComp(i, l).size() > 0)
            btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0), ilDesDesComp(i, l), btas::shape(indexCi+2), 1.0, lDesCompScr);
      }

      // Cl x DlComp
      if(lDesCompScr.size() > 0)
         btas::Contract( 1.0, blockDl->Cre(lx), btas::shape(0), lDesCompScr, btas::shape(indexDl), 1.0, hamScr);

      // Dl x ClComp
      if(lCreCompScr.size() > 0)
         btas::Contract( 1.0, blockDl->Cre(lx), btas::shape(1), lCreCompScr, btas::shape(indexDl), 1.0, hamScr);
   }

   btas::TVector<size_t, N> reorder;
   IndexPermute(reorder, indexDk, indexCj, indexCi, indexDl);
   PermutedAxpy( 1.0, hamScr, reorder, hamCiCjDkDl);

   TTNS_DEBUG("HamCiCjDkDlS : Finished.");
}

//
// H := Cl x [ Ci x [ vijkl [ Ck x [ Cj x bra ] ] ] ] ( all is in large block )
//
template<size_t N, size_t K, class Q>
void HamCiCjDkDlL (
      const btas::QSDArray<N, Q>& bra,
            Block* blockCi,
      const size_t& indexCi,
      const std::vector<size_t>& i_loop_index,
            Block* blockCj,
      const size_t& indexCj,
      const std::vector<size_t>& j_loop_index,
            Block* blockDk,
      const size_t& indexDk,
      const std::vector<size_t>& k_loop_index,
            Block* blockDl,
      const size_t& indexDl,
      const std::vector<size_t>& l_loop_index,
      const btas::IVector<K>& order,
            btas::QSDArray<N, Q>& hamCiCjDkDl)
{
   TTNS_DEBUG("HamCiCjDkDlL : Started.");

   // must be indexCi < indexCj < indexDk < indexDl
   {
      std::vector<size_t> originalOrder;
      originalOrder.push_back(indexCi);
      originalOrder.push_back(indexCj);
      originalOrder.push_back(indexDk);
      originalOrder.push_back(indexDl);

      std::vector<size_t> sortedOrder(originalOrder);
      std::sort(sortedOrder.begin(), sortedOrder.end());

      TTNS_ASSERT(std::equal(sortedOrder.begin(), sortedOrder.end(), originalOrder.begin()),
      "HamCiCjDkDlL : contraction order must be sorted before calling");
   }

   // parity operation
   btas::QSDArray<N, Q> braCopy(bra);
   std::vector<int> p1;
   {
      std::vector<size_t> bounds;
      bounds.push_back(order[indexCi]);
      bounds.push_back(order[indexCj]);
      bounds.push_back(order[indexDk]);
      bounds.push_back(order[indexDl]);
      std::sort(bounds.begin(), bounds.end());

      for(size_t iBlock = 0; iBlock < K; ++iBlock)
      {
         if(bounds[0] < order[iBlock] && order[iBlock] <= bounds[1])
            p1.push_back(iBlock);
         if(bounds[2] < order[iBlock] && order[iBlock] <= bounds[3])
            p1.push_back(iBlock);
      }
   }
   braCopy.parity(p1);

   btas::SDArray<N> hamScr;

   size_t ni = i_loop_index.size();
   size_t nj = j_loop_index.size();
   size_t nk = k_loop_index.size();
   size_t nl = l_loop_index.size();

   btas::TArray<btas::SDArray<N>, 2> ilCreCreComp(ni, nl);
   btas::TArray<btas::SDArray<N>, 2> ilCreDesComp(ni, nl);
   btas::TArray<btas::SDArray<N>, 2> ilDesCreComp(ni, nl);
   btas::TArray<btas::SDArray<N>, 2> ilDesDesComp(ni, nl);

   double opSign = 1.0;
   if(order[indexCi] > order[indexCj]) opSign = -opSign;
   if(order[indexCi] > order[indexDk]) opSign = -opSign;
   if(order[indexCi] > order[indexDl]) opSign = -opSign;
   if(order[indexCj] > order[indexDk]) opSign = -opSign;
   if(order[indexCj] > order[indexDl]) opSign = -opSign;
   if(order[indexDk] > order[indexDl]) opSign = -opSign;

   for(size_t j = 0; j < nj; ++j)
   {
      size_t jx = j_loop_index[j];
      int jSpin = 1.0; if(jx % 2 != 0) jSpin = -1.0;
      btas::SDArray<N> jCreScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0), braCopy, btas::shape(indexCj), 1.0, jCreScr);

      btas::SDArray<N> jDesScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(1), braCopy, btas::shape(indexCj), 1.0, jDesScr);

      for(size_t k = 0; k < nk; ++k)
      {
         size_t kx = k_loop_index[k];
         int kSpin = 1.0; if(kx % 2 != 0) kSpin = -1.0;
         btas::SDArray<N> jkCreCreScr;
         btas::Contract( 1.0, blockDk->Cre(kx), btas::shape(0), jCreScr, btas::shape(indexDk), 1.0, jkCreCreScr);

         btas::SDArray<N> jkCreDesScr;
         btas::Contract( 1.0, blockDk->Cre(kx), btas::shape(1), jCreScr, btas::shape(indexDk), 1.0, jkCreDesScr);

         btas::SDArray<N> jkDesCreScr;
         btas::Contract( 1.0, blockDk->Cre(kx), btas::shape(0), jDesScr, btas::shape(indexDk), 1.0, jkDesCreScr);

         btas::SDArray<N> jkDesDesScr;
         btas::Contract( 1.0, blockDk->Cre(kx), btas::shape(1), jDesScr, btas::shape(indexDk), 1.0, jkDesDesScr);

         for(size_t i = 0; i < ni; ++i)
         {
            size_t ix = i_loop_index[i];
            int iSpin = 1.0; if(ix % 2 != 0) iSpin = -1.0;
            for(size_t l = 0; l < nl; ++l)
            {
               size_t lx = l_loop_index[l];
               int lSpin = 1.0; if(lx % 2 != 0) lSpin = -1.0;

               double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
               if(fabs(vikjl) >= THRE_INT && (iSpin == jSpin) && (kSpin == lSpin))
               {
                  // Pil
                  btas::Axpy(-opSign*vikjl, jkCreCreScr, ilCreCreComp(i, l));
                  btas::Axpy(-opSign*vikjl, jkDesDesScr, ilDesDesComp(i, l));
                  // Qil
                  btas::Axpy(+opSign*vikjl, jkDesCreScr, ilCreDesComp(i, l));
                  btas::Axpy(+opSign*vikjl, jkCreDesScr, ilDesCreComp(i, l));
               }

               double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
               if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
               {
                  // Pil
                  btas::Axpy(+opSign*vijkl, jkCreCreScr, ilCreCreComp(i, l));
                  btas::Axpy(+opSign*vijkl, jkDesDesScr, ilDesDesComp(i, l));
                  // Qil
                  btas::Axpy(-opSign*vijkl, jkCreDesScr, ilCreDesComp(i, l));
                  btas::Axpy(-opSign*vijkl, jkDesCreScr, ilDesCreComp(i, l));
               }

               double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
               if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
               {
                  // Qil
                  btas::Axpy(+opSign*vijlk, jkCreDesScr, ilCreDesComp(i, l));
                  btas::Axpy(-opSign*vijlk, jkDesCreScr, ilCreDesComp(i, l));
                  btas::Axpy(+opSign*vijlk, jkDesCreScr, ilDesCreComp(i, l));
                  btas::Axpy(-opSign*vijlk, jkCreDesScr, ilDesCreComp(i, l));
               }
            }
         }
      }
   }

   for(size_t i = 0; i < ni; ++i)
   {
      size_t ix = i_loop_index[i];
      btas::SDArray<N> iCreCompScr;
      btas::SDArray<N> iDesCompScr;

      for(size_t l = 0; l < nl; ++l)
      {
         size_t lx = l_loop_index[l];
         // CiComp
         if(ilDesCreComp(i, l).size() > 0)
            btas::Contract( 1.0, blockDl->Cre(lx), btas::shape(0), ilDesCreComp(i, l), btas::shape(indexDl), 1.0, iCreCompScr);
         if(ilCreCreComp(i, l).size() > 0)
            btas::Contract( 1.0, blockDl->Cre(lx), btas::shape(1), ilCreCreComp(i, l), btas::shape(indexDl), 1.0, iCreCompScr);

         // DiComp
         if(ilCreDesComp(i, l).size() > 0)
            btas::Contract( 1.0, blockDl->Cre(lx), btas::shape(1), ilCreDesComp(i, l), btas::shape(indexDl), 1.0, iDesCompScr);
         if(ilDesDesComp(i, l).size() > 0)
            btas::Contract( 1.0, blockDl->Cre(lx), btas::shape(0), ilDesDesComp(i, l), btas::shape(indexDl), 1.0, iDesCompScr);
      }

      // Ci x DiComp
      if(iDesCompScr.size() > 0)
         btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0), iDesCompScr, btas::shape(indexCi+3), 1.0, hamScr);

      // Di x CiComp
      if(iCreCompScr.size() > 0)
         btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(1), iCreCompScr, btas::shape(indexCi+3), 1.0, hamScr);
   }

   btas::TVector<size_t, N> reorder;
   IndexPermute(reorder, indexCj, indexDk, indexDl, indexCi);
   PermutedAxpy( 1.0, hamScr, reorder, hamCiCjDkDl);

   TTNS_DEBUG("HamCiCjDkDlL : Finished.");
}

} // namespace ttns

#endif // __TTNS_OP_COMPONENTS_HAM_H

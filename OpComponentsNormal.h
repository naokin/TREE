#ifndef __TTNS_OP_COMPONENTS_NORMAL_H
#define __TTNS_OP_COMPONENTS_NORMAL_H

#include "OpComponentsUtil.h"
#include "ttns_assert.h"

namespace ttns
{

#ifndef THRE_INT
#define THRE_INT 1.0e-16
#endif

  //==========================================================================================//
  // Compute Normal Blocks for renormalization.                                        //
  //==========================================================================================//

  //
  // Cre Block
  //
  template<size_t N, class Q>
  void CiBlockCi(const btas::QSDArray<N, Q>& bra, const btas::QSDArray<N, Q>& ket,
                       Block* blockCi, const size_t& indexCi, const std::vector<size_t>& i_loop_index,
                 const btas::IVector<N >& order,
                       Block* newBlock)
  {
    const size_t K = N-1;

    btas::IVector<K> braContracts;
    braContracts[indexCi] = 0;
    for(size_t i = 0;           i < indexCi; ++i) braContracts[i] = i + 1;
    for(size_t i = indexCi + 1; i < K;  ++i) braContracts[i] = i;

    btas::IVector<K> ketContracts;
    for(size_t i = 0; i < K; ++i) ketContracts[i] = i;

    // parity operation
    btas::QSDArray<N, Q> braCopy(bra); 
    std::vector<int> p1;
    {
      size_t lowerBound = order[indexCi];
  
      for(size_t i = 0; i < K; ++i)
        if(lowerBound < order[i]) p1.push_back(i);
    }
    braCopy.parity(p1);
  
    for(size_t i = 0; i < i_loop_index.size(); ++i)
    {
      size_t ix = i_loop_index[i];

      // Ci
      btas::SDArray<N> iCreScr;
      btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0),
                            braCopy, btas::shape(indexCi), 1.0, iCreScr);
      btas::Contract( 1.0, iCreScr, braContracts, ket, ketContracts,
                       1.0, newBlock->Cre(ix));
    }
  }

  //
  // CreCre & CreDes Block
  //
  template<size_t N, class Q>
  void CijBlockCij(const btas::QSDArray<N, Q>& bra, const btas::QSDArray<N, Q>& ket,
                         Block* blockCij, const size_t& indexCij, const std::vector<size_t>& ij_loop_index,
                   const btas::IVector<N >& order,
                         Block* newBlock)
  {
    const size_t K = N-1;

    btas::IVector<K> braContracts;
    braContracts[indexCij] = 0;
    for(size_t i = 0;            i < indexCij; ++i) braContracts[i] = i + 1;
    for(size_t i = indexCij + 1; i < K;   ++i) braContracts[i] = i;

    btas::IVector<K> ketContracts;
    for(size_t i = 0; i < K; ++i) ketContracts[i] = i;

    for(size_t i = 0; i < ij_loop_index.size(); ++i)
    {
      size_t ix = ij_loop_index[i];
      // Bii
      btas::SDArray<N> iiCreDesScr;
      btas::Contract( 1.0, blockCij->CreDes(ix, ix), btas::shape(0),
                            bra, btas::shape(indexCij), 1.0, iiCreDesScr);
      btas::Contract( 1.0, iiCreDesScr, braContracts, ket, ketContracts,
                       1.0, newBlock->CreDes(ix, ix));

      for(size_t j = i + 1; j < ij_loop_index.size(); ++j)
      {
        size_t jx = ij_loop_index[j];
        // Aij
        btas::SDArray<N> ijCreCreScr;
        btas::Contract( 1.0, blockCij->CreCre(ix, jx), btas::shape(0),
                              bra, btas::shape(indexCij), 1.0, ijCreCreScr);
        btas::Contract( 1.0, ijCreCreScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreCre(ix, jx));
        ScaledCopy(-1.0, newBlock->CreCre(ix, jx),
                                newBlock->CreCre(jx, ix));
        // Bij
        btas::SDArray<N> ijCreDesScr;
        btas::Contract( 1.0, blockCij->CreDes(ix, jx), btas::shape(0),
                              bra, btas::shape(indexCij), 1.0, ijCreDesScr);
        btas::Contract( 1.0, ijCreDesScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreDes(ix, jx));
        btas::Permute(newBlock->CreDes(ix, jx), btas::shape(1, 0),
                       newBlock->CreDes(jx, ix));
      }
    }
  }

  template<size_t N, class Q>
  void CijBlockCiCjS(const btas::QSDArray<N, Q>& bra, const btas::QSDArray<N, Q>& ket,
                           Block* blockCi, const size_t& indexCi, const std::vector<size_t>& i_loop_index,
                           Block* blockCj, const size_t& indexCj, const std::vector<size_t>& j_loop_index,
                     const btas::IVector<N >& order,
                           Block* newBlock)
  {
    if(indexCi >= indexCj)
      TTNS_THROW(false, "CijBlockCiCjS; must be called with indexCi < indexCj");

    const size_t K = N-1;

    btas::IVector<K> braContracts;
    braContracts[indexCi] = 0;
    braContracts[indexCj] = 1;
    for(size_t i = 0;           i < indexCi; ++i) braContracts[i] = i + 2;
    for(size_t i = indexCi + 1; i < indexCj; ++i) braContracts[i] = i + 1;
    for(size_t i = indexCj + 1; i < K;  ++i) braContracts[i] = i;

    btas::IVector<K> ketContracts;
    for(size_t i = 0; i < K; ++i) ketContracts[i] = i;

    // parity operation
    btas::QSDArray<N, Q> braCopy(bra); 
    double opSign = 1.0; std::vector<int> p1;
    {
      size_t lowerBound = order[indexCi];
      size_t upperBound = order[indexCj];
      if(lowerBound > upperBound)
      {
        std::swap(lowerBound, upperBound);
        opSign = -opSign;
      }
  
      for(size_t i = 0; i < K; ++i)
        if(lowerBound < order[i] && order[i] <= upperBound)
          p1.push_back(i);
    }
    braCopy.parity(p1);

    for(size_t j = 0; j < j_loop_index.size(); ++j)
    {
      size_t jx = j_loop_index[j];

      btas::SDArray<N> jCreScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0),
                            braCopy, btas::shape(indexCj), 1.0, jCreScr);

      for(size_t i = 0; i < i_loop_index.size(); ++i)
      {
        size_t ix = i_loop_index[i];


        // Aij
        btas::SDArray<N> ijCreCreScr;
        // -opSign: CjCi -> CiCj
        btas::Contract(-opSign, blockCi->Cre(ix), btas::shape(0),
                                 jCreScr, btas::shape(indexCi+1), 1.0, ijCreCreScr);
        btas::Contract( 1.0, ijCreCreScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreCre(ix, jx));
        ScaledCopy(-1.0, newBlock->CreCre(ix, jx),
                                newBlock->CreCre(jx, ix));
        // Bij
        btas::SDArray<N> jiCreDesScr;
        // -opSign: CjCi -> CiCj
        btas::Contract(-opSign, blockCi->Cre(ix), btas::shape(1),
                                 jCreScr, btas::shape(indexCi+1), 1.0, jiCreDesScr);
        btas::Contract( 1.0, jiCreDesScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreDes(jx, ix));
        btas::Permute(newBlock->CreDes(jx, ix), btas::shape(1, 0),
                       newBlock->CreDes(ix, jx));
      }
    }
  }

  template<size_t N, class Q>
  void CijBlockCiCjM(const btas::QSDArray<N, Q>& bra, const btas::QSDArray<N, Q>& ket,
                           Block* blockCi, const size_t& indexCi, const std::vector<size_t>& i_loop_index,
                           Block* blockCj, const size_t& indexCj, const std::vector<size_t>& j_loop_index,
                     const btas::IVector<N >& order,
                           Block* newBlock)
  {
    if(indexCi >= indexCj)
      TTNS_THROW(false, "CijBlockCiCjL; must be called with indexCi < indexCj");

    const size_t K = N-1;

    btas::IVector<K> braContracts;
    braContracts[indexCi] = 0;
    braContracts[indexCj] = 1;
    for(size_t i = 0;           i < indexCi; ++i) braContracts[i] = i + 2;
    for(size_t i = indexCi + 1; i < indexCj; ++i) braContracts[i] = i + 1;
    for(size_t i = indexCj + 1; i < K;  ++i) braContracts[i] = i;

    btas::IVector<K> ketContracts;
    for(size_t i = 0; i < K; ++i) ketContracts[i] = i;

    // parity operation
    btas::QSDArray<N, Q> braCopy(bra); 
    double opSign = 1.0; std::vector<int> p1;
    {
      size_t lowerBound = order[indexCi];
      size_t upperBound = order[indexCj];
      if(lowerBound > upperBound)
      {
        std::swap(lowerBound, upperBound);
        opSign = -opSign;
      }
  
      for(size_t i = 0; i < K; ++i)
        if(lowerBound < order[i] && order[i] <= upperBound)
          p1.push_back(i);
    }
    braCopy.parity(p1);

    for(size_t j = 0; j < j_loop_index.size(); ++j)
    {
      size_t jx = j_loop_index[j];

      btas::SDArray<N> jCreScr;
      btas::Contract( 1.0, blockCj->Cre(jx), btas::shape(0),
                            braCopy, btas::shape(indexCj), 1.0, jCreScr);

      for(size_t i = 0; i < i_loop_index.size(); ++i)
      {
        size_t ix = i_loop_index[i];

        // Aij
        btas::SDArray<N> ijCreCreScr;
        btas::Contract(+opSign, blockCi->Cre(ix), btas::shape(0),
                                 jCreScr, btas::shape(indexCi+1), 1.0, ijCreCreScr);
        btas::Contract( 1.0, ijCreCreScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreCre(ix, jx));
        ScaledCopy(-1.0, newBlock->CreCre(ix, jx),
                                newBlock->CreCre(jx, ix));
        // Bij
        btas::SDArray<N> ijDesCreScr;
        btas::Contract(-opSign, blockCi->Cre(ix), btas::shape(1),
                                 jCreScr, btas::shape(indexCi+1), 1.0, ijDesCreScr);
        btas::Contract( 1.0, ijDesCreScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreDes(jx, ix));
        btas::Permute(newBlock->CreDes(jx, ix), btas::shape(1, 0),
                       newBlock->CreDes(ix, jx));
      }
    }
  }

  template<size_t N, class Q>
  void CijBlockCiCjL(const btas::QSDArray<N, Q>& bra, const btas::QSDArray<N, Q>& ket,
                           Block* blockCi, const size_t& indexCi, const std::vector<size_t>& i_loop_index,
                           Block* blockCj, const size_t& indexCj, const std::vector<size_t>& j_loop_index,
                     const btas::IVector<N >& order,
                           Block* newBlock)
  {
    if(indexCi >= indexCj)
      TTNS_THROW(false, "CijBlockCiCjL; must be called with indexCi < indexCj");

    const size_t K = N-1;

    btas::IVector<K> braContracts;
    braContracts[indexCi] = 1;
    braContracts[indexCj] = 0;
    for(size_t i = 0;           i < indexCi; ++i) braContracts[i] = i + 2;
    for(size_t i = indexCi + 1; i < indexCj; ++i) braContracts[i] = i + 1;
    for(size_t i = indexCj + 1; i < K;  ++i) braContracts[i] = i;

    btas::IVector<K> ketContracts;
    for(size_t i = 0; i < K; ++i) ketContracts[i] = i;

    // parity operation
    btas::QSDArray<N, Q> braCopy(bra); 
    double opSign = 1.0; std::vector<int> p1;
    {
      size_t lowerBound = order[indexCi];
      size_t upperBound = order[indexCj];
      if(lowerBound > upperBound)
      {
        std::swap(lowerBound, upperBound);
        opSign = -opSign;
      }
  
      for(size_t i = 0; i < K; ++i)
        if(lowerBound < order[i] && order[i] <= upperBound)
          p1.push_back(i);
    }
    braCopy.parity(p1);

    for(size_t i = 0; i < i_loop_index.size(); ++i)
    {
      size_t ix = i_loop_index[i];

      btas::SDArray<N> iCreScr;
      btas::Contract( 1.0, blockCi->Cre(ix), btas::shape(0),
                            braCopy, btas::shape(indexCi), 1.0, iCreScr);

      for(size_t j = 0; j < j_loop_index.size(); ++j)
      {
        size_t jx = j_loop_index[j];
        // Aij
        btas::SDArray<N> ijCreCreScr;
        btas::Contract(+opSign, blockCj->Cre(jx), btas::shape(0),
                                 iCreScr, btas::shape(indexCj), 1.0, ijCreCreScr);
        btas::Contract( 1.0, ijCreCreScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreCre(ix, jx));
        ScaledCopy(-1.0, newBlock->CreCre(ix, jx),
                                newBlock->CreCre(jx, ix));
        // Bij
        btas::SDArray<N> ijCreDesScr;
        btas::Contract(+opSign, blockCj->Cre(jx), btas::shape(1),
                                 iCreScr, btas::shape(indexCj), 1.0, ijCreDesScr);
        btas::Contract( 1.0, ijCreDesScr, braContracts, ket, ketContracts,
                         1.0, newBlock->CreDes(ix, jx));
        btas::Permute(newBlock->CreDes(ix, jx), btas::shape(1, 0),
                       newBlock->CreDes(jx, ix));
      }
    }
  }

} // namespace ttns

#endif // __TTNS_OP_COMPONENTS_NORMAL_H

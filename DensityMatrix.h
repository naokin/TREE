#ifndef __TTNS_DENSITY_MATRIX_H
#define __TTNS_DENSITY_MATRIX_H

#include <btas/DENSE/DArray.h>

#include "OpComponentsUtil.h"
#include "ttns_assert.h"

namespace ttns
{

btas::DArray<2> onepdm_alpha;
btas::DArray<2> onepdm_beta;

template<size_t Z, class Q>
void UpdateOnePDM(Site<Z, Q>* sys_dot)
{
   TTNS_DEBUG("UpdateOnePDM : Started.");

   Block* sysOps = sys_dot->Dot();

   if(!sysOps) return;

   const btas::QSDArray<Z+1, Q>& sysTns = sys_dot->Bra();

   size_t iphys = Z-1; if(sys_dot->Depth() == 0) iphys = Z;

   ComputeOnePDM(sysTns, sysOps, iphys);

   size_t n = sys_dot->Branches();

   for(size_t i = 0; i < n; ++i)
   {
      ComputeOnePDM(sysTns, sys_dot->Branch(i), i, sysOps, iphys);

      for(size_t j = i + 1; j < n; ++j)
      {
         ComputeOnePDM(sysTns, sys_dot->Branch(i), i, sys_dot->Branch(j), j);
      }
   }

   TTNS_DEBUG("UpdateOnePDM : Finished.");
}

template<size_t N, class Q>
void ComputeOnePDM (const btas::QSDArray<N, Q>& tns, const Block* block, const size_t& index)
{
   TTNS_DEBUG("ComputeOnePDM : Dii : Started.");

   std::vector<size_t> loop_index(block->inside());

   for(size_t i = 0; i < loop_index.size(); ++i)
   {
      size_t ix = loop_index[i];
      {
         btas::SDArray<N> DiiScr;
         btas::Contract(1.0, block->CreDes(ix, ix), btas::shape(0), tns, btas::shape(index), 1.0, DiiScr);

         btas::IVector<N> reorder;
         IndexPermute(reorder, index);
         btas::QSDArray<N, Q> Dii(tns.q(), tns.qshape());
         btas::Permute(DiiScr, reorder, Dii);

         size_t i2 = ix/2;
         if(ix % 2 == 0)
            onepdm_alpha(i2, i2) = btas::Dotc(Dii, tns);
         else
            onepdm_beta (i2, i2) = btas::Dotc(Dii, tns);
      }

      for(size_t j = i+1; j < loop_index.size(); ++j)
      {
         size_t jx = loop_index[j];

         if((ix % 2) != (jx % 2)) continue;

         btas::SDArray<N> DijScr;
         btas::Contract(1.0, block->CreDes(ix, jx), btas::shape(0), tns, btas::shape(index), 1.0, DijScr);

         btas::IVector<N> reorder;
         IndexPermute(reorder, index);
         btas::QSDArray<N, Q> Dij(tns.q(), tns.qshape());
         btas::Permute(DijScr, reorder, Dij);

         size_t i2 = ix/2;
         size_t j2 = jx/2;
         if(ix % 2 == 0)
         {
            onepdm_alpha(i2, j2) = btas::Dotc(Dij, tns);
            onepdm_alpha(j2, i2) = onepdm_alpha(i2, j2);
         }
         else
         {
            onepdm_beta (i2, j2) = btas::Dotc(Dij, tns);
            onepdm_beta (j2, i2) = onepdm_beta(i2, j2);
         }
      }
   }

   TTNS_DEBUG("ComputeOnePDM : Dii : Finished.");
}

template<size_t N, class Q>
void ComputeOnePDM(const btas::QSDArray<N, Q>& tns, const Block* block_i, const size_t& index_i, const Block* block_j, const size_t& index_j)
{
   TTNS_DEBUG("ComputeOnePDM : Dij : Started.");
   TTNS_ASSERT(index_i < index_j, "ComputeOnePDM; Must be index_i < index_j");

   std::vector<size_t> loop_index_i(block_i->inside());
   std::vector<size_t> loop_index_j(block_j->inside());
   for(size_t i = 0; i < loop_index_i.size(); ++i)
   {
      size_t ix = loop_index_i[i];
      for(size_t j = 0; j < loop_index_j.size(); ++j)
      {
         size_t jx = loop_index_j[j];
         if((ix % 2) != (jx % 2)) continue;

         btas::SDArray<N> iCre;
         btas::Contract( 1.0, block_i->Cre(ix), btas::shape(1), tns, btas::shape(index_i), 1.0, iCre);

         btas::SDArray<N> DijScr;
         btas::Contract( 1.0, block_j->Cre(jx), btas::shape(0), iCre, btas::shape(index_j), 1.0, DijScr);

         btas::IVector<N> reorder;
         IndexPermute(reorder, index_i, index_j);
         btas::QSDArray<N, Q> Dij(tns.q(), tns.qshape());
         btas::Permute(DijScr, reorder, Dij);

         size_t i2 = ix / 2;
         size_t j2 = jx / 2;
         if(ix % 2 == 0)
         {
            onepdm_alpha(i2, j2) = btas::Dotc(Dij, tns);
            onepdm_alpha(j2, i2) = onepdm_alpha(i2, j2);
         }
         else
         {
            onepdm_beta(i2, j2) = btas::Dotc(Dij, tns);
            onepdm_beta(j2, i2) = onepdm_beta(i2, j2);
         }
      }
   }

   TTNS_DEBUG("ComputeOnePDM : Dij : Finished.");
}

} // namespace ttns

#endif // __TTNS_DENSITY_MATRIX_H

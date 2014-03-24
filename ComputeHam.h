#ifndef __TTNS_COMPUTE_HAM_H
#define __TTNS_COMPUTE_HAM_H

#include "Block.h"
#include "OpComponentsHam.h"

namespace ttns
{

/// Compute H * C : local H
template<size_t K, size_t N, class Q>
void ComputeHamHam (
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::QSDArray<N, Q>& sgv)
{
   TTNS_DEBUG("ComputeHamHam : Started.");

   // V += Sum_{op1} H[op1] * 1
   for(long op1 = 0; op1 < K; ++op1)
   {
      if(!bs[op1]) continue;

      HamHam(psi, bs[op1], op1, sgv);
   }

   TTNS_DEBUG("ComputeHamHam : Finished.");
}

/// Compute H * C : Sum_{i} Ci * Ri
template<size_t K, size_t N, class Q>
void ComputeHamCiRi (
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::QSDArray<N, Q>& sgv)
{
   TTNS_DEBUG("ComputeHamCiRi : Started.");

   // Check # of small blocks, i.e. null, dot, and leaf sites
   long N_small_blocks = 0;
   for(size_t i = 0; i < K-1; ++i)
   {
      if(bs[i] && bs[i]->size() > 2) break; // theoretical
      ++N_small_blocks;
   }

   for(long op1 = 0; op1 < N_small_blocks; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         HamCiCjS(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, op2_loop_index, order, sgv);
      }
   }

   for(long op1 = N_small_blocks; op1 < static_cast<long>(K)-1; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         HamCiCjL(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, op2_loop_index, order, sgv);
      }
   }

   TTNS_DEBUG("ComputeHamCiRi : Finished.");
}

/// Compute H * C : Sum_{ij} Cij * Rij
template<size_t K, size_t N, class Q>
void ComputeHamCijRij (
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::QSDArray<N, Q>& sgv)
{
   TTNS_DEBUG("ComputeHamCijRij : Started.");

   // Check # of small blocks, i.e. null, dot, and leaf sites
   long N_small_blocks = 0;
   if(K > 2)
   {
      for(size_t i = 0; i < K-1; ++i)
      {
//       if(bs[i] && bs[i]->size() > 6) break; // theoretical
         if(bs[i] && bs[i]->size() > 2) break; // practical
         ++N_small_blocks;
      }
   }

   for(long op1 = 0; op1 < N_small_blocks; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         HamCijRijS(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, order, sgv);
      }
   }

   for(long op1 = N_small_blocks; op1 < static_cast<long>(K)-1; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < K; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         HamCijRijL(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, order, sgv);
      }
   }

   TTNS_DEBUG("ComputeHamCijRij : Finished.");
}

/********** 
template<size_t K, size_t N, class Q>
void ComputeHamCiCjRij (
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K>& order,
      const btas::QSDArray<K, Q>& psi,
            btas::QSDArray<K, Q>& sgv)
{
   // Check # of small blocks, i.e. null, dot, and leaf sites
   long N_small_blocks = 0;
   for(size_t i = 0; i < K-1; ++i)
   {
      if(bs[i] && bs[i]->size() > 2) break; // theoretical
      ++N_small_blocks;
   }

   //
   // A, B, C; small block
   //
   for(long op1 = 0; op1 < N_small_blocks-2; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < N_small_blocks-1; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = op2+1; op3 < N_small_blocks; ++op3)
         {
            if(!bs[op3]) continue;
            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Ci(a) x Cj(b) x Pij(c)
            HamCiCjRijS(psi, bs[op2], op2, op2_loop_index, bs[op1], op1, op1_loop_index, bs[op3], op3, order, sgv);

            // Ci(a) x Pij(b) x Cj(c)
            HamCiCjRijS(psi, bs[op3], op3, op3_loop_index, bs[op1], op1, op1_loop_index, bs[op2], op2, order, sgv);

            // Pij(a) x Ci(b) x Cj(c)
            HamCiCjRijS(psi, bs[op2], op2, op2_loop_index, bs[op3], op3, op3_loop_index, bs[op1], op1, order, sgv);
         }
      }
   }

   //
   // A, B; small block / C; large block
   //
   for(long op1 = 0; op1 < N_small_blocks-1; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < N_small_blocks; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = N_small_blocks; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;
            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Ci(a) x Cj(b) x Pij(c)
            HamCiCjRijS(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, op2_loop_index, bs[op3], op3, order, sgv);

            // Ci(a) x Pij(b) x Cj(c)
            HamCiCjRijS(psi, bs[op1], op1, op1_loop_index, bs[op3], op3, op3_loop_index, bs[op2], op2, order, sgv);

            // Pij(a) x Ci(b) x Cj(c)
            HamCiCjRijS(psi, bs[op2], op2, op2_loop_index, bs[op3], op3, op3_loop_index, bs[op1], op1, order, sgv);
         }
      }
   }

   //
   // A; small block / B, C; large block
   //
   for(long op1 = 0; op1 < N_small_blocks; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = N_small_blocks; op2 < static_cast<long>(K)-1; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = op2+1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;
            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Ci(a) x Cj(b) x Pij(c)
            HamCiCjRijS(psi, bs[op2], op2, op2_loop_index, bs[op1], op1, op1_loop_index, bs[op3], op3, order, sgv);

            // Ci(a) x Pij(b) x Cj(c)
            HamCiCjRijS(psi, bs[op3], op3, op3_loop_index, bs[op1], op1, op1_loop_index, bs[op2], op2, order, sgv);

            // Pij(a) x Ci(b) x Cj(c)
            HamCiCjRijS(psi, bs[op2], op2, op2_loop_index, bs[op3], op3, op3_loop_index, bs[op1], op1, order, sgv);
         }
      }
   }

   //
   // A, B, C; large block
   //
   for(long op1 = N_small_blocks; op1 < static_cast<long>(K)-2; ++op1)
   {
      if(!bs[op1]) continue;
      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < static_cast<long>(K)-1; ++op2)
      {
         if(!bs[op2]) continue;
         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = op2+1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;
            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Ci(a) x Cj(b) x Pij(c)
            HamCiCjRijL(psi, bs[op2], op2, op2_loop_index, bs[op1], op1, op1_loop_index, bs[op3], op3, order, sgv);

            // Ci(a) x Pij(b) x Cj(c)
            HamCiCjRijL(psi, bs[op3], op3, op3_loop_index, bs[op1], op1, op1_loop_index, bs[op2], op2, order, sgv);

            // Pij(a) x Ci(b) x Cj(c)
            HamCiCjRijL(psi, bs[op3], op3, op3_loop_index, bs[op2], op2, op2_loop_index, bs[op1], op1, order, sgv);
         }
      }
   }
}
**********/

/// Compute H * C : Sum_{ij} Ci * Cj * Rij
/// Algorithm with pre-contraction
template<size_t K, size_t N, class Q>
void ComputeHamCiCjRij (
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::QSDArray<N, Q>& sgv)
{
   TTNS_DEBUG("ComputeHamCiCjRij : Started.");

   // Check # of small blocks, i.e. null, dot, and leaf sites
   long N_small_blocks = 0;
   for(size_t i = 0; i < K-1; ++i)
   {
      if(bs[i] && bs[i]->size() > 2) break; // theoretical
      ++N_small_blocks;
   }

   //
   // A [ B, C ] & A [ C, B ]; A, B in small
   //
   for(long op1 = 0; op1 < N_small_blocks-1; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      btas::TArray<btas::SDArray<N>, 1> CreCompRi;
      btas::TArray<btas::SDArray<N>, 1> DesCompRi;

      for(long op2 = op1+1; op2 < N_small_blocks; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = op2+1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;

            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Cj(b) x Pij(c)
            CreCompCjRijS(psi, op1, op1_loop_index, bs[op2], op2, op2_loop_index, bs[op3], op3, order, CreCompRi, DesCompRi);

            // Cj(c) x Pij(b)
            CreCompCjRijS(psi, op1, op1_loop_index, bs[op3], op3, op3_loop_index, bs[op2], op2, order, CreCompRi, DesCompRi);
         }
      }

      HamCiRi(CreCompRi, DesCompRi, bs[op1], op1, op1_loop_index, order, sgv);
   }

   //
   // B [ C, A ]; A, B in small
   //
   for(long op2 = 1; op2 < N_small_blocks; ++op2)
   {
      if(!bs[op2]) continue;

      std::vector<size_t> op2_loop_index(bs[op2]->inside());

      btas::TArray<btas::SDArray<N>, 1> CreCompRi;
      btas::TArray<btas::SDArray<N>, 1> DesCompRi;

      for(long op1 = 0; op1 < op2; ++op1)
      {
         if(!bs[op1]) continue;

         std::vector<size_t> op1_loop_index(bs[op1]->inside());

         for(long op3 = op2+1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;

            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Cj(c) x Pij(a)
            CreCompCjRijS(psi, op2, op2_loop_index, bs[op3], op3, op3_loop_index, bs[op1], op1, order, CreCompRi, DesCompRi);
         }
      }

      HamCiRi(CreCompRi, DesCompRi, bs[op2], op2, op2_loop_index, order, sgv);
   }

   //
   // B [ A, C ] & B [ C, A ]; A in small
   //
   for(long op2 = N_small_blocks; op2 < static_cast<long>(K)-1; ++op2)
   {
      if(!bs[op2]) continue;

      std::vector<size_t> op2_loop_index(bs[op2]->inside());

      btas::TArray<btas::SDArray<N>, 1> CreCompRi;
      btas::TArray<btas::SDArray<N>, 1> DesCompRi;

      for(long op1 = 0; op1 < N_small_blocks; ++op1)
      {
         if(!bs[op1]) continue;

         std::vector<size_t> op1_loop_index(bs[op1]->inside());

         for(long op3 = op2+1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;

            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Cj(a) x Pij(c)
            CreCompCjRijS(psi, op2, op2_loop_index, bs[op1], op1, op1_loop_index, bs[op3], op3, order, CreCompRi, DesCompRi);

            // Cj(c) x Pij(a)
            CreCompCjRijS(psi, op2, op2_loop_index, bs[op3], op3, op3_loop_index, bs[op1], op1, order, CreCompRi, DesCompRi);
         }
      }

      HamCiRi(CreCompRi, DesCompRi, bs[op2], op2, op2_loop_index, order, sgv);
   }

   //
   // C [ A, B ]; A in small
   //
   for(long op3 = N_small_blocks+1; op3 < K; ++op3)
   {
      if(!bs[op3]) continue;

      std::vector<size_t> op3_loop_index(bs[op3]->inside());

      btas::TArray<btas::SDArray<N>, 1> CreCompRi;
      btas::TArray<btas::SDArray<N>, 1> DesCompRi;

      for(long op2 = N_small_blocks; op2 < op3; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op1 = 0; op1 < N_small_blocks; ++op1)
         {
            if(!bs[op1]) continue;

            std::vector<size_t> op1_loop_index(bs[op1]->inside());

            // Cj(a) x Pij(b)
            CreCompCjRijS(psi, op3, op3_loop_index, bs[op1], op1, op1_loop_index, bs[op2], op2, order, CreCompRi, DesCompRi);
         }
      }

      HamCiRi(CreCompRi, DesCompRi, bs[op3], op3, op3_loop_index, order, sgv);
   }

   //
   // B [ A, C ]; all in large
   //
   for(long op2 = N_small_blocks+1; op2 < static_cast<long>(K)-1; ++op2)
   {
      if(!bs[op2]) continue;

      std::vector<size_t> op2_loop_index(bs[op2]->inside());

      btas::TArray<btas::SDArray<N>, 1> CreCompRi;
      btas::TArray<btas::SDArray<N>, 1> DesCompRi;

      for(long op1 = N_small_blocks; op1 < op2; ++op1)
      {
         if(!bs[op1]) continue;

         std::vector<size_t> op1_loop_index(bs[op1]->inside());

         for(long op3 = op2+1; op3 < K; ++op3)
         {
            if(!bs[op3]) continue;

            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            // Cj(a) x Pij(c)
            CreCompCjRijS(psi, op2, op2_loop_index, bs[op1], op1, op1_loop_index, bs[op3], op3, order, CreCompRi, DesCompRi);
         }
      }

      HamCiRi(CreCompRi, DesCompRi, bs[op2], op2, op2_loop_index, order, sgv);
   }

   //
   // C [ A, B ] & C [ B, A ]; all in large
   //
   for(long op3 = N_small_blocks+2; op3 < K; ++op3)
   {
      if(!bs[op3]) continue;

      std::vector<size_t> op3_loop_index(bs[op3]->inside());

      btas::TArray<btas::SDArray<N>, 1> CreCompRi;
      btas::TArray<btas::SDArray<N>, 1> DesCompRi;

      for(long op2 = N_small_blocks+1; op2 < op3; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op1 = N_small_blocks; op1 < op2; ++op1)
         {
            if(!bs[op1]) continue;

            std::vector<size_t> op1_loop_index(bs[op1]->inside());

            // Cj(a) x Pij(b)
            CreCompCjRijS(psi, op3, op3_loop_index, bs[op1], op1, op1_loop_index, bs[op2], op2, order, CreCompRi, DesCompRi);

            // Cj(b) x Pij(a)
            CreCompCjRijS(psi, op3, op3_loop_index, bs[op2], op2, op2_loop_index, bs[op1], op1, order, CreCompRi, DesCompRi);
         }
      }

      HamCiRi(CreCompRi, DesCompRi, bs[op3], op3, op3_loop_index, order, sgv);
   }

   TTNS_DEBUG("ComputeHamCiCjRij : Finished.");
}

/// Compute H * C : Sum_{ij} Ci * Cj * Rij
template<size_t K, size_t N, class Q>
void ComputeHamCiCjDkDl (
      const btas::TVector<Block*, K>& bs,
      const btas::IVector<K>& order,
      const btas::QSDArray<N, Q>& psi,
            btas::QSDArray<N, Q>& sgv)
{
   TTNS_DEBUG("ComputeHamCiCjDkDl : Started.");

   // Check # of small blocks, i.e. null, dot, and leaf sites
   long N_small_blocks = 0;
   for(size_t i = 0; i < static_cast<long>(K)-2; ++i)
   {
      if(bs[i] && bs[i]->size() > 2) break; // theoretical
      ++N_small_blocks;
   }

   for(long op1 = 0; op1 < N_small_blocks-1; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < N_small_blocks; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = op2+1; op3 < static_cast<long>(K)-1; ++op3)
         {
            if(!bs[op3]) continue;

            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            for(long op4 = op3+1; op4 < K; ++op4)
            {
               if(!bs[op4]) continue;

               std::vector<size_t> op4_loop_index(bs[op4]->inside());

               HamCiCjDkDlS(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, op2_loop_index,
                                 bs[op3], op3, op3_loop_index, bs[op4], op4, op4_loop_index, order, sgv);
            }
         }
      }
   }

   N_small_blocks = std::max(N_small_blocks-1, 0l);

   for(long op1 = N_small_blocks; op1 < static_cast<long>(K)-3; ++op1)
   {
      if(!bs[op1]) continue;

      std::vector<size_t> op1_loop_index(bs[op1]->inside());

      for(long op2 = op1+1; op2 < static_cast<long>(K)-2; ++op2)
      {
         if(!bs[op2]) continue;

         std::vector<size_t> op2_loop_index(bs[op2]->inside());

         for(long op3 = op2+1; op3 < static_cast<long>(K)-1; ++op3)
         {
            if(!bs[op3]) continue;

            std::vector<size_t> op3_loop_index(bs[op3]->inside());

            for(long op4 = op3+1; op4 < K; ++op4)
            {
               if(!bs[op4]) continue;

               std::vector<size_t> op4_loop_index(bs[op4]->inside());

               HamCiCjDkDlL(psi, bs[op1], op1, op1_loop_index, bs[op2], op2, op2_loop_index,
                                 bs[op3], op3, op3_loop_index, bs[op4], op4, op4_loop_index, order, sgv);
            }
         }
      }
   }

   TTNS_DEBUG("ComputeHamCiCjDkDl : Finished.");
}

} // namespace ttns

#endif // __TTNS_COMPUTE_HAM_H

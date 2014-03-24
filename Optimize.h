#ifndef __TTNS_OPTIMIZE_H
#define __TTNS_OPTIMIZE_H

#include <iostream>
#include <iomanip>
#include <functional>

#include <time_stamp.h>

#include "Site.h"
#include "ComputeSigmaVector.h"
#include "ComputeDiagonals.h"
#include "CreateHalfDot.h"
#include "Davidson.h"
#include "Canonicalize.h"
#include "DensityMatrix.h"
#include "Randomize.h"

// Site<Z>* Site<Z>::next()
// -----------------------------------------------------------------
// increases bond index then returns pointer to next site
// in case of branch site,
//   if bond index > branches, set bond index to trunk (= -1)
//   if bond index is trunk, set bond index to first branch (= 0)
// in case of trunk site,
//   if bond index > branches, returns NULL pointer (set index = -1)

namespace ttns
{

template<size_t N, class Q>
double OneDotShared (
            btas::QSDArray<N, Q>& tns,
      const btas::TVector<Block*, N>& bs,
      const SweepParameters& param)
{
   TTNS_DEBUG("OneDotShared : Started.");

   btas::QSDArray<N, Q> hDiag(tns.q(), tns.qshape(), tns.dshape(), false);
   ComputeDiagonals(bs, hDiag);

   // contraction ordering
   btas::IVector<N> contractOrder;
   for(size_t i = 0; i < N; ++i) contractOrder[i] = i;

   davidson::Functor<N, Q> f_contract = std::bind(ComputeSigmaVector<N, Q>, bs, contractOrder, std::placeholders::_1, std::placeholders::_2);

   double energy = davidson::diagonalize(f_contract, hDiag, tns, param);

// Randomize(param.Noise(), tns);

   TTNS_DEBUG("OneDotShared : Finished.");

   return energy;
}

template<size_t NA, size_t NB, class Q>
double TwoDotShared (
            btas::QSDArray<NA, Q>& tnsA,
      const btas::TVector<Block*, NA-1>& bsA,
            btas::QSDArray<NB, Q>& tnsB,
      const btas::TVector<Block*, NB-1>& bsB,
      const size_t& bIndex,
      const size_t& bDepth,
      const SweepParameters& param,
      const Direction& direction)
{
   TTNS_DEBUG("TwoDotShared : Started.");

   btas::QSDArray<NA+NB-2, Q> psi;
   btas::Contract(1.0, tnsA, btas::shape(NA-1), tnsB, btas::shape(bIndex), 1.0, psi);

   btas::QSDArray<NA+NB-2, Q> hDiag(psi.q(), psi.qshape(), psi.dshape(), false);
   ComputeDiagonals(bsA, bsB, hDiag);

   btas::TVector<Block*, NA+NB-2> bs;
   for(size_t i = 0; i < NA-1; ++i) bs[i]      = bsA[i];
   for(size_t i = 0; i < NB-1; ++i) bs[NA-1+i] = bsB[i];

   // contraction ordering
   btas::IVector<NA+NB-2> ordering;
   for(size_t i = 0;      i < bIndex; ++i) ordering[NA-1+i] = i;
   for(size_t i = 0;      i < NA - 1; ++i) ordering[i]      = i + bIndex;
   for(size_t i = bIndex; i < NB - 1; ++i) ordering[NA-1+i] = i + NA - 1;

   davidson::Functor<NA+NB-2, Q> f_contract = std::bind(ComputeSigmaVector<NA+NB-2, Q>, bs, ordering, std::placeholders::_1, std::placeholders::_2);

   auto psi_dshape = psi.dshape();
   psi_dshape[NA-2] = btas::Dshapes(4, 1);
   if(bDepth == 0)
      psi_dshape[NA+NB-3] = btas::Dshapes(4, 1);
   else
      psi_dshape[NA+NB-4] = btas::Dshapes(4, 1);

   Randomize(param.Noise(), psi, psi_dshape);

   double energy = davidson::diagonalize(f_contract, hDiag, psi, param);

// Randomize(param.Noise(), psi);

   // canonicalization
   CanonicalizeTwoDot(psi, tnsA, tnsB, bIndex, param.M(), direction);

   TTNS_DEBUG("TwoDotShared : Finished.");

   return energy;
}

template<size_t NA, size_t NB, class Q>
double TwoDotShared (
            btas::QSDArray<NA+NB-2, Q>& psi,
            btas::QSDArray<NA, Q>& tnsA,
      const btas::TVector<Block*, NA-1>& bsA,
            btas::QSDArray<NB, Q>& tnsB,
      const btas::TVector<Block*, NB-1>& bsB,
      const size_t& bIndex,
      const size_t& bDepth,
      const SweepParameters& param,
      const Direction& direction)
{
   TTNS_DEBUG("TwoDotShared(with Guess) : Started.");

   // compute guess
   if(psi.size() == 0)
      btas::Contract(1.0, tnsA, btas::shape(NA-1), tnsB, btas::shape(bIndex), 1.0, psi);

   btas::QSDArray<NA+NB-2, Q> hDiag(psi.q(), psi.qshape(), psi.dshape(), false);
   ComputeDiagonals(bsA, bsB, hDiag);

   btas::TVector<Block*, NA+NB-2> bs;
   for(size_t i = 0; i < NA-1; ++i) bs[i]      = bsA[i];
   for(size_t i = 0; i < NB-1; ++i) bs[NA-1+i] = bsB[i];

   // contraction ordering
   btas::TVector<size_t, NA+NB-2> ordering;
   for(size_t i = 0;      i < bIndex; ++i) ordering[NA-1+i] = i;
   for(size_t i = 0;      i < NA - 1; ++i) ordering[i]      = i + bIndex;
   for(size_t i = bIndex; i < NB - 1; ++i) ordering[NA-1+i] = i + NA - 1;

   davidson::Functor<NA+NB-2, Q> f_contract = std::bind(ComputeSigmaVector<NA+NB-2, Q>, bs, ordering, std::placeholders::_1, std::placeholders::_2);

   auto psi_dshape = psi.dshape();
   psi_dshape[NA-2] = btas::Dshapes(4, 1);
   if(bDepth == 0)
      psi_dshape[NA+NB-3] = btas::Dshapes(4, 1);
   else
      psi_dshape[NA+NB-4] = btas::Dshapes(4, 1);

   Randomize(param.Noise(), psi, psi_dshape);

   double energy = davidson::diagonalize(f_contract, hDiag, psi, param);

// Randomize(param.Noise(), psi);

   // canonicalization
   CanonicalizeTwoDot(psi, tnsA, tnsB, bIndex, param.M(), direction);

   TTNS_DEBUG("TwoDotShared(with Guess) : Finished.");

   return energy;
}

template<size_t Z, class Q>
double Optimize (Site<Z, Q>* trunk, const SweepParameters& param, const std::string& prefix = ".")
{
   time_stamp ts;

   size_t iter = 0;
   bool converged = false;
   std::vector<double> sweepEnergies;
   double energySave = 0.0;
   Site<Z, Q>* sysDot = trunk;
   Site<Z, Q>* envDot = trunk->Current();
   if(!envDot) envDot = trunk->Next();

   ts.start();
   size_t istep = 0;
   while(iter < param.MaxIter())
   {
std::cout << "\t\t================================================================================" << std::endl;
std::cout << "\t\t\t\t\t    Sweep Iteration [ " << std::setw(4) << istep++ << " ] " << std::endl;
std::cout << "\t\t================================================================================" << std::endl;

std::cout << "\t\t\tLoading environment blocks...";
      envDot->LoadSite(prefix);
std::cout << "done" << std::endl;

      // swapping direction
      if(param.IsTwoDot() && envDot->Branches() == 0) std::swap(sysDot, envDot);

      double energy = integrals::moeri::core_energy;
      if(param.IsTwoDot())
         energy += TwoDot(sysDot, envDot, param);
      else
         energy += OneDot(sysDot, param);

std::cout << "\t\t\tSaving system block...";
      sysDot->SaveSite(prefix);
std::cout << "done" << std::endl;

      sweepEnergies.push_back(energy);

      sysDot = envDot;         /* next sys */
      envDot = sysDot->Next(); /* next env */

      if(!envDot) /* one sweep finished */
      {
         ++iter;
         istep = 0;

         if(!param.IsTwoDot()) UpdateOnePDM(sysDot);

         double minE = *std::min_element(sweepEnergies.begin(), sweepEnergies.end());
         double maxE = *std::max_element(sweepEnergies.begin(), sweepEnergies.end());
         double aveE =  std::accumulate (sweepEnergies.begin(), sweepEnergies.end(), 0.0)
                                                             / sweepEnergies.size();

         double residual = std::fabs(minE - energySave);
         energySave      = minE;

         converged = (iter > 1 && residual < param.ToleEnergy());

std::cout.setf(std::ios::fixed, std::ios::floatfield);
std::cout.precision(8);
std::cout << "\t\t################################################################################" << std::endl;
std::cout << "\t\t\tOptimization Cycle [ " << iter << " ] finished. " << std::endl;
std::cout << "\t\t\t\tEnergy: " << energySave << " ( Residual " << residual << " ) " << std::endl;
std::cout << "\t\t\t\tE(min): " << minE << ", E(max): "         << maxE              << std::endl;
std::cout << "\t\t\t\tE(ave): " << aveE << ", Error : "         << maxE - minE       << std::endl;
std::cout << "\t\t\t\tCpus  : " << std::setprecision(2) << ts.lap() << " sec. "      << std::endl;
if(converged)
std::cout << "\t\t\t\tStatus: converged" << std::endl;
else
std::cout << "\t\t\t\tStatus: not converged" << std::endl;
std::cout << "\t\t################################################################################" << std::endl;

         envDot = sysDot->Next(); /* next env */
         if(converged) break;
         sweepEnergies.clear();
      }
   }

std::cout.precision(2);
std::cout << "\t\t================================================================================" << std::endl;
if(converged)
std::cout << "\t\t\tOptimization converged after " << iter << " cycles. ( Cpu_Total: " << ts.elapsed() << " sec. ) " << std::endl;
else
std::cout << "\t\t\tWarning: not converged until " << iter << " cycles. ( Cpu_Total: " << ts.elapsed() << " sec. ) " << std::endl;
std::cout << "\t\t================================================================================" << std::endl;

   return energySave;
}

template<size_t ZS, class Q>
double OneDot (Site<ZS, Q>* sysDot, const SweepParameters& param)
{
   time_stamp ts;

std::cout << "\t\t--------------------------------------------------------------------------------" << std::endl;
std::cout << "\t\t\tsysDot : address   = " << sysDot << std::endl;
std::cout << "\t\t\t       : depth     = " << sysDot->Depth() << std::endl;
if(sysDot->Trunk())
{
std::cout << "\t\t\t       : trunk     = " << sysDot->Trunk()
                           << " ( # sites = " << sysDot->Trunk()->size() << " ) " << std::endl;
}
for(size_t i = 0; i < sysDot->Branches(); ++i)
{
std::cout << "\t\t\t       : branch[" << std::setw(1) << i << "] = " << sysDot->Branch(i)
                           << " ( # sites = " << sysDot->Branch(i)->size() << " ) " << std::endl;
}
std::cout << "\t\t--------------------------------------------------------------------------------" << std::endl;

   double energy = 0.0;

   btas::QSDArray<ZS+1, Q>& sysTns = sysDot->Bra();

   btas::TVector<Block*, ZS+1> sysBlocks;
   sysDot->LoadBlocks(sysBlocks);

   ts.start();
   energy = OneDotShared(sysTns, sysBlocks, param);

std::cout << "\t\t\tFinished onesite optimization ( Cpu: " << std::fixed << std::setprecision(2) << ts.lap() << " sec. ) " << std::endl;

   if(sysDot->Current() == sysDot->Trunk()) UpdateOnePDM(sysDot);

   sysDot->Canonicalize(param.M());
   sysDot->Renormalize();

std::cout << "\t\t\tFinished renormalization      ( Cpu: " << std::fixed << std::setprecision(2) << ts.lap() << " sec. ) " << std::endl;

   return energy;
}

template<size_t ZS, size_t ZE, class Q>
double TwoDot (Site<ZS, Q>* sysDot, Site<ZE, Q>* envDot, const SweepParameters& param)
{
   time_stamp ts;

std::cout << "\t\t--------------------------------------------------------------------------------" << std::endl;
std::cout << "\t\t\tsysDot : address   = " << sysDot << std::endl;
std::cout << "\t\t\t       : depth     = " << sysDot->Depth() << std::endl;
if(sysDot->Trunk())
{
std::cout << "\t\t\t       : trunk     = " << sysDot->Trunk()
                           << " ( # sites = " << sysDot->Trunk()->size() << " ) " << std::endl;
}
for(size_t i = 0; i < sysDot->Branches(); ++i)
{
std::cout << "\t\t\t       : branch[" << std::setw(1) << i << "] = " << sysDot->Branch(i)
                           << " ( # sites = " << sysDot->Branch(i)->size() << " ) " << std::endl;
}
std::cout << "\t\t--------------------------------------------------------------------------------" << std::endl;
std::cout << "\t\t\tenvDot : address   = " << envDot << std::endl;
std::cout << "\t\t\t       : depth     = " << envDot->Depth() << std::endl;
if(envDot->Trunk())
{
std::cout << "\t\t\t       : trunk     = " << envDot->Trunk()
                           << " ( # sites = " << envDot->Trunk()->size() << " ) " << std::endl;
}
for(size_t i = 0; i < envDot->Branches(); ++i)
{
std::cout << "\t\t\t       : branch[" << std::setw(1) << i << "] = " << envDot->Branch(i)
                           << " ( # sites = " << envDot->Branch(i)->size() << " ) " << std::endl;
}
std::cout << "\t\t--------------------------------------------------------------------------------" << std::endl;

   double energy = 0.0;

   // loading sysDot info
   btas::QSDArray<ZS+1, Q>& sysTns = sysDot->Bra();
   size_t sysIndex; /* bond to env */
   btas::TVector<Block*, ZS> sysBlocks;
   sysIndex = sysDot->LoadBlocks(sysBlocks);

   // loading envDot info
   btas::QSDArray<ZE+1, Q>& envTns = envDot->Bra();
   size_t envIndex; /* bond to sys */
   btas::TVector<Block*, ZE> envBlocks;
   envIndex = envDot->LoadBlocks(envBlocks);

   ts.start();

   bool doSysHalf = (ZS > 2 && param.IsSysHalf());
   bool doEnvHalf = (ZE > 2 && param.IsEnvHalf());

   if(sysDot->Depth() > envDot->Depth()) /* forward (leaves to trunk) */
   {

std::cout << "\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
std::cout << "\t\t                              FORWARD OPTIMIZATION                              " << std::endl;
std::cout << "\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

      btas::QSDArray<4, Q>  psi;

      btas::QSDArray<3, Q>  sysMps;
      btas::QSDArray<ZS, Q> sysRes;
      Block sysOpsHalf;
      btas::TVector<Block*, 2> sysBlocksHalf;

      if(doSysHalf)
      {
         btas::TVector<Block*, ZS-1> sysBlocksScr;
         for(size_t i = 0; i < ZS-1; ++i) sysBlocksScr[i] = sysBlocks[i];
         CreateHalfDot(sysTns, sysBlocksScr, sysDot->Depth(), sysIndex, sysMps, sysRes, &sysOpsHalf);
         sysBlocksHalf[0] = &sysOpsHalf;
         sysBlocksHalf[1] = sysBlocks[ZS-1];
      }

      btas::QSDArray<3, Q>  envMps;
      btas::QSDArray<ZE, Q> envRes;
      Block envOpsHalf;
      btas::TVector<Block*, 2> envBlocksHalf;

      if(doEnvHalf)
      {
         btas::TVector<Block*, ZE-1 > envBlocksScr;
         if(envDot->Depth() == 0)
         {
            for(size_t i = 0; i < ZE-1; ++i) envBlocksScr[i] = envBlocks[i];
         }
         else
         {
            for(size_t i = 0; i < ZE-2; ++i) envBlocksScr[i] = envBlocks[i];
            envBlocksScr[ZE-2] = envBlocks[ZE-1];
         }

         if(doEnvHalf && param.Noise() > 1.0e-12)
         {
            CreateRandHalfDot(envTns, envBlocksScr, envDot->Depth(), envIndex, sysMps, param, psi, envRes, &envOpsHalf);
         }
         else
         {
            // Load from previous sweep
            if(sysDot->Branches() > 0 && envDot->HalfBlock().size() > 0)
            {
               envMps     = envDot->HalfMps();
               envRes     = envDot->HalfTns();
               envOpsHalf = envDot->HalfBlock();
            }
            else
            {
               CreateHalfDot(envTns, envBlocksScr, envDot->Depth(), envIndex, envMps, envRes, &envOpsHalf);
            }
         }

         if(envDot->Depth() == 0)
         {
            envBlocksHalf[0] = &envOpsHalf;
            envBlocksHalf[1] = envBlocks[ZE-1];
         }
         else
         {
            envBlocksHalf[0] = envBlocks[ZE-2];
            envBlocksHalf[1] = &envOpsHalf;
         }
      }

      if(doSysHalf && doEnvHalf)
         energy = TwoDotShared(psi, sysMps, sysBlocksHalf, envMps, envBlocksHalf, 0, envDot->Depth(), param, Forward);
      if(doSysHalf && !doEnvHalf)
         energy = TwoDotShared(sysMps, sysBlocksHalf, envTns, envBlocks, envIndex, envDot->Depth(), param, Forward);
      if(!doSysHalf && doEnvHalf)
         energy = TwoDotShared(sysTns, sysBlocks, envMps, envBlocksHalf, 0, envDot->Depth(), param, Forward);
      if(!doSysHalf && !doEnvHalf)
         energy = TwoDotShared(sysTns, sysBlocks, envTns, envBlocks, envIndex, envDot->Depth(), param, Forward);

std::cout << "\t\t\tFinished optimization    ( Cpu: " << std::fixed << std::setprecision(2) << ts.lap() << " sec. ) " << std::endl;

      if(doSysHalf)
      {
         sysDot->HalfMps()   = sysMps;
         sysDot->HalfTns()   = sysRes;
         sysDot->HalfBlock() = sysOpsHalf;

         Blocking(sysMps, sysMps, sysBlocksHalf, 2, sysDot);
         BackTransform(sysMps, sysRes, sysDot->Depth(), sysIndex, sysTns);
      }
      else
      {
         Blocking(sysTns, sysTns, sysBlocks, ZS, sysDot);
      }
      if(doEnvHalf)
      {
         BackTransform(envMps, envRes, envDot->Depth(), envIndex, envTns);
      }
   }
   if(sysDot->Depth() < envDot->Depth()) /* backward(trunk to leaves) */
   {

std::cout << "\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
std::cout << "\t\t                              BACKWARD OPTIMIZATION                             " << std::endl;
std::cout << "\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

      btas::QSDArray<4, Q>  psi;

      btas::QSDArray<3, Q>  sysMps;
      btas::QSDArray<ZS, Q> sysRes;
      Block sysOpsHalf;
      btas::TVector<Block*, 2 > sysBlocksHalf;
      if(doSysHalf)
      {
         btas::TVector<Block*, ZS-1 > sysBlocksScr;
         if(sysDot->Depth() == 0)
         {
            for(size_t i = 0; i < ZS-1; ++i) sysBlocksScr[i] = sysBlocks[i];
         }
         else
         {
            for(size_t i = 0; i < ZS-2; ++i) sysBlocksScr[i] = sysBlocks[i];
            sysBlocksScr[ZS-2] = sysBlocks[ZS-1];
         }

         CreateHalfDot(sysTns, sysBlocksScr, sysDot->Depth(), sysIndex, sysMps, sysRes, &sysOpsHalf);

         if(sysDot->Depth() == 0)
         {
            sysBlocksHalf[0] = &sysOpsHalf;
            sysBlocksHalf[1] = sysBlocks[ZS-1];
         }
         else
         {
            sysBlocksHalf[0] = sysBlocks[ZS-2];
            sysBlocksHalf[1] = &sysOpsHalf;
         }
      }

      btas::QSDArray<3, Q>  envMps;
      btas::QSDArray<ZE, Q> envRes;
      Block envOpsHalf;
      btas::TVector<Block*, 2 > envBlocksHalf;
      if(doEnvHalf)
      {
         btas::TVector<Block*, ZE-1 > envBlocksScr;
         for(size_t i = 0; i < ZE-1; ++i) envBlocksScr[i] = envBlocks[i];

         if(doEnvHalf && param.Noise() > 1.0e-12)
         {
            CreateRandHalfDot(envTns, envBlocksScr, envDot->Depth(), envIndex, sysMps, param, psi, envRes, &envOpsHalf);
         }
         else
         {
            if(sysDot->Branches() > 0 && envDot->HalfBlock().size() > 0)
            {
               envMps     = envDot->HalfMps();
               envRes     = envDot->HalfTns();
               envOpsHalf = envDot->HalfBlock();
            }
            else
            {
               CreateHalfDot(envTns, envBlocksScr, envDot->Depth(), envIndex, envMps, envRes, &envOpsHalf);
            }
         }

         envBlocksHalf[0] = &envOpsHalf;
         envBlocksHalf[1] = envBlocks[ZE-1];
      }

      if(doSysHalf && doEnvHalf)
         energy = TwoDotShared(psi, envMps, envBlocksHalf, sysMps, sysBlocksHalf, 0, sysDot->Depth(), param, Backward);
      if(doSysHalf && !doEnvHalf)
         energy = TwoDotShared(envTns, envBlocks, sysMps, sysBlocksHalf, 0, sysDot->Depth(), param, Backward);
      if(!doSysHalf && doEnvHalf)
         energy = TwoDotShared(envMps, envBlocksHalf, sysTns, sysBlocks, sysIndex, sysDot->Depth(), param, Backward);
      if(!doSysHalf && !doEnvHalf)
         energy = TwoDotShared(envTns, envBlocks, sysTns, sysBlocks, sysIndex, sysDot->Depth(), param, Backward);

std::cout << "\t\t\tFinished optimization    ( Cpu: " << std::fixed << std::setprecision(2) << ts.lap() << " sec. ) " << std::endl;

      if(doSysHalf)
      {
         sysDot->HalfMps()   = sysMps;
         sysDot->HalfTns()   = sysRes;
         sysDot->HalfBlock() = sysOpsHalf;

         Blocking(sysMps, sysMps, sysBlocksHalf, 0, sysDot);
         BackTransform(sysMps, sysRes, sysDot->Depth(), sysIndex, sysTns);
      }
      else
      {
         Blocking(sysTns, sysTns, sysBlocks, sysIndex, sysDot);
      }
      if(doEnvHalf)
      {
         BackTransform(envMps, envRes, envDot->Depth(), envIndex, envTns);
      }
   }

   sysDot->Ket() = sysTns;
   envDot->Ket() = envTns;

std::cout << "\t\t\tFinished renormalization ( Cpu: " << std::fixed << std::setprecision(2) << ts.lap() << " sec. ) " << std::endl;

   return energy;
}

} // namespace ttns

#endif // __TTNS_OPTIMIZE_H

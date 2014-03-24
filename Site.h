#ifndef __TTNS_SITE_H
#define __TTNS_SITE_H

#include <vector>
#include <random>
#include <functional>

#include <sstream>
#include <fstream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "Block.h"
#include "Canonicalize.h"
#include "Renormalize.h"
#include "integrals.h"
#include "ttns_assert.h"

namespace ttns
{

#ifndef THRE_INT
#define THRE_INT 1.0e-16
#endif

enum SITE_TYPE { SPIN, FERMION, BOSON };

template<size_t Z, class Q>
class Site : public Block
{

private:

   //
   //  Member Variables
   //

   size_t m_depth; ///< Depth of Tree Graph

   int m_current; ///< Current index being renormalized, -1 for trunk

   int m_label; ///< Site label

   Site<Z, Q>* m_trunk; ///< pointer to my trunk site, if null, I'm root

   std::vector<Site<Z, Q>*> m_branch; ///< pointers to branch sites which conneted with me

   btas::Qshapes<Q> m_physIndex; ///< physical quanta index

   //
   //  Storages
   //  TNS indexing
   //  trunk : 0,1,2,... branches
   //        : Z         physical
   //  branch: 0,1,2,... branches
   //        : Z-1       physical
   //        : Z         trunk
   //

   btas::QSDArray<Z+1, Q> m_bra; ///< tns storage
   btas::QSDArray<Z+1, Q> m_ket; ///< ket = bra during opt.

   Block m_ops; ///< site operator

   Block m_opsHalf; ///< half-renormalized operator
   btas::QSDArray<3, Q> m_mpsHalf; ///< half-site mps
   btas::QSDArray<Z, Q> m_tnsHalf; ///< res tensor

   /// update TNS by res matrix
   void UpDate (const btas::QSDArray<2, Q>& res)
   {
      TTNS_DEBUG("Site::UpDate : Started.");

      btas::QSDArray<Z+1, Q> braScr;
      // called from trunk
      if(m_current < 0)
      {
         btas::Contract(1.0, m_bra, btas::shape(Z), res, btas::shape(1), 1.0, braScr);
         btas::Copy(braScr, m_bra);
      }
      // called from branch
      else
      {
         btas::Contract(1.0, res, btas::shape(1), m_bra, btas::shape(m_current), 1.0, braScr);
         btas::IVector<Z+1> reorder;
         reorder[m_current] = 0;
         for(size_t i = 0; i < m_current; ++i)     reorder[i] = i+1;
         for(size_t i = m_current+1; i < Z+1; ++i) reorder[i] = i;
         btas::Permute(braScr, reorder, m_bra);

         // Parity for Q_total
         if(m_bra.q().parity())
         {
            std::vector<int> p1;
            for(size_t i = 0; i < m_current; ++i) p1.push_back(i);
            m_bra.parity(p1);
         }
      }
      m_ket = m_bra;

      TTNS_DEBUG("Site::UpDate : Finished.");
   }

   // create quantum numbers
   btas::Qshapes<Q> CreateQIndices ()
   {
      // this must be called from branch site
      btas::Qshapes<Q> q0, qTrunk, qBranch;
      btas::TVector<btas::Qshapes<Q>, Z+1> q_shape;

      if(m_bra.qshape(Z).size() > 0)
      {
         qTrunk = -m_bra.qshape(Z);
      }
      else
      {
         q0.push_back(Q::zero());
         std::fill(q_shape.begin(), q_shape.end(), q0);

         // physical q_shape
         q_shape[Z-1] = m_physIndex;

         qTrunk = q0;
         for(size_t i = 0; i < m_branch.size(); ++i)
         {
            qBranch = m_branch[i]->CreateQIndices();
            qTrunk  = qTrunk & qBranch;
            q_shape[i] = qBranch;
         }
         qTrunk = qTrunk & q_shape[Z-1];

//
// has been moved to CreateRandomTns
//
//       size_t qSize   = qTrunk.size();
//       size_t qMxSize = 10; if(qSize % 2 == 1) --qMxSize;
//       if(qSize > qMxSize)
//       {
//         size_t offSet = (qSize - qMxSize) / 2;
//         qTrunk = btas::Qshapes<Q>(qTrunk.begin() + offSet, qTrunk.begin() + offSet + qMxSize);
//       }

         q_shape[Z] = -qTrunk;

         m_bra.resize(Q::zero(), q_shape);
      }

      return qTrunk;
   }

   std::string GetFileName (const std::string& prefix = ".")
   {
      std::ostringstream oss;

      oss << prefix << "/site.";
      if(m_current < 0)
         oss << m_label << "-tr.tmp";
      else
         oss << m_label << "-b" << m_current << ".tmp";

      return oss.str();
   }

public:

   Site ()
   :  m_depth (0), m_current (-1), m_label (0), m_trunk (nullptr)
   { }

  ~Site ()
   { }

   // Case 0: No Orbitals
   Site (Site* trunk, const int& lab)
   :  m_depth (0), m_current (-1), m_label (lab), m_trunk (nullptr)
   {
      if(trunk)
      {
         m_depth = trunk->Depth()+1;
         trunk->Connect(this);
      }

      m_physIndex = btas::Qshapes<Q>(1, Q::zero());
   }

   // Case 2: Fermion - 1 Spacial Orbital
   Site (Site* trunk, const int& lab, const size_t& iOrb, const size_t& jOrb)
   :  Block (), m_depth (0), m_current (-1), m_label (lab), m_trunk (nullptr)
   {
      if(trunk)
      {
         m_depth = trunk->Depth()+1;
         trunk->Connect(this);
      }

      std::vector<size_t> i_orbs = { iOrb, jOrb };

      this->resize(FERMION, i_orbs);
   }

   // *SITE_TYPE is not available ( only FERMION )
   void resize (const SITE_TYPE& type, const std::vector<size_t>& i_orbs)
   {
      // resizing
      m_ops.resize(i_orbs);

      TTNS_DEBUG("Site::resize : Address : " << this << " ( block ) / " << &m_ops << " ( site )");
      TTNS_DEBUG("Site::resize : Orbital : " << std::endl << m_ops);

      // quantum numbers of physical index ( FERMION )
      // TODO: make physical quanta for generic purpose
      btas::Qshapes<Q> q;

      q.push_back(Q( 0,  0));
      q.push_back(Q(+1, +1));
      q.push_back(Q(+1, -1));
      q.push_back(Q(+2,  0));

      m_physIndex = q;

      btas::TVector<btas::Qshapes<Q>, 2> qv = { q, -q };
      btas::DArray<2> block(1, 1);

      std::vector<size_t> loop_index(m_ops.inside());

      {
//       m_ops.Ham().resize(Q(0, 0), qv);
         m_ops.Ham().resize(btas::make_array(btas::Dshapes(4,1), btas::Dshapes(4,1)), false);

         size_t ix = loop_index[0]; // alpha
         size_t jx = loop_index[1]; // beta

         double tii = integrals::moeri::OneInt(ix, ix);
         block = tii;
         if(fabs(tii) >= THRE_INT) m_ops.Ham().insert(btas::shape(1, 1), block);

         double tjj = integrals::moeri::OneInt(jx, jx);
         block = tjj;
         if(fabs(tjj) >= THRE_INT) m_ops.Ham().insert(btas::shape(2, 2), block);

         double vij = tii + tjj + integrals::moeri::TwoInt(ix, jx, ix, jx);
//       double vij = integrals::hubbard::U;
         block = vij;
         if(fabs(vij) >= THRE_INT) m_ops.Ham().insert(btas::shape(3, 3), block);
      }
//    TTNS_DEBUG("Site::resize : m_ops.Ham() = " << m_ops.Ham());

      // normal operators
      for(size_t i = 0; i < loop_index.size(); ++i)
      {
         size_t ix = loop_index[i];
         if(ix % 2 == 0)
         {
//          m_ops.Cre(ix).resize(Q(1, 1), qv);
            m_ops.Cre(ix).resize(btas::make_array(btas::Dshapes(4,1), btas::Dshapes(4,1)), false);
            block = 1.0;
            m_ops.Cre(ix).insert(btas::shape(1, 0), block);
            block =-1.0;
            m_ops.Cre(ix).insert(btas::shape(3, 2), block); // parity
         }
         else
         {
//          m_ops.Cre(ix).resize(Q(1,-1), qv);
            m_ops.Cre(ix).resize(btas::make_array(btas::Dshapes(4,1), btas::Dshapes(4,1)), false);
            block = 1.0;
            m_ops.Cre(ix).insert(btas::shape(2, 0), block);
            m_ops.Cre(ix).insert(btas::shape(3, 1), block);
         }
      }

      for(size_t i = 0; i < loop_index.size(); ++i)
      {
         size_t ix = loop_index[i];
         // Bii
         btas::Contract( 1.0, m_ops.Cre(ix), btas::shape(1), m_ops.Cre(ix), btas::shape(1), 1.0, m_ops.CreDes(ix, ix));

         for(size_t j = i + 1; j < loop_index.size(); ++j)
         {
            size_t jx = loop_index[j];
            // Aij
            btas::Contract( 1.0, m_ops.Cre(ix), btas::shape(1), m_ops.Cre(jx), btas::shape(0), 1.0, m_ops.CreCre(ix, jx));
            m_ops.CreCre(jx, ix) = m_ops.CreCre(ix, jx);
            btas::Scal(-1.0, m_ops.CreCre(jx, ix));

            // Bij
            btas::Contract( 1.0, m_ops.Cre(ix), btas::shape(1), m_ops.Cre(jx), btas::shape(1), 1.0, m_ops.CreDes(ix, jx));
            btas::Permute(m_ops.CreDes(ix, jx), btas::shape(1, 0), m_ops.CreDes(jx, ix));
         }
      }

      // complementary operators
      std::vector<size_t> comp_index(m_ops.outside());
      for(size_t i = 0; i < comp_index.size(); ++i)
      {
         size_t ix = comp_index[i];
         int iSpin = 1; if(ix % 2 != 0) iSpin = -1;
         for(size_t j = 0; j < loop_index.size(); ++j)
         {
            size_t jx = loop_index[j];
            int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

            // 1-particle term: t_ij x Cj
            double hij = integrals::moeri::OneInt(ix, jx);
            if(fabs(hij) >= THRE_INT && (iSpin == jSpin))
               btas::Axpy(0.5*hij, m_ops.Cre(jx), m_ops.CreComp(ix));
         }
         // NOTE: loop order should be changed to j -> k -> l -> i ?
         for(size_t j = 0; j < loop_index.size(); ++j)
         {
            size_t jx = loop_index[j];
            int jSpin = 1; if(jx % 2 != 0) jSpin = -1;
            if(iSpin != jSpin) continue;

            for(size_t k = 0; k < loop_index.size(); ++k)
            {
               size_t kx = loop_index[k];
               int kSpin = 1; if(kx % 2 != 0) kSpin = -1;
               if(jx == kx) continue;

               for(size_t l = 0; l < loop_index.size(); ++l)
               {
                  size_t lx = loop_index[l];
                  int lSpin = 1; if(lx % 2 != 0) lSpin = -1;
                  if(kSpin != lSpin) continue;

                  double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
                  if(fabs(vikjl) < THRE_INT) continue;

                  // 2-particle term: v_ikjl x Cj x Ck x Dl
                  btas::Contract(vikjl, m_ops.Cre(jx), btas::shape(1), m_ops.CreDes(kx, lx), btas::shape(0), 1.0, m_ops.CreComp(ix));
               }
            }
         }
      }

      // Pij
      for(size_t k = 0; k < loop_index.size(); ++k)
      {
         size_t kx = loop_index[k];
         int kSpin = 1; if(kx % 2 != 0) kSpin = -1;
         for(size_t l = k + 1; l < loop_index.size(); ++l)
         {
            size_t lx = loop_index[l];
            int lSpin = 1; if(lx % 2 != 0) lSpin = -1;

            for(size_t i = 0; i < comp_index.size(); ++i)
            {
               size_t ix = comp_index[i];
               int iSpin = 1; if(ix % 2 != 0) iSpin = -1;
               for(size_t j = i + 1; j < comp_index.size(); ++j)
               {
                  size_t jx = comp_index[j];
                  int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

                  double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
                  if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
                     btas::Axpy(-vijkl, m_ops.CreCre(kx, lx), m_ops.CreCreComp(ix, jx));

                  double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
                  if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
                     btas::Axpy( vijlk, m_ops.CreCre(kx, lx), m_ops.CreCreComp(ix, jx));
               }
            }
         }
      }

      for(size_t i = 0; i < comp_index.size(); ++i)
      {
         size_t ix = comp_index[i];
         for(size_t j = i + 1; j < comp_index.size(); ++j)
         {
            size_t jx = comp_index[j];
            if(m_ops.CreCreComp(ix, jx).size() == 0) continue;

//          btas::Axpy(-1.0, m_ops.CreCreComp(ix, jx), m_ops.CreCreComp(jx, ix));
            m_ops.CreCreComp(jx, ix) = m_ops.CreCreComp(ix, jx);
            btas::Scal(-1.0, m_ops.CreCreComp(jx, ix));
         }
      }

      // Qij
      for(size_t k = 0; k < loop_index.size(); ++k)
      {
         size_t kx = loop_index[k];
         int kSpin = 1; if(kx % 2 != 0) kSpin = -1;
         for(size_t i = 0; i < comp_index.size(); ++i)
         {
            size_t ix = comp_index[i];
            int iSpin = 1; if(ix % 2 != 0) iSpin = -1;
            // Qii := (v_ikik - v_iikk) Ck x Dk
            double vikik = integrals::moeri::TwoInt(ix, kx, ix, kx);
            if(fabs(vikik) >= THRE_INT) // (iSpin == iSpin) && (kSpin == kSpin) is trivial
               btas::Axpy( vikik, m_ops.CreDes(kx, kx), m_ops.CreDesComp(ix, ix));

            double viikk = integrals::moeri::TwoInt(ix, ix, kx, kx);
            if(fabs(viikk) >= THRE_INT && (iSpin == kSpin))
               btas::Axpy(-viikk, m_ops.CreDes(kx, kx), m_ops.CreDesComp(ix, ix));

            for(size_t j = i + 1; j < comp_index.size(); ++j)
            {
               size_t jx = comp_index[j];
               int jSpin = 1; if(jx % 2 != 0) jSpin = -1;
               // Qij := (v_ikjk - v_ijkk) Ck x Dk
               double vikjk = integrals::moeri::TwoInt(ix, kx, jx, kx);
               if(fabs(vikjk) >= THRE_INT && (iSpin == jSpin)) // (kSpin == kSpin) is trivial
                  btas::Axpy( vikjk, m_ops.CreDes(kx, kx), m_ops.CreDesComp(ix, jx));

               double vijkk = integrals::moeri::TwoInt(ix, jx, kx, kx);
               if(fabs(vijkk) >= THRE_INT && (iSpin == kSpin) && (jSpin == kSpin))
                  btas::Axpy(-vijkk, m_ops.CreDes(kx, kx), m_ops.CreDesComp(ix, jx));
            }
         }

         for(size_t l = k + 1; l < loop_index.size(); ++l)
         {
            size_t lx = loop_index[l];
            int lSpin = 1; if(lx % 2 != 0) lSpin = -1;

            for(size_t i = 0; i < comp_index.size(); ++i)
            {
               size_t ix = comp_index[i];
               int iSpin = 1; if(ix % 2 != 0) iSpin = -1;
               // Qii := v_ikil Ck x Dl
               double vikil = integrals::moeri::TwoInt(ix, kx, ix, lx);
               if(fabs(vikil) >= THRE_INT && (kSpin == lSpin)) // (iSpin == iSpin) is trivial
               {
                  btas::Axpy( vikil, m_ops.CreDes(kx, lx), m_ops.CreDesComp(ix, ix));
                  btas::Axpy( vikil, m_ops.CreDes(lx, kx), m_ops.CreDesComp(ix, ix));
               }

               double viikl = integrals::moeri::TwoInt(ix, ix, kx, lx);
               if(fabs(viikl) >= THRE_INT && (iSpin == kSpin) && (iSpin == lSpin))
               {
                  btas::Axpy(-viikl, m_ops.CreDes(kx, lx), m_ops.CreDesComp(ix, ix));
                  btas::Axpy(-viikl, m_ops.CreDes(lx, kx), m_ops.CreDesComp(ix, ix));
               }

               for(size_t j = i + 1; j < comp_index.size(); ++j)
               {
                  size_t jx = comp_index[j];
                  int jSpin = 1; if(jx % 2 != 0) jSpin = -1;

                  double vikjl = integrals::moeri::TwoInt(ix, kx, jx, lx);
                  if(fabs(vikjl) >= THRE_INT && (iSpin == jSpin) && (kSpin == lSpin))
                  {
                     btas::Axpy( vikjl, m_ops.CreDes(kx, lx), m_ops.CreDesComp(ix, jx));
                     btas::Axpy( vikjl, m_ops.CreDes(lx, kx), m_ops.CreDesComp(ix, jx));
                  }

                  double vijkl = integrals::moeri::TwoInt(ix, jx, kx, lx);
                  if(fabs(vijkl) >= THRE_INT && (iSpin == kSpin) && (jSpin == lSpin))
                  {
                     btas::Axpy(-vijkl, m_ops.CreDes(lx, kx), m_ops.CreDesComp(ix, jx));
                  }

                  double vijlk = integrals::moeri::TwoInt(ix, jx, lx, kx);
                  if(fabs(vijlk) >= THRE_INT && (iSpin == lSpin) && (jSpin == kSpin))
                  {
                     btas::Axpy(-vijlk, m_ops.CreDes(kx, lx), m_ops.CreDesComp(ix, jx));
                  }
               }
            }
         }
      }

      for(size_t i = 0; i < comp_index.size(); ++i)
      {
         size_t ix = comp_index[i];
         for(size_t j = i + 1; j < comp_index.size(); ++j)
         {
            size_t jx = comp_index[j];
            if(m_ops.CreDesComp(ix, jx).size() == 0) continue;

            // Qji = (Qij)'
            btas::Permute(m_ops.CreDesComp(ix, jx), btas::shape(1, 0), m_ops.CreDesComp(jx, ix));
         }
      }

//    debug::PrintBlockOps(&m_ops, "CreCreComp");
   }

   const size_t& Depth () const { return m_depth; }

   size_t Branches () const { return m_branch.size(); }

   Site* Trunk (void) const { return m_trunk; }

   Site* Branch (size_t i) const { return m_branch[i]; }

   const btas::QSDArray<Z+1, Q>& Bra () const { return m_bra; }
         btas::QSDArray<Z+1, Q>& Bra ()       { return m_bra; }

   const btas::QSDArray<Z+1, Q>& Ket () const { return m_ket; }
         btas::QSDArray<Z+1, Q>& Ket ()       { return m_ket; }

   Block* Dot () { return (m_ops.size() == 0) ? nullptr : &m_ops; }

   const Block& HalfBlock(void) const { return m_opsHalf; }
         Block& HalfBlock(void)       { return m_opsHalf; }

   const btas::QSDArray<3, Q>& HalfMps(void) const { return m_mpsHalf; }
         btas::QSDArray<3, Q>& HalfMps(void)       { return m_mpsHalf; }

   const btas::QSDArray<Z, Q>& HalfTns(void) const { return m_tnsHalf; }
         btas::QSDArray<Z, Q>& HalfTns(void)       { return m_tnsHalf; }

   void LoadSite (const std::string& prefix = ".")
   {
      TTNS_DEBUG("Site::LoadSite : Started.");

      std::ifstream fst(GetFileName(prefix).c_str());
      boost::archive::binary_iarchive iar(fst);

      for(int i = 0; i < m_branch.size(); ++i)
         if(m_branch[i] && i != m_current) m_branch[i]->load(iar, 0);

      if(m_trunk && m_current >= 0)
         m_trunk->load(iar, 0);

      iar >> m_bra;
      iar >> m_ket;
      iar >> m_mpsHalf;
      iar >> m_tnsHalf;
      iar >> m_opsHalf;

      TTNS_DEBUG("Site::LoadSite : Finished.");
   }

   void SaveSite (const std::string& prefix = ".")
   {
      TTNS_DEBUG("Site::SaveSite : Started.");

      std::ofstream fst(GetFileName(prefix).c_str());
      boost::archive::binary_oarchive oar(fst);

      for(int i = 0; i < m_branch.size(); ++i)
         if(m_branch[i] && i != m_current) m_branch[i]->save(oar, 0);

      if(m_trunk && m_current >= 0)
         m_trunk->save(oar, 0);

      oar << m_bra;
      oar << m_ket;
      oar << m_mpsHalf;
      oar << m_tnsHalf;
      oar << m_opsHalf;

      m_bra.clear();
      m_ket.clear();
      m_mpsHalf.clear();
      m_tnsHalf.clear();
      m_opsHalf.clear();

      TTNS_DEBUG("Site::SaveSite : Finished.");
   }

   // loading blocks for onedot optimization
   void LoadBlocks(btas::TVector<Block*, Z+1>& blocks)
   {
      std::fill(blocks.begin(), blocks.end(), nullptr);

      for(size_t i = 0; i < m_branch.size(); ++i)
         blocks[i] = m_branch[i];

      if(m_depth == 0)
      {
         blocks[Z]   = this->Dot();
      }
      else
      {
         blocks[Z-1] = this->Dot();
         blocks[Z]   = m_trunk;
      }
   }

   // loading blocks for twodot optimization
   size_t LoadBlocks(btas::TVector<Block*, Z>& blocks)
   {
      size_t current;

      size_t nBlock = 0;
      std::fill(blocks.begin(), blocks.end(), nullptr);

      for(size_t i = 0; i < m_branch.size(); ++i)
         if(i != m_current) blocks[nBlock++] = m_branch[i];

      if(m_depth == 0)
      {
         current = static_cast<size_t>(m_current);
         blocks[Z-1] = this->Dot();
      }
      else
      {
         if(m_current == -1)
         {
            current = Z;
            blocks[Z-1] = this->Dot();
         }
         else
         {
            current = static_cast<size_t>(m_current);
            blocks[Z-2] = this->Dot();
            blocks[Z-1] = m_trunk;
         }
      }
      return current;
   }

   Site* Current () const
   {
      if(m_current == -1)
         return m_trunk;
      else if(m_current < m_branch.size())
         return m_branch[m_current];
      else
         return nullptr;
   }

   Site* Next ()
   {
      ++m_current;
      if(m_current < m_branch.size())
         return m_branch[m_current];

      m_current = -1;
      return m_trunk;
   }

   void Connect(Site* other)
   {
      if((m_depth + 1) == other->Depth()) // this is trunk
      {
         std::set<Site*> branchSet(m_branch.begin(), m_branch.end());
         if(branchSet.find(other) == branchSet.end())
         {
            m_branch.push_back(other);
            other->Connect(this);
         }
      }
      else if((m_depth - 1) == other->Depth()) // this is branch
      {
         if(m_trunk) m_trunk->DisConnect(this);
         m_trunk = other;
      }
      else
      {
         TTNS_ASSERT(false, "Site::Connect; cannot connect due to mismatched depths.");
      }
   }

   // disconnect all branches and trunk
   void DisConnect ()
   {
      while(m_branch.size() > 0) m_branch[0]->DisConnect(this);
   }

   // disconnect specified site
   void DisConnect (Site* other)
   {
      if((m_depth + 1) == other->Depth()) // this is trunk
      {
         for(auto it = m_branch.begin(); it != m_branch.end(); ++it)
         {
            if(*it == other)
            {
               m_branch.erase(it);
               break;
            }
         }
      }
      else if((m_depth - 1) == other->Depth()) // this is branch
      {
         if(m_trunk == other)
         {
            other->DisConnect(this);
            m_trunk = nullptr;
         }
      }
      else
      {
         TTNS_ASSERT(false, "Site::DisConnect; cannot disconnect due to mismatched depths.");
      }
   }

   //
   // initialization
   //
   void Initialize (const size_t& D, const std::string& prefix = ".")
   {
      TTNS_DEBUG("Site::Initialize : Started.");
//    std::cout << "\t\t\tinitialize ( " << this << " ) " << std::endl;

      // m_current is used to determine the branch calling update(...)
      m_current = -1;
      while(++m_current < m_branch.size()) m_branch[m_current]->Initialize(D, prefix);

      // initial canonicalization
      m_current = -1;
      if(m_trunk)
      {
         double norm = btas::Dotc(m_bra, m_bra);
         btas::Scal(1.0/sqrt(norm), m_bra);
         this->Canonicalize(D);
      }
      // normalize trunk tns
      else
      {
         double norm = btas::Dotc(m_bra, m_bra);
         btas::Scal(1.0/sqrt(norm), m_bra);
         m_ket = m_bra;
      }

      // storage initialization
      if(m_trunk)
      {
         this->Renormalize();
         this->SaveSite(prefix);
      }

      TTNS_DEBUG("Site::Initialize : Finished.");
   }

   //
   // canonicalization
   //
   void Canonicalize (const size_t& D)
   {
      TTNS_DEBUG("Site::Canonicalize : Started.");

      if(m_current == -1)
      {
         if(m_trunk)
         {
            btas::QSDArray<2, Q> res;
            CanonicalizeOneDot(m_bra, Z, D, res);
            m_trunk->UpDate(res);
         }
         else
         {
            TTNS_ASSERT(false, "failed Site::Canonicalize(const int&)");
         }
      }
      else
      {
         btas::QSDArray<2, Q> res;
         CanonicalizeOneDot(m_bra, m_current, D, res);
         m_branch[m_current]->UpDate(res);
      }

      m_ket = m_bra;

      TTNS_DEBUG("Site::Canonicalize : Finished.");
   }

   //
   // renormalization
   //
   void Renormalize ()
   {
      TTNS_DEBUG("Site::Renormalize : Started.");
//    std::cout << "\t\t\trenormalize ( " << this << " ) : " << m_current << std::endl;

      btas::TVector<Block*, Z> blocks;

      size_t nBlock = 0;
      std::fill(blocks.begin(), blocks.end(), nullptr);

      // branch index
      for(size_t i = 0; i < m_branch.size(); ++i)
         if(i != m_current) blocks[nBlock++] = m_branch[i];

      size_t current;
      if(m_depth == 0)
      {
         // physical index
         blocks[Z-1] = this->Dot();
         current = static_cast<size_t>(m_current);
      }
      else
      {
         if(m_current == -1)
         {
            current = Z;
            blocks[Z-1]  = this->Dot();
         }
         else
         {
            current = static_cast<size_t>(m_current);
            blocks[Z-2]  = this->Dot();
            blocks[Z-1]  = m_trunk;
         }
      }

      Blocking(m_bra, m_ket, blocks, current, this);

      TTNS_DEBUG("Site::Renormalize : Finished.");
   }

   //
   // creating random tns
   //

   // for trunk site
   void CreateRandomTns (const Q& qTotal, const int& D)
   {
      if(m_depth != 0)
         TTNS_ASSERT(false, "Site::CreateRandomTns ( for TRUNK ) is called from branch");

      btas::Qshapes<Q> q0;
      btas::TVector<btas::Qshapes<Q>, Z+1> q_shape;
      btas::TVector<btas::Dshapes,    Z+1> d_shape;

      q0.push_back(Q::zero());

      std::fill(q_shape.begin(), q_shape.end(), q0);

      q_shape[Z] = m_physIndex;
      d_shape[Z] = btas::Dshapes(q_shape[Z].size(), 1);

      // create quantum q_shape
      size_t nBranches = m_branch.size();
      for(size_t i = 0; i < nBranches; ++i)
         q_shape[i] = m_branch[i]->CreateQIndices();

      // reducing quantum q_shape which doesn't contribute
      Q lhs;
      size_t iCount, nCount;
      btas::IVector<Z+1> index;
      for(size_t i = 0; i < nBranches; ++i)
      {
         for(auto it = q_shape[i].begin(); it != q_shape[i].end();)
         {
            nCount = 1;
            for(size_t j = 0; j < Z+1; ++j) if(j != i) nCount *= q_shape[j].size();

            index = btas::uniform<size_t, Z+1>(0ul);
            iCount = 0;
            while(iCount < nCount)
            {
               Q lhs = *it;
               for(size_t j = 0; j < Z+1; ++j) if(j != i) lhs = lhs + q_shape[j][index[j]];
               if(lhs == qTotal) break;

               for(size_t j = Z; j >= 0; --j)
               {
                  if(j != i)
                  {
                     if(++index[j] < q_shape[j].size()) break;
                     index[j] = 0;
                  }
               }
               ++iCount;
            }

            if(iCount == nCount) q_shape[i].erase(it);
            else                 it++;
         }
      }

      // further truncation
//    for(size_t i = 0; i < nBranches; ++i)
//    {
//       size_t qSize   = q_shape[i].size();
//       size_t qMxSize = 10; if(qSize % 2 == 1) --qMxSize;
//       if(qSize > qMxSize)
//       {
//          size_t offSet = (qSize - qMxSize) / 2;
//          q_shape[i] = btas::Qshapes<Q>(q_shape[i].begin() + offSet, q_shape[i].begin() + offSet + qMxSize);
//       }
//    }

      // resizing
      // branch
      for(size_t i = 0; i < nBranches; ++i)
         d_shape[i] = btas::Dshapes(q_shape[i].size(), D);
      // dummy
      for(size_t i = nBranches; i < Z; ++i)
         d_shape[i] = btas::Dshapes(q_shape[i].size(), 1);

      std::mt19937 rgen;
      std::uniform_real_distribution<double> dist(-1.0, 1.0);

      m_bra.resize(qTotal, q_shape, d_shape, std::bind(dist, rgen));
//    m_bra.resize(qTotal, q_shape, d_shape, 1.0);

      for(size_t i = 0; i < nBranches; ++i)
         m_branch[i]->CreateRandomTns(q_shape[i], D);
   }

   // for branch site
   void CreateRandomTns (const btas::Qshapes<Q>& qTrunk, const int& D)
   {
      if(m_depth == 0)
         TTNS_ASSERT(false, "Site::CreateRandomTns ( for BRANCH ) is called from trunk");

      btas::TVector<btas::Qshapes<Q>, Z+1> q_shape;
      btas::TVector<btas::Dshapes,    Z+1> d_shape;

      q_shape      = m_bra.qshape();
      d_shape[Z-1] = btas::Dshapes(q_shape[Z-1].size(), 1);

      size_t nBranches = m_branch.size();

      // reducing quantum q_shape which doesn't contribute
//    if(m_bra.qshape()[Z].size() != qTrunk.size()) // sometimes, gives error
      {
         q_shape[Z] = -qTrunk;

         Q lhs;
         size_t iCount, nCount;
         btas::IVector<Z+1> index;
         for(size_t i = 0; i < nBranches; ++i)
         {
            for(auto it = q_shape[i].begin(); it != q_shape[i].end();)
            {
               nCount = 1;
               for(size_t j = 0; j < Z; ++j) if(j != i) nCount *= q_shape[j].size();

               index = btas::uniform<size_t, Z+1>(0ul);
               iCount = 0;
               while(iCount < nCount)
               {
                  lhs = *it;
                  for(size_t j = 0; j < Z; ++j) if(j != i) lhs = lhs + q_shape[j][index[j]];

                  bool found = false;
                  for(size_t k = 0; k < qTrunk.size(); ++k)
                  {
                     if(lhs == qTrunk[k])
                     {
                        found = true;
                        break;
                     }
                  }
                  if(found) break;

                  for(size_t j = Z-1; j >= 0; --j)
                  {
                     if(j != i)
                     {
                        if(++index[j] < q_shape[j].size()) break;
                        index[j] = 0;
                     }
                  }
                  ++iCount;
               }

               if(iCount == nCount)
                  q_shape[i].erase(it);
               else
                  it++;
            }
         }
      }

      // further truncation
      for(size_t i = 0; i < nBranches; ++i)
      {
         size_t qSize   = q_shape[i].size();
         size_t qMxSize = 20; if(qSize % 2 == 1) --qMxSize;
         if(qSize > qMxSize)
         {
            size_t offSet = (qSize - qMxSize) / 2;
            q_shape[i] = btas::Qshapes<Q>(q_shape[i].begin()+offSet, q_shape[i].begin()+offSet+qMxSize);
         }
      }

      // resizing
      // branch
      for(size_t i = 0; i < nBranches; ++i)
         d_shape[i] = btas::Dshapes(q_shape[i].size(), D);
      // dummy
      for(size_t i = nBranches; i < Z-1; ++i)
         d_shape[i] = btas::Dshapes(q_shape[i].size(), 1);
      // trunk
         d_shape[Z] = btas::Dshapes(q_shape[Z].size(), D);

      std::mt19937 rgen;
      std::uniform_real_distribution<double> dist(-1.0, 1.0);

      m_bra.resize(Q::zero(), q_shape, d_shape, std::bind(dist, rgen));
//    m_bra.resize(Q::zero(), q_shape, d_shape, 1.0);

      for(size_t i = 0; i < nBranches; ++i)
         m_branch[i]->CreateRandomTns(q_shape[i], D);
   }

   //
   // HF config. as initial
   //

   Q CreateHFConf (const std::map<Site*, Q>& occ)
   {
      btas::TVector<btas::Qshapes<Q>, Z+1> q_shape;
      btas::TVector<btas::Dshapes,    Z+1> d_shape;

      size_t nBranches = m_branch.size();

      Q qTotal = Q::zero();

      for(size_t i = 0; i < nBranches; ++i)
      {
         Q qBranch;
         qBranch = m_branch[i]->CreateHFConf(occ);

//       size_t ne = qBranch.Particles();
//       int ns = qBranch.Spins();
         btas::Qshapes<Q> iq;
//       iq.push_back(Q(ne-2, ns));
         iq.push_back(qBranch);
//       iq.push_back(Q(ne+2, ns));

         q_shape[i] = iq;
         qTotal = qTotal + qBranch;
      }

      qTotal = qTotal + occ.find(this)->second;

      if(m_depth == 0)
      {
         for(size_t i = nBranches; i < Z; ++i)
            q_shape[i] = btas::Qshapes<Q>(1, Q::zero());

         q_shape[Z]   = m_physIndex;
      }
      else
      {
         for(size_t i = nBranches; i < Z-1; ++i)
            q_shape[i] = btas::Qshapes<Q>(1, Q::zero());

//       size_t ne = qTotal.Particles();
//       int ns = qTotal.Spins();
         btas::Qshapes<Q> iq;
//       iq.push_back(Q(ne-2, ns));
         iq.push_back(qTotal);
//       iq.push_back(Q(ne+2, ns));

         q_shape[Z-1] = m_physIndex;
         q_shape[Z]   = -iq;
      }

      for(size_t iRank = 0; iRank < Z+1; ++iRank)
         d_shape[iRank] = btas::Dshapes(q_shape[iRank].size(), 1);

      if(m_depth == 0)
         m_bra.resize(qTotal, q_shape, d_shape, 1.0);
      else
         m_bra.resize(Q::zero(), q_shape, d_shape, 1.0);

      return qTotal;
   }

}; // class Site

} // namespace ttns

#endif // __TTNS_SITE_H

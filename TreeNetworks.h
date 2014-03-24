#ifndef __TTNS_TREE_NETWORKS_H
#define __TTNS_TREE_NETWORKS_H

#include "Input.h"
#include "Site.h"
#include "Optimize.h"
#include "Fermion.h"

namespace ttns
{

//
// total system which consists of site objects
//
template<size_t Z>
class TreeNetworks
{

private:

   std::vector<int> m_connect; ///< connectivity

   std::vector<Site<Z, Fermion>*> m_sites; ///< site objects

   double m_energy; ///< energy

   bool m_isInit; ///< whether to initialize tree

   // copying function
   void DeepCopy (const TreeNetworks& other)
   {
      m_connect = other.m_connect;
      m_energy  = other.m_energy;
      m_isInit  = other.m_isInit;

      this->Create(m_connect);

      for(size_t iSite = 0; iSite < m_sites.size(); ++iSite) *m_sites[iSite] = *other.m_sites[iSite];
   }

   // destruct
   void Clean ()
   {
      m_energy = 0.0;
      m_isInit = false;

      for(size_t iSite = 0; iSite < m_sites.size(); ++iSite) delete m_sites[iSite];

      m_sites.clear();
   }

public:

   TreeNetworks ()
   :  m_energy (0.0), m_isInit (false)
   { }

  ~TreeNetworks ()
   {
      this->Clean();
   }

   TreeNetworks (const std::vector<int>& connect)
   {
      this->Create(connect);
   }

   TreeNetworks (const Input& input)
   {
      this->Compute(input);
   }

   TreeNetworks (const TreeNetworks& other)
   {
      this->DeepCopy(other);
   }

   // energy
   double Energy () const { return m_energy; }

   // sites
   size_t Size () const { return m_sites.size(); }
   Site<Z, Fermion>* operator() (const size_t& iSite) const { return m_sites[iSite]; }
   Site<Z, Fermion>* operator[] (const size_t& iSite) const { return m_sites[iSite]; }

   TreeNetworks& operator= (const TreeNetworks& other)
   {
      this->DeepCopy(other);
      return *this;
   }

   void Create (const std::vector<int>& connect, const std::vector<size_t>& orbitals)
   {
      Block::orbitals() = std::accumulate(orbitals.begin(), orbitals.end(), 0);

      this->Clean();
      m_connect = connect;

      size_t nSpinOrb = 0;

      TTNS_DEBUG("TreeNetworks::Create : Create Site[0]");

      if(orbitals[0] > 0)
         m_sites.push_back(new Site<Z, Fermion>(NULL, 0, nSpinOrb++, nSpinOrb++));
      else
         m_sites.push_back(new Site<Z, Fermion>(NULL, 0));

      for(size_t iSite = 1; iSite < connect.size(); ++iSite)
      {
         TTNS_DEBUG("TreeNetworks::Create : Create Site[" << iSite << "]");

         if(orbitals[iSite] > 0)
            m_sites.push_back(new Site<Z, Fermion>(m_sites[connect[iSite]], iSite, nSpinOrb++, nSpinOrb++));
         else
            m_sites.push_back(new Site<Z, Fermion>(m_sites[connect[iSite]], iSite));
      }
   }

   void Initialize (const Fermion& qTotal, const Input& input)
   {
//    m_sites[0]->CreateRandomTns(qTotal, input[0].M());

      std::vector<int> guessConfig(input.GuessConfig());
      if(guessConfig.size() > 0)
      {
         size_t nSites = m_sites.size();
         size_t nElectrons = qTotal.p();
         int nSpins = qTotal.s();

         size_t nDoublyOcc = (nElectrons - nSpins) / 2;
         size_t nOcc = nDoublyOcc + nSpins;

         std::map< size_t, size_t > orbitalMap(input.OrbitalMap());
         std::vector< size_t > offSet;
         {
            size_t nDummy = 0;
            std::vector< size_t > orbitals(input.Orbitals());
            for(size_t i = 0; i < nSites; ++i)
            {
               if(orbitals[i] > 0)
                  offSet.push_back(nDummy);
               else
                  nDummy++;
            }
         }

         btas::Qshapes<Fermion> hfConfig(nSites, Fermion::zero());

         if(guessConfig.size() == nSites)
         {
            for(size_t i = 0; i < guessConfig.size(); ++i)
            {
               size_t iSite = orbitalMap.find(i)->second;
               switch(guessConfig[i])
               {
                  case  1:
                     hfConfig[iSite+offSet[iSite]] = Fermion(1, 1);
                     break;
                  case -1:
                     hfConfig[iSite+offSet[iSite]] = Fermion(1,-1);
                     break;
                  case  2:
                     hfConfig[iSite+offSet[iSite]] = Fermion(2, 0);
                     break;
               }
            }
         }
         else
         {
            for(size_t i = 0; i < nDoublyOcc; ++i)
            {
               size_t iSite = orbitalMap.find(i)->second;
               hfConfig[iSite+offSet[iSite]] = Fermion(2, 0);
            }
            for(size_t i = nDoublyOcc; i < nOcc; ++i)
            {
               size_t iSite = orbitalMap.find(i)->second;
               hfConfig[iSite+offSet[iSite]] = Fermion(1, 1);
            }
         }

         std::map< Site<Z, Fermion>*, Fermion > occ;
         for(size_t i = 0; i < hfConfig.size(); ++i)
            occ.insert(std::make_pair(m_sites[i], hfConfig[i]));
         m_sites[0]->CreateHFConf(occ);
      }
      else
      {
         m_sites[0]->CreateRandomTns(qTotal, 1);
      }

      m_sites[0]->Initialize(0, input.prefix());
      m_isInit = true;
   }

   TreeNetworks& Optimize(const Input& input)
   {
      if(!m_isInit)
      {
         this->Create(input.Connection(), input.Orbitals());
         this->Initialize(input.QTotal(), input);
      }

      for(size_t iTask = 0; iTask < input.Schedules(); ++iTask)
      {
         std::cout << "\t\t************************* Task [ " << std::setw(3) << iTask << " ] *************************" << std::endl;
         std::cout << input[iTask] << std::endl;
         m_energy = ttns::Optimize(m_sites[0], input[iTask], input.prefix());
         // printing results here
      }

      return *this;
   }

   TreeNetworks& Compute(const Input& input)
   {
      TTNS_DEBUG("TreeNetworks::Compute : Calling TreeNetworks::Create(...)");
      this->Create(input.Connection(), input.Orbitals());
      TTNS_DEBUG("TreeNetworks::Compute : Calling TreeNetworks::Initialize(...)");
      this->Initialize(input.QTotal(), input);
      TTNS_DEBUG("TreeNetworks::Compute : Calling TreeNetworks::Optimize(...)");
      return this->Optimize(input);
   }
};

} // namespace ttns

#endif // __TTNS_TREE_NETWORKS_H

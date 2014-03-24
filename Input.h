#ifndef __TTNS_INPUT_H
#define __TTNS_INPUT_H

#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include "Fermion.h"
#include "integrals.h"

namespace ttns
{

class SweepParameters
{
//private:
public:
   size_t m_mxRedSt; ///< max number of reduced states
   size_t m_mxItOpt; ///< max iterations for sweep
   size_t m_mxItDav; ///< max iterations for Davidson
   size_t m_mxNDav;  ///< max dimension of Davidson
   double m_toleE;   ///< tolerance of energy for sweep
   double m_toleD;   ///< tolerance of Davidson's eigensolver
   double m_noise;   ///< scaling factor of random noise
   bool   m_twoDot;  ///< flag of twodot algorithm
   bool   m_sysHalf; ///< flag of half-dot decomposition for system-dot
   bool   m_envHalf; ///< flag of half-dot decomposition for environment-dot
public:
   SweepParameters(void)
   {
      m_mxRedSt = 1;
      m_mxItOpt = 100;
      m_mxItDav = 20;
      m_mxNDav  = 20;
      m_toleE   = 1.0e-7;
      m_toleD   = 1.0e-8;
      m_noise   = 0.0;
      m_twoDot  = false;
      m_sysHalf = false;
      m_envHalf = false;
   }

   SweepParameters(const SweepParameters& other)
   {
      m_mxRedSt = other.m_mxRedSt;
      m_mxItOpt = other.m_mxItOpt;
      m_mxItDav = other.m_mxItDav;
      m_mxNDav  = other.m_mxNDav;
      m_toleE   = other.m_toleE;
      m_toleD   = other.m_toleD;
      m_noise   = other.m_noise;
      m_twoDot  = other.m_twoDot;
      m_sysHalf = other.m_sysHalf;
      m_envHalf = other.m_envHalf;
   }

  ~SweepParameters(void) { }

   size_t M(void)               const { return m_mxRedSt; }
   size_t MaxIter(void)         const { return m_mxItOpt; }
   size_t MaxIterDavidson(void) const { return m_mxItDav; }
   size_t MaxDimDavidson(void)  const { return m_mxNDav;  }
   double ToleEnergy(void)      const { return m_toleE;   }
   double ToleDavidson(void)    const { return m_toleD;   }
   double Noise(void)           const { return m_noise;   }
   bool   IsTwoDot(void)        const { return m_twoDot;  }
   bool   IsSysHalf(void)       const { return m_sysHalf; }
   bool   IsEnvHalf(void)       const { return m_envHalf; }

// friend std::ostream& ::operator<< (std::ostream& ost, const SweepParameters& param);

   friend class Input;
};

class Input
{
//private:
public:

   size_t                         m_rank;       // maximum bonds on one-site
   std::vector< SweepParameters > m_schedule;   // sets of parameters
   std::vector< int >             m_connect;    // connectivity
   std::vector< size_t >          m_orbitals;   // #s of orbitals inside site
   std::vector< int >             m_guessConf;  // guess configuration
   std::map< size_t, size_t >     m_orbitalMap; // orbital ordering to be used integral sort
   Fermion                        m_qTotal;     // total quantum #
   std::string                    m_int1;       // one-electron integral file
   std::string                    m_int2;       // two-electron integral file
   std::string                    m_fcidump;    // molpro FCIDUMP fule
   std::string                    m_prefix;     // prefix for scratch files
   bool                           m_usefcidump;

public:
   Input(void)
   {
      m_rank = 2;
      m_int1 = "qcdmrg.int1";
      m_int2 = "qcdmrg.int2";
      m_fcidump = "FCIDUMP";
      m_prefix = ".";
      m_usefcidump = true;
      m_qTotal = Fermion(0, 0);
   }
  ~Input(void) { }

   Input(const std::string& config)
   {
      m_rank = 2;
      m_int1 = "qcdmrg.int1";
      m_int2 = "qcdmrg.int2";
      m_fcidump = "FCIDUMP";
      m_prefix = ".";
      m_usefcidump = true;
      m_qTotal = Fermion(0, 0);

      std::ifstream fconf(config.c_str());
      ReadConfigFile(fconf);
   }

   Input(std::ifstream& fconf)
   {
      m_rank = 2;
      m_int1 = "qcdmrg.int1";
      m_int2 = "qcdmrg.int2";
      m_fcidump = "FCIDUMP";
      m_prefix = ".";
      m_usefcidump = true;
      m_qTotal = Fermion(0, 0);

      ReadConfigFile(fconf);
   }

   Input(const Input& input)
   {
      m_rank       = input.m_rank;
      m_int1       = input.m_int1;
      m_int2       = input.m_int2;
      m_fcidump    = input.m_fcidump;
      m_prefix     = input.m_prefix;
      m_qTotal     = input.m_qTotal;
      m_schedule   = input.m_schedule;
      m_connect    = input.m_connect;
      m_orbitalMap = input.m_orbitalMap;
      m_qTotal     = input.m_qTotal;
   }

   void ReadConfigFile(const std::string& config)
   {
      std::ifstream fconf(config.c_str());
      ReadConfigFile(fconf);
   }

   // read parameters from config file
   void ReadConfigFile(std::ifstream& fconf)
   {
      size_t nAlpha, nBeta;

      std::string entry;
      while(fconf >> entry)
      {
        if(entry == "graph")
        {
            BTAS_DEBUG("Input::ReadConfigFile : Read graph keyword");
          m_connect.clear();

          size_t nSites; fconf >> nSites;
          size_t nOrb = 0;
          std::vector< size_t > nBonds(nSites, 0);
          for(size_t iSite = 0; iSite < nSites; ++iSite)
          {
            int orbLabel; fconf >> orbLabel;
            if(orbLabel < 0)
            {
              m_orbitals.push_back(0);
            }
            else
            {
              m_orbitals.push_back(2);
              m_orbitalMap.insert(std::make_pair(orbLabel, nOrb++));
            }
            // trunk : connect = -1
            // branch: connect = trunk site to be connected
            int connect;  fconf >> connect;
            m_connect.push_back(connect);
            if(connect >= 0)
            {
              ++nBonds[iSite];
              ++nBonds[connect];
            }
          }
          m_rank = *std::max_element(nBonds.begin(), nBonds.end());
        }

        if(entry == "alpha") fconf >> nAlpha;
        if(entry == "beta" ) fconf >> nBeta;

        if(entry == "hubbard") fconf >> integrals::hubbard::U;

        if(entry == "guess")
        {
            BTAS_DEBUG("Input::ReadConfigFile : Read guess keyword");
          std::string guessConf;
          fconf >> guessConf;
          for(std::string::iterator it = guessConf.begin(); it != guessConf.end(); ++it)
          {
            char occ = *it;
            if(occ == '0')
              m_guessConf.push_back( 0);
            if(occ == 'u')
              m_guessConf.push_back( 1);
            if(occ == 'd')
              m_guessConf.push_back(-1);
            if(occ == '2')
              m_guessConf.push_back( 2);
          }
        }

        if(entry == "schedule")
        {
            BTAS_DEBUG("Input::ReadConfigFile : Read schedule keyword");
          size_t nTasks; fconf >> nTasks;

          // M, IsTwoDot, IsSysHalf, IsEnvHalf, ToleE, Noise
          m_schedule = std::vector<SweepParameters>(nTasks, SweepParameters());
          for(size_t iTask = 0; iTask < nTasks; ++iTask)
          {
            fconf >> m_schedule[iTask].m_mxRedSt;
            fconf >> m_schedule[iTask].m_mxItOpt;

            std::string algo;

            fconf >> algo;
            if(algo == "twodot") m_schedule[iTask].m_twoDot = true;
            else                 m_schedule[iTask].m_twoDot = false;

            fconf >> algo;
            if(algo == "h" || algo == "half") m_schedule[iTask].m_sysHalf = true;
            else                              m_schedule[iTask].m_sysHalf = false;
            fconf >> algo;
            if(algo == "h" || algo == "half") m_schedule[iTask].m_envHalf = true;
            else                              m_schedule[iTask].m_envHalf = false;

            fconf >> m_schedule[iTask].m_toleE;
//          fconf >> m_schedule[iTask].m_toleD;
            fconf >> m_schedule[iTask].m_noise;
                     m_schedule[iTask].m_toleD = m_schedule[iTask].m_toleE / 10;
          }
        }

        if(entry == "fcidump")
        {
          m_usefcidump = true;
          fconf >> m_fcidump;
        }
        if(entry == "prefix")
         {
            fconf >> m_prefix;
         }
        if(entry == "int1")
        {
          m_usefcidump = false;
          fconf >> m_int1;
        }
        if(entry == "int2")
        {
          m_usefcidump = false;
          fconf >> m_int2;
        }
      }

      m_qTotal = Fermion(nAlpha + nBeta, nAlpha - nBeta);

      return;
   }

   // rank
   size_t Rank(void) const { return m_rank; }

   // schedules
   size_t            Schedules  (void)             const { return m_schedule.size(); }
   SweepParameters   operator() (const size_t& iTask) const { return m_schedule[iTask]; }
   SweepParameters   operator[] (const size_t& iTask) const { return m_schedule[iTask]; }

   // connectivity
   std::vector< int > Connection(void) const { return m_connect; }

   // #s of orbitals
   std::vector< size_t > Orbitals(void) const { return m_orbitals; }

   // orbital map
   std::vector< int > GuessConfig(void) const { return m_guessConf; }

   // orbital map
   std::map< size_t, size_t > OrbitalMap(void) const { return m_orbitalMap; }

   // q total
   Fermion QTotal(void) const { return m_qTotal; }

   // integral files
   const char* FileInt1(void) const { return m_int1.c_str(); }
   const char* FileInt2(void) const { return m_int2.c_str(); }
   const char* FciDump(void) const { return m_fcidump.c_str(); }
   inline bool UseFciDump(void) const { return m_usefcidump; }
   const std::string& prefix () const { return m_prefix; }

   // printing
// friend std::ostream& ::operator<< (std::ostream& ost, const Input& input);
};

} // namespace ttns

inline std::ostream& operator<< (std::ostream& ost, const ttns::SweepParameters& param)
{
ost.precision(2);
ost << "\t\t\tSweep Parameters: " << std::endl;
ost << "\t\t\tGeneral Information " << std::endl;
ost << "\t\t\t\tMax States  = " << param.m_mxRedSt << std::endl;
ost << "\t\t\t\tNoise       = " << std::scientific << param.m_noise << std::endl;
ost << "\t\t\t\tConvergence = " << std::scientific << param.m_toleE << std::endl;
ost << "\t\t\t\tIterations  = " << param.m_mxItOpt << std::endl;
ost << "\t\t\tDavidson Eigensolver " << std::endl;
ost << "\t\t\t\tRitz Dims.  = " << param.m_mxNDav << std::endl;
ost << "\t\t\t\tConvergence = " << std::scientific << param.m_toleD << std::endl;
ost << "\t\t\t\tIterations  = " << param.m_mxItDav << std::endl;
ost << "\t\t\tAlgorithm " << std::endl;
if(param.m_twoDot)  ost << "\t\t\t\tTwo-Dot";
else                ost << "\t\t\t\tOne-Dot";
if(param.m_sysHalf) ost << " / Half-Site ( Sys: On";
else                ost << " / Half-Site ( Sys: Off";
if(param.m_envHalf) ost << ", Env: On )";
else                ost << ", Env: Off )";
ost << std::endl;
return ost;
}

inline std::ostream& operator<< (std::ostream& ost, const ttns::Input& input)
{
ost << "\t================================================================================================" << std::endl;
ost << "\t\t\t\t\t\tSITE INFORMATION " << std::endl;
ost << "\t------------------------------------------------------------------------------------------------" << std::endl;
ost << "\t\t\t# sites : " << input.m_connect.size() << " / Global quanta = " << input.m_qTotal << std::endl;
for(size_t iTask = 0; iTask < input.Schedules(); ++iTask)
{
ost << "\t\t************************* Task [ " << std::setw(3) << iTask << " ] *************************" << std::endl;
ost << input.m_schedule[iTask] << std::endl;
}
return ost;
}

#endif

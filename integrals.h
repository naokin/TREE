#ifndef __TTNS_INTEGRALS_H
#define __TTNS_INTEGRALS_H

#include <vector>
#include <fstream>
#include <map>
#include <boost/algorithm/string.hpp>

#include <btas/DENSE/DArray.h>

std::vector<std::string> gettoken(std::ifstream& fin)
{
  // read msg from fin
  std::string msg;
  std::getline(fin, msg);
  boost::trim(msg);
  // split msg into token
  std::vector<std::string> tok;
  boost::split(tok, msg, boost::is_any_of("=, \t"), boost::token_compress_on);
  return tok;
}

namespace integrals
{
  namespace heisenberg
  {
    double J  = 1.0;
    double Jz = 1.0;
  };

  namespace hubbard
  {
    double t  = 1.0;
    double U  = 1.0;
  };

  namespace moeri
  {
    double core_energy = 0.0;
    btas::DArray<2> oneint;
    btas::DArray<4> twoint;

    void ReadFromFciDump(std::ifstream& fdump, const std::map<size_t, size_t>& OrbMap)
    {
      std::vector<std::string> tok;

      size_t norbs = 0;
      // get first line
      tok = gettoken(fdump);
      for(size_t i = 0; i < tok.size(); ++i) {
        if(tok[i] == "NORB" ) norbs = atoi(tok[++i].c_str());
      }

      oneint.resize(norbs, norbs);
      twoint.resize(norbs, norbs, norbs, norbs);

      // find &END control
      while(1) {
        tok = gettoken(fdump);
        if(tok[0] == "&END" || tok[0] == "/") break;
      }

      // read integrals
      while((tok = gettoken(fdump)).size() > 1)
      {
        double value = atof(tok[0].c_str());
        int p = atoi(tok[1].c_str()) - 1;
        int q = atoi(tok[2].c_str()) - 1;
        int r = atoi(tok[3].c_str()) - 1;
        int s = atoi(tok[4].c_str()) - 1;
        if(p <  0 && q <  0) {
          core_energy = value;
        }
        else if(r <  0 && s <  0) {
          size_t i = OrbMap.find(p)->second;
          size_t j = OrbMap.find(q)->second;
          oneint(i, j) = value;
          oneint(j, i) = value;
        }
        else {
          size_t i = OrbMap.find(p)->second;
          size_t j = OrbMap.find(q)->second;
          size_t k = OrbMap.find(r)->second;
          size_t l = OrbMap.find(s)->second;
          twoint(i, k, j, l) = value;
          twoint(i, l, j, k) = value;
          twoint(j, k, i, l) = value;
          twoint(j, l, i, k) = value;
          twoint(k, i, l, j) = value;
          twoint(k, j, l, i) = value;
          twoint(l, i, k, j) = value;
          twoint(l, j, k, i) = value;
        }
      }

      return;
    }

    void ReadOneInt(std::ifstream& fin, const std::map<size_t, size_t>& OrbMap)
    {
      size_t nSites;
      size_t p,q;
      size_t i,j;
      double value;

      fin >> nSites;
      oneint.resize(nSites, nSites);
      oneint = 0.0;

      while(!fin.eof())
      {
        fin >> p >> q >> value;
        i = OrbMap.find(p)->second;
        j = OrbMap.find(q)->second;
        oneint(i, j) = value;
        oneint(j, i) = value;
      }

      return;
    }

    void ReadTwoInt(std::ifstream& fin, const std::map<size_t, size_t>& OrbMap)
    {
      size_t nSites;
      size_t p,q,r,s;
      size_t i,j,k,l;
      double value;

      fin >> nSites;
      twoint.resize(nSites, nSites, nSites, nSites);
      twoint = 0.0;

      while(!fin.eof())
      {
        fin >> p >> q >> r >> s >> value;
        i = OrbMap.find(p)->second;
        j = OrbMap.find(q)->second;
        k = OrbMap.find(r)->second;
        l = OrbMap.find(s)->second;

        twoint(i,j,k,l) = value;
        twoint(i,l,k,j) = value;
        twoint(j,i,l,k) = value;
        twoint(j,k,l,i) = value;
        twoint(k,j,i,l) = value;
        twoint(k,l,i,j) = value;
        twoint(l,i,j,k) = value;
        twoint(l,k,j,i) = value;
      }

      return;
    }

    void DumpOneInt(std::ofstream& fdump, const size_t& nSites)
    {
      fdump << nSites << std::endl;
      fdump.setf(std::ios::scientific, std::ios::floatfield); fdump.precision(20);
      for(size_t i = 0; i < nSites; ++i)
      {
        if(fabs(oneint(i,i)) >= 1.0e-16) fdump << i << " " << i << " "
                                               << std::scientific << oneint(i,i) << std::endl;
        for(size_t j = 0; j < i; ++j)
        {
          if(fabs(oneint(i,j)) >= 1.0e-16) fdump << i << " " << j << " "
                                                 << std::scientific << oneint(i,j) << std::endl;
        }
      }
    }

    void DumpTwoInt(std::ofstream& fdump, const size_t& nSites)
    {
      fdump << nSites << std::endl;
      fdump.setf(std::ios::scientific, std::ios::floatfield); fdump.precision(20);
      for(size_t i = 0; i < nSites; ++i)
      {
        if(fabs(twoint(i,i,i,i)) >= 1.0e-16) fdump << i << " " << i << " " << i << " " << i << " "
                                                   << std::scientific << twoint(i,i,i,i) << std::endl;
        for(size_t j = 0; j < i; ++j)
        {
          if(fabs(twoint(i,i,i,j)) >= 1.0e-16) fdump << i << " " << i << " " << i << " " << j << " "
                                                     << std::scientific << twoint(i,i,i,j) << std::endl;
          if(fabs(twoint(i,j,j,j)) >= 1.0e-16) fdump << i << " " << j << " " << j << " " << j << " "
                                                     << std::scientific << twoint(i,j,j,j) << std::endl;
          if(fabs(twoint(i,i,j,j)) >= 1.0e-16) fdump << i << " " << i << " " << j << " " << j << " "
                                                     << std::scientific << twoint(i,i,j,j) << std::endl;
          if(fabs(twoint(i,j,i,j)) >= 1.0e-16) fdump << i << " " << j << " " << i << " " << j << " "
                                                     << std::scientific << twoint(i,j,i,j) << std::endl;
          for(size_t k = 0; k < j; ++k)
          {
            if(fabs(twoint(i,i,j,k)) >= 1.0e-16) fdump << i << " " << i << " " << j << " " << k << " "
                                                       << std::scientific << twoint(i,i,j,k) << std::endl;
            if(fabs(twoint(i,j,i,k)) >= 1.0e-16) fdump << i << " " << j << " " << i << " " << k << " "
                                                       << std::scientific << twoint(i,j,i,k) << std::endl;
            if(fabs(twoint(i,j,j,k)) >= 1.0e-16) fdump << i << " " << j << " " << j << " " << k << " "
                                                       << std::scientific << twoint(i,j,j,k) << std::endl;
            if(fabs(twoint(i,j,k,j)) >= 1.0e-16) fdump << i << " " << j << " " << k << " " << j << " "
                                                       << std::scientific << twoint(i,j,k,j) << std::endl;
            if(fabs(twoint(i,j,k,k)) >= 1.0e-16) fdump << i << " " << j << " " << k << " " << k << " "
                                                       << std::scientific << twoint(i,j,k,k) << std::endl;
            if(fabs(twoint(i,k,j,k)) >= 1.0e-16) fdump << i << " " << k << " " << j << " " << k << " "
                                                       << std::scientific << twoint(i,k,j,k) << std::endl;
            for(size_t l = 0; l < k; ++l)
            {
              if(fabs(twoint(i,j,k,l)) >= 1.0e-16) fdump << i << " " << j << " " << k << " " << l << " "
                                                         << std::scientific << twoint(i,j,k,l) << std::endl;
              if(fabs(twoint(i,j,l,k)) >= 1.0e-16) fdump << i << " " << j << " " << l << " " << k << " "
                                                         << std::scientific << twoint(i,j,l,k) << std::endl;
              if(fabs(twoint(i,k,j,l)) >= 1.0e-16) fdump << i << " " << k << " " << j << " " << l << " "
                                                         << std::scientific << twoint(i,k,j,l) << std::endl;
            }
          }
        }
      }
    }

    inline double OneInt(const size_t& iLabel, const size_t& jLabel)
    {
      // convert spin orbital label to spacial one
      size_t i = iLabel / 2;
      size_t j = jLabel / 2;
      return oneint(i, j);
    }

    inline double TwoInt(const size_t& iLabel, const size_t& jLabel, const size_t& kLabel, const size_t& lLabel)
    {
      // convert spin orbital label to spacial one
      size_t i = iLabel / 2;
      size_t j = jLabel / 2;
      size_t k = kLabel / 2;
      size_t l = lLabel / 2;
      return twoint(i, j, k, l);
    }
  };
};

#endif

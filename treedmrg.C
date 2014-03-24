#include <iostream>
#include <fstream>
#include <cstring>

#include <time_stamp.h>

#include "integrals.h"
#include "Input.h"
#include "TreeNetworks.h"
#include "DensityMatrix.h"
using namespace ttns;

int main(int argc, char* argv[])
{
   if(argc < 2)
   {
      std::cout << "Error: No Config File" << std::endl;
      return 1;
   }

   std::string config(argv[1]);
   Input ttnsInput(config);

   if(!ttnsInput.UseFciDump())
   {
      BTAS_DEBUG("main : Reading 1- and 2-integrals separately");
      std::ifstream fint1(ttnsInput.FileInt1());
      integrals::moeri::ReadOneInt(fint1, ttnsInput.OrbitalMap());
      std::ifstream fint2(ttnsInput.FileInt2());
      integrals::moeri::ReadTwoInt(fint2, ttnsInput.OrbitalMap());
   }
   else
   {
      BTAS_DEBUG("main : Reading integrals from FCIDUMP");
      std::ifstream fdump(ttnsInput.FciDump());
      integrals::moeri::ReadFromFciDump(fdump, ttnsInput.OrbitalMap());
   }

   BTAS_DEBUG("main : Setting density matrices");
   size_t nOrbs = integrals::moeri::oneint.extent(0);
   ttns::onepdm_alpha.resize(integrals::moeri::oneint.shape());
   ttns::onepdm_beta.resize(integrals::moeri::oneint.shape());

   BTAS_DEBUG("main : Doing TTNS optimization");

   time_stamp ts;

   const size_t z = ttnsInput.Rank();

   if(z == 2)
   {
      TreeNetworks<2> tree2(ttnsInput);
   }
   else if(z == 3)
   {
      TreeNetworks<3> tree3(ttnsInput);
   }
   else if(z == 4)
   {
      TreeNetworks<4> tree4(ttnsInput);
   }
   else
   {
      std::cout << "\tError: Tree with Z > 4 has not supported." << std::endl;
   }
   std::cout << "\t\tOptimization finished: " << std::fixed << std::setprecision(2) << ts.elapsed() << " sec. " << std::endl;

   std::cout << "\t\t\tSaving 1-Particle Reduced Density Matrix (Alpha)" << std::endl;
   std::ofstream ofsAlpha("onepdm_alpha.dat");
   ofsAlpha << nOrbs << std::endl;
   ofsAlpha.setf(std::ios::scientific, std::ios::floatfield); ofsAlpha.precision(20);
   for(size_t i = 0; i < nOrbs; ++i)
   {
      if(fabs(ttns::onepdm_alpha(i,i)) >= 1.0e-16) ofsAlpha << i << " " << i << " " << std::scientific << ttns::onepdm_alpha(i,i) << std::endl;
      for(size_t j = 0; j < i; ++j)
      {
         ofsAlpha << i << " " << j << " " << std::scientific << ttns::onepdm_alpha(i,j) << std::endl;
      }
   }

   std::cout << "\t\t\tSaving 1-Particle Reduced Density Matrix (Beta)" << std::endl;
   std::ofstream ofsBeta("onepdm_beta.dat");
   ofsBeta << nOrbs << std::endl;
   ofsBeta.setf(std::ios::scientific, std::ios::floatfield); ofsBeta.precision(20);
   for(size_t i = 0; i < nOrbs; ++i)
   {
      if(fabs(ttns::onepdm_beta(i,i)) >= 1.0e-16) ofsBeta << i << " " << i << " " << std::scientific << ttns::onepdm_beta(i,i) << std::endl;
      for(size_t j = 0; j < i; ++j)
      {
         ofsBeta << i << " " << j << " " << std::scientific << ttns::onepdm_beta(i,j) << std::endl;
      }
   }

   return 0;
}


#include <iostream>
#include <fstream>
#include <cstring>
#include "integrals.h"
#include "Input.h"
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

  int nSites = ttnsInput.Connection().size();

  std::string name1(ttnsInput.FileInt1());
  std::string name2(ttnsInput.FileInt2());
  {
    std::ifstream fint1(name1.c_str());
    integrals::moeri::ReadOneInt(fint1, ttnsInput.OrbitalMap());
    std::ifstream fint2(name2.c_str());
    integrals::moeri::ReadTwoInt(fint2, ttnsInput.OrbitalMap());
  }

  name1 += ".sorted";
  name2 += ".sorted";
  {
    std::ofstream fint1(name1.c_str());
    integrals::moeri::DumpOneInt(fint1, nSites);
    std::ofstream fint2(name2.c_str());
    integrals::moeri::DumpTwoInt(fint2, nSites);
  }

  return 0;
}

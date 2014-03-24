#include <iostream>
#include <fstream>
#include <random/uniform.h>

ranlib::Uniform< double > uniform;

void CreateIntegral(std::ostream& ost, const double& alpha,
                    const int& i, const int& j)
{
  double integral = alpha * uniform.random();
  if(fabs(integral) >= 1.0e-8)
    ost << i << " " << j << " " << integral << std::endl;
  return;
}

void CreateIntegral(std::ostream& ost, const double& alpha,
                    const int& i, const int& j, const int& k, const int& l)
{
  double integral = alpha * uniform.random();
  if(fabs(integral) >= 1.0e-8)
    ost << i << " " << j << " " << k << " " << l << " " << integral << std::endl;
  return;
}

int main(void)
{
  std::ofstream fint1("testTreeNetworks.int1");
  std::ofstream fint2("testTreeNetworks.int2");

  int nSites; std::cin >> nSites;

  fint1 << nSites << std::endl;
  fint2 << nSites << std::endl;
  for(int i = 0; i < nSites; ++i)
  {
    CreateIntegral(fint1,-0.50, i, i);
    CreateIntegral(fint2, 1.00, i, i, i, i);
    for(int j = i + 1; j < nSites; ++j)
    {
      double scale = 1.0/(j - i);
      CreateIntegral(fint1, 1.00*scale, i, j);
      CreateIntegral(fint2, 0.01*scale, i, i, i, j);
      CreateIntegral(fint2, 0.01*scale, i, j, j, j);
      CreateIntegral(fint2, 0.04*scale, i, j, i, j);
      CreateIntegral(fint2,-0.02*scale, i, i, j, j);
      for(int k = j + 1; k < nSites; ++k)
      {
        scale = 1.0/(j + k - 2 * i);
        CreateIntegral(fint2, 0.01*scale, i, i, j, k);
        CreateIntegral(fint2, 0.02*scale, i, j, i, k);
        scale = 1.0/(k - i);
        CreateIntegral(fint2, 0.01*scale, i, j, j, k);
        CreateIntegral(fint2, 0.02*scale, i, j, k, j);
        scale = 1.0/(2 * k - i - j);
        CreateIntegral(fint2, 0.01*scale, i, j, k, k);
        CreateIntegral(fint2, 0.02*scale, i, k, j, k);
        for(int l = k + 1; l < nSites; ++l)
        {
          scale = 1.0/(k + l - i - j);
          CreateIntegral(fint2, 0.10*scale, i, j, l, k);
          CreateIntegral(fint2, 0.10*scale, i, j, k, l);
          CreateIntegral(fint2, 0.10*scale, i, k, j, l);
        }
      }
    }
  }

  return 0;
}

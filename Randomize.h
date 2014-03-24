#ifndef __TTNS_RANDOMIZE_H
#define __TTNS_RANDOMIZE_H

#include <functional> // std::bind
#include <random> // std::mt19937, std::distribution

#include <btas/QSPARSE/QSDArray.h>

namespace ttns
{

// adding noise
template<size_t N, class Q>
void Randomize (const double& scale, btas::QSDArray<N, Q>& x)
{
   if(scale < 1.0e-12) return;

   std::mt19937 rgen;
   std::uniform_real_distribution<double> dist(-1.0, 1.0);

   btas::QSDArray<N, Q> noise(x.q(), x.qshape(), x.dshape(), std::bind(dist, rgen));

   btas::Axpy(scale, noise, x);
   btas::Normalize(x);
}

template<size_t N, class Q>
void Randomize (const double& scale, btas::QSDArray<N, Q>& x, const btas::TVector<btas::Dshapes, N>& dshape)
{
   if(scale < 1.0e-12) return;

   std::mt19937 rgen;
   std::uniform_real_distribution<double> dist(-1.0, 1.0);

   btas::QSDArray<N, Q> noise(x.q(), x.qshape(), dshape, std::bind(dist, rgen));

   btas::Axpy(scale, noise, x);
   btas::Normalize(x);
}

} // namespace ttns

#endif // __TTNS_RANDOMIZE_H

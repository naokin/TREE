#ifndef __TTNS_FERMIONIC_QUANTUM_H
#define __TTNS_FERMIONIC_QUANTUM_H 1

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

class Fermion
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_nptcl & m_nspin; }

public:
  const static Fermion zero() { return Fermion(0, 0); }

  Fermion() : m_nptcl(0), m_nspin(0) { }
  Fermion(int nptcl, int nspin) : m_nptcl(nptcl), m_nspin(nspin) { }

  bool operator== (const Fermion& other) const { return (m_nptcl == other.m_nptcl && m_nspin == other.m_nspin); }
  bool operator!= (const Fermion& other) const { return (m_nptcl != other.m_nptcl || m_nspin != other.m_nspin); }
  bool operator<  (const Fermion& other) const { return (m_nptcl == other.m_nptcl) ? (m_nspin < other.m_nspin) : (m_nptcl < other.m_nptcl); }
  bool operator>  (const Fermion& other) const { return (m_nptcl == other.m_nptcl) ? (m_nspin > other.m_nspin) : (m_nptcl > other.m_nptcl); }

  Fermion operator* (const Fermion& other) const { return Fermion(m_nptcl+other.m_nptcl, m_nspin+other.m_nspin); }
  Fermion operator+ (const Fermion& other) const { return Fermion(m_nptcl+other.m_nptcl, m_nspin+other.m_nspin); }

  Fermion operator+ () const { return Fermion(+m_nptcl, +m_nspin); }
  Fermion operator- () const { return Fermion(-m_nptcl, -m_nspin); }

  bool   parity () const { return m_nptcl & 1; }
  double clebsch() const { return 1.0; }

  const int& p() const { return m_nptcl; }
  const int& s() const { return m_nspin; }


private:
  int
    m_nptcl;
  int
    m_nspin;
};

inline std::ostream& operator<< (std::ostream& ost, const Fermion& q)
{
   return ost << "{" << std::setw(2) << q.p() << "," << std::setw(2) << q.s() << "}";
}

#endif // __TTNS_FERMIONIC_QUANTUM_H

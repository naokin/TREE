#ifndef __TTNS_BLOCK_H
#define __TTNS_BLOCK_H

#include <iostream>
#include <vector>
#include <map>

#include <btas/SPARSE/SDArray.h>

#include "ttns_assert.h"

namespace ttns
{

typedef btas::SDArray<2> DmrgOperator;

typedef btas::TArray<DmrgOperator, 1> OpArray1d;

typedef btas::TArray<DmrgOperator, 2> OpArray2d;

/// Block operators
/// For some complications in OpComponents*.h, use SDArray, i.e. not quantum number-based array
/// This is not fully object-oriented, and thus, not really safe,
/// while, it can reduce overheads from quantum number contraction.
class Block
{

private:

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive& ar, const unsigned int version)
   {
      ar & i_index_ & o_index_;
      ar & op_Ham_;
      ar & op_Cre_ & op_Cre_Comp_;
      ar & op_Cre_Cre_ & op_Cre_Cre_Comp_;
      ar & op_Cre_Des_ & op_Cre_Des_Comp_;
   }

protected:

   std::map<size_t, size_t> i_index_; ///< index belongs to this block

   std::map<size_t, size_t> o_index_; ///< complement of i_index_

   DmrgOperator op_Ham_;

   OpArray1d op_Cre_;
   OpArray1d op_Cre_Comp_;

   OpArray2d op_Cre_Cre_;
   OpArray2d op_Cre_Cre_Comp_;

   OpArray2d op_Cre_Des_;
   OpArray2d op_Cre_Des_Comp_;

   template<class Archive>
   void load (Archive& ar, const unsigned int version)
   {
      ar >> i_index_ >> o_index_;
      ar >> op_Ham_;
      ar >> op_Cre_ >> op_Cre_Comp_;
      ar >> op_Cre_Cre_ >> op_Cre_Cre_Comp_;
      ar >> op_Cre_Des_ >> op_Cre_Des_Comp_;
   }

   template<class Archive>
   void save (Archive& ar, const unsigned int version)
   {
      ar << i_index_ << o_index_;
      ar << op_Ham_;
      ar << op_Cre_ << op_Cre_Comp_;
      ar << op_Cre_Cre_ << op_Cre_Cre_Comp_;
      ar << op_Cre_Des_ << op_Cre_Des_Comp_;

      op_Ham_.clear();
      op_Cre_.clear();
      op_Cre_Comp_.clear();
      op_Cre_Cre_.clear();
      op_Cre_Cre_Comp_.clear();
      op_Cre_Des_.clear();
      op_Cre_Des_Comp_.clear();
   }

private:

   static size_t n_orbs_;

   template<typename... Args>
   void __resize_by_args(std::vector<size_t>& i_orbs, const size_t& i, const Args&... ids)
   {
      i_orbs.push_back(i);
      __resize_by_args(i_orbs, ids...);
   }

   inline void __resize_by_args(std::vector<size_t>& i_orbs, const size_t& i)
   {
      i_orbs.push_back(i);
      this->resize(i_orbs);
   }

public:

   static size_t& orbitals () { return n_orbs_; }

   /// Default constructor
   Block ();

   /// Destructor
  ~Block ();

   /// Copy constructor
   Block (const Block& x)
   :  i_index_ (x.i_index_),
      o_index_ (x.o_index_),
      op_Ham_ (x.op_Ham_),
      op_Cre_ (x.op_Cre_),
      op_Cre_Comp_ (x.op_Cre_Comp_),
      op_Cre_Cre_ (x.op_Cre_Cre_),
      op_Cre_Cre_Comp_ (x.op_Cre_Cre_Comp_),
      op_Cre_Des_ (x.op_Cre_Des_),
      op_Cre_Des_Comp_ (x.op_Cre_Des_Comp_)
   { }

   /// Copy assignment
   Block& operator= (const Block& x);

   /// Construct by arguments
   template<typename... Args>
   Block (const size_t& i, const Args&... ids)
   {
      std::vector<size_t> i_orbs;
      __resize_by_args(i_orbs, i, ids...);
   }

   /// resize by orbitals inside
   Block (const std::vector<size_t>& i_orbs);

   /// resize by arguments
   template<typename... Args>
   void resize (const size_t& i, const Args&... ids)
   {
      std::vector<size_t> i_orbs;
      __resize_by_args(i_orbs, i, ids...);
   }

   /// resize by orbitals inside
   void resize (const std::vector<size_t>& i_orbs);

   /// return size
   inline size_t size () const { return i_index_.size(); }

   /// return index inside
   std::vector<size_t> inside () const;

   /// return index outside
   std::vector<size_t> outside () const;

   /// return Ham
   inline DmrgOperator& Ham ()
   { return op_Ham_; }

   /// return Ham
   inline const DmrgOperator& Ham () const
   { return op_Ham_; }

   /// return Cre
   inline DmrgOperator& Cre (size_t i)
   { return op_Cre_(i_index_.find(i)->second); }

   /// return Cre
   inline const DmrgOperator& Cre (size_t i) const
   { return op_Cre_(i_index_.find(i)->second); }

   /// return Cre Comp
   inline DmrgOperator& CreComp (size_t i)
   { return op_Cre_Comp_(o_index_.find(i)->second); }

   /// return Cre Comp
   inline const DmrgOperator& CreComp (size_t i) const
   { return op_Cre_Comp_(o_index_.find(i)->second); }

   /// return Cre Cre
   inline DmrgOperator& CreCre (size_t i, size_t j)
   { return op_Cre_Cre_(i_index_.find(i)->second, i_index_.find(j)->second); }

   /// return Cre Cre
   inline const DmrgOperator& CreCre (size_t i, size_t j) const
   { return op_Cre_Cre_(i_index_.find(i)->second, i_index_.find(j)->second); }

   /// return Cre Cre Comp
   inline DmrgOperator& CreCreComp (size_t i, size_t j)
   { return op_Cre_Cre_Comp_(o_index_.find(i)->second, o_index_.find(j)->second); }

   /// return Cre Cre Comp
   inline const DmrgOperator& CreCreComp (size_t i, size_t j) const
   { return op_Cre_Cre_Comp_(o_index_.find(i)->second, o_index_.find(j)->second); }

   /// return Cre Des
   inline DmrgOperator& CreDes (size_t i, size_t j)
   { return op_Cre_Des_(i_index_.find(i)->second, i_index_.find(j)->second); }

   /// return Cre Des
   inline const DmrgOperator& CreDes (size_t i, size_t j) const
   { return op_Cre_Des_(i_index_.find(i)->second, i_index_.find(j)->second); }

   /// return Cre Des Comp
   inline DmrgOperator& CreDesComp (size_t i, size_t j)
   { return op_Cre_Des_Comp_(o_index_.find(i)->second, o_index_.find(j)->second); }

   /// return Cre Des Comp
   inline const DmrgOperator& CreDesComp (size_t i, size_t j) const
   { return op_Cre_Des_Comp_(o_index_.find(i)->second, o_index_.find(j)->second); }

   void clear ();
};

} // namespace ttns

/// printing function
std::ostream& operator<< (std::ostream& ost, const ttns::Block& block);

#endif // __TTNS_BLOCK_H

#include "Block.h"

#include <algorithm>

/// Default constructor
ttns::Block::Block () { }

/// Destructor
ttns::Block::~Block () { }

/// Copy assignment
ttns::Block& ttns::Block::operator= (const ttns::Block& x)
{
   i_index_ = x.i_index_;
   o_index_ = x.o_index_;

   op_Ham_ = x.op_Ham_;

   op_Cre_ = x.op_Cre_;
   op_Cre_Comp_ = x.op_Cre_Comp_;

   op_Cre_Cre_ = x.op_Cre_Cre_;
   op_Cre_Cre_Comp_ = x.op_Cre_Cre_Comp_;
   op_Cre_Des_ = x.op_Cre_Des_;
   op_Cre_Des_Comp_ = x.op_Cre_Des_Comp_;

   return *this;
}

/// resize by orbitals inside
ttns::Block::Block (const std::vector<size_t>& i_orbs)
{
   this->resize(i_orbs);
}

/// resize by orbitals inside
void ttns::Block::resize (const std::vector<size_t>& i_orbs)
{
   std::vector<size_t> i_sort_orbs(i_orbs);
   std::sort(i_sort_orbs.begin(), i_sort_orbs.end());

   i_index_.clear();

   size_t ni = 0;

   for(size_t i = 0; i < i_sort_orbs.size(); ++i)
   {
      i_index_.insert(std::make_pair(i_sort_orbs[i], ni++));
   }

   o_index_.clear();

   size_t no = 0;

   for(size_t i = 0; i < n_orbs_; ++i)
   {
      if(i_index_.find(i) != i_index_.end()) continue;

      o_index_.insert(std::make_pair(i, no++));
   }

   op_Ham_.clear();

   op_Cre_.clear();
   op_Cre_Comp_.clear();

   op_Cre_Cre_.clear();
   op_Cre_Cre_Comp_.clear();

   op_Cre_Des_.clear();
   op_Cre_Des_Comp_.clear();

   op_Cre_.resize(ni);
   op_Cre_Comp_.resize(no);

   op_Cre_Cre_.resize(ni, ni);
   op_Cre_Cre_Comp_.resize(no, no);

   op_Cre_Des_.resize(ni, ni);
   op_Cre_Des_Comp_.resize(no, no);
}

/// return index inside
std::vector<size_t> ttns::Block::inside () const
{
   std::vector<size_t> i_orbs; i_orbs.reserve(i_index_.size());

   for(auto it = i_index_.begin(); it != i_index_.end(); ++it) i_orbs.push_back(it->first);

   return i_orbs;
}

/// return index outside
std::vector<size_t> ttns::Block::outside () const
{
   std::vector<size_t> o_orbs; o_orbs.reserve(o_index_.size());

   for(auto it = o_index_.begin(); it != o_index_.end(); ++it) o_orbs.push_back(it->first);

   return o_orbs;
}

/// Deallocate storages
void ttns::Block::clear ()
{
   i_index_.clear();
   o_index_.clear();

   op_Ham_.clear();
   op_Cre_.clear();
   op_Cre_Comp_.clear();
   op_Cre_Cre_.clear();
   op_Cre_Cre_Comp_.clear();
   op_Cre_Des_.clear();
   op_Cre_Des_Comp_.clear();
}

/// printing function
std::ostream& operator<< (std::ostream& ost, const ttns::Block& x)
{
   auto loop_index = x.inside();
   ost << "orbitals (inside): "; for(auto it = loop_index.begin(); it != loop_index.end(); ++it) ost << *it << ","; ost << std::endl;

   auto comp_index = x.outside();
   ost << "orbitals (outside): "; for(auto it = comp_index.begin(); it != comp_index.end(); ++it) ost << *it << ","; ost << std::endl;

   return ost;
}

size_t ttns::Block::n_orbs_ = 0ul;


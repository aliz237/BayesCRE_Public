// Relation.cpp
#include "headers.h"

//--------------------------------------------------------------------------------------------------------------
Relation :: Relation(const string& relf, const string& dimf)
{
  ifstream fin(relf);
  if (!fin) ferror("init_rels can't open file: " + relf + "\n");

  Index i;
  while (fin >> i)
    rels.push_back(i);
  
  fin.close();

  fin.open(dimf);
  if (!fin) ferror("init_rels can't open file: " + dimf + "\n");

  for (i=0; i<4; ++i)
    fin >> unique_key_dim[i];
  
  fin.close();
}
//--------------------------------------------------------------------------------------------------------------

Hashtab Relation :: group_by (enum rel_idx idx) const
{
  Index sz = rels.size();
  
  Hashtab h(unique_dim(idx));
    
  for (Index i=0; i<sz; i += R_N)
      h.insert(rels[i+idx], i);
  
  return h;
}
//--------------------------------------------------------------------------------------------------------------
Index Relation :: unique_dim (enum rel_idx idx) const
{
  if (idx == R_MESHID) return unique_key_dim[2];
  if (idx == R_APPID ) return unique_key_dim[3];
  return unique_key_dim[idx];
}
  
//--------------------------------------------------------------------------------------------------------------
Index Relation :: current_child_trg(Index idx) const
{
  return rels[idx + R_TRGIND];
}
//--------------------------------------------------------------------------------------------------------------
Index Relation :: current_child_app(Index idx) const
{
  return rels[idx + R_APPIND];
}
//--------------------------------------------------------------------------------------------------------------
set<Index> Relation :: markov_blanket(const vector<Index>& v, enum rel_idx idx) const
{
  set<Index> res;
  
  for (const Index& i : v)
    res.insert(rels[i+idx]);
  
  return res;
}
//--------------------------------------------------------------------------------------------------------------
vector<Index> Relation :: get_idx(const vector<Index>& v, enum rel_idx idx) const
{
  vector<Index> s;
  
  for (const Index& i : v)
    s.push_back(rels[i+idx]);

  return s;
}
//--------------------------------------------------------------------------------------------------------------

#include "headers.h"

Hashtab::Hashtab(Index sz):
    tab {vector<vector<Node>>(next_prime(sz))} {}
//----------------------------------------------------------------------------------------------------------
Index Hashtab::hash_map(Index i) {return i % tab.size();}
//----------------------------------------------------------------------------------------------------------
void Hashtab :: insert (Index key, Index val)
{
  Index i = hash_map (key);

  auto pos = std::find_if(tab[i].begin(), tab[i].end(), [key](const Node& it)
			  { return it.key == key; });
			  
  if (pos != tab[i].end())
    pos->idx.push_back(val);
  
  else
    {
      Node n;
      n.idx.push_back(val);
      n.key = key;
      tab[i].push_back(n);
    }
}
//----------------------------------------------------------------------------------------------------------

vector<Index> Hashtab :: find (Index key)
{
  Index i = hash_map (key);

  auto pos = find_if (tab[i].begin(), tab[i].end(), [key](const Node& n)
		      { return n.key == key; });
  
  if (pos == tab[i].end())
    herror("can't find" + std::to_string(key) + "in hash table.\n");

  return pos->idx;
}

//----------------------------------------------------------------------------------------------------------

Index Hashtab :: next_prime(Index i)
{
  if (i < 2) return 2;
  if (i % 2 == 0) ++i;

  while (!is_prime(i)) i += 2;

  return i;
}

//----------------------------------------------------------------------------------------------------------

bool Hashtab :: is_prime (Index i)
{
  // this is called with i > 2
  Index end = ceil(sqrt(i));
  
  for (Index j=3; j<=end; j += 2)
    if ( i % j == 0)
      return false;
  
  return true;
}    


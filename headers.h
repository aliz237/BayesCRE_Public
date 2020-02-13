#ifndef GIBBS_SAMPLING_H
#define GIBBS_SAMPLING_H

#include <vector>
#include <string>
#include <set>

#include <fstream>
#include <iostream>
//#include <stdint.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/special_functions/binomial.hpp>

/*---------------------------------------------------------------------------------------------------------*/

using std :: vector;
using std :: array;
using std :: string;
using std :: set;
using std :: ifstream;

typedef boost :: numeric :: ublas :: matrix<string,boost :: numeric :: ublas :: row_major, vector<string>> Matrix2s;
typedef boost :: numeric :: ublas :: matrix<double> Matrix2d;
typedef boost :: numeric :: ublas :: zero_matrix<double> Zmatrix2d;
typedef boost :: numeric :: ublas :: matrix_column<Matrix2s> Col1s;
typedef int Index; // only of sizeof(int) == 4

/*---------------------------------------------------------------------------------------------------------*/
struct Node
{
  vector<Index> idx;
  Index key;
};

class Hashtab
{
  vector<vector<Node>> tab;
  Index next_prime(Index);
  bool is_prime(Index);
  Index hash_map(Index i);
  //Index binary_search(Index key, Index i);

 public:
  
  Hashtab(Index sz);
  void insert(Index key, Index val);
  vector<Index> find (Index key);
  // don't need delete, rehash, etc ...(fixed sized dataset)
};
/*---------------------------------------------------------------------------------------------------------*/

enum rel_idx : char
  {
      R_SRCUID,
      R_TRGUID,
      R_TYPE,
      R_MESHID,
      R_APPID,
      R_SRCIND,
      R_TRGIND,
      R_MESHIND,
      R_APPIND,
      R_TYPE_SRCIND,
      R_N
  };

class Relation
{

  vector<Index> rels;
  
  // this holds the number of unique elements in
  // SRC, TRG, MESH, APP columns
  Index unique_key_dim[4]; 
  Index unique_dim (enum rel_idx) const;
  
  public:
  
  Relation (const string& relf, const string& dimf);
  Hashtab group_by (enum rel_idx) const;
  Index current_child_trg (Index i) const;
  Index current_child_app(Index  i) const;
  set<Index> markov_blanket(const vector<Index>& v, enum rel_idx idx) const;
  vector<Index> get_idx(const vector<Index>& v, enum rel_idx idx) const;  
  friend vector<Index> applicable_edges (const vector<Index>& v);
};

/*---------------------------------------------------------------------------------------------------------*/

enum ents_idx : char
{
  E_UID,
  E_NAME,
  E_ID,
  E_TYPE
};

struct Entry
{
  vector<Index> pc_idx;
  vector<Index> m_idx;
  vector<Index> a_idx;
  vector<Index> uid;
  // vector<Index> non_zero_mesh;
};

/*---------------------------------------------------------------------------------------------------------*/

struct Prob{

  double pc; // probability of the edge being wrong for nonzero edges
  double pa; // probability of the edge being wrong for zero edges
  double pm; // prior probability of the ture state of the genes being -1
  double pz; // prior probability of the ture state of the genes being 0
  double pp; // prior probability of the ture state of the genes being 1
  double hz; //prior probability of the hypothesis being off
  double alpha; // false positive rate
  double beta; // false negative rate
  double w; // weight of the sigmoid
  double cp;
  
  Prob(const string& name);
};

/*---------------------------------------------------------------------------------------------------------*/

struct Hash_error
{
string name;
Hash_error(string s) : name{s} {}
Hash_error(const char * s) : name{s} {}
};
inline void herror (string p){ throw Hash_error(p); }
//std :: ostream& operator << (std :: ostream& o, const rel& r);
/*---------------------------------------------------------------------------------------------------------*/

struct File_error
{
  string name;
  File_error(const char* n) : name(n) {}
  File_error(string s) : name(s) {}
};

inline void ferror (const char* p){ throw File_error(p); }
inline void ferror (string p){ throw File_error(p); }

std :: pair<int, int> dim(const string& name);
void read_matrix(const string& name, Matrix2s& m);
void init_matrix(const string& name, Matrix2s& m);

/*---------------------------------------------------------------------------------------------------------*/
/*
 * These are the ... */
void init();
std::pair<double, double> comb(int np, int nm, int k, int l, double pc);
void comp_prob_H(int nm, int np, double c);
void comp_pZ(int z);
double comp_ch_p(int nm, int np, double pac, int val);
std::pair<Index, Index> comp_nm_np(const vector<Index>& Val, const vector<Index>& pa_lab);
std::pair<Index, Index> comp_nm_np(const vector<Index>& n);
#endif

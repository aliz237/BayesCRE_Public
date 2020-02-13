#include "headers.h"

Prob :: Prob(const string& name)
{
  ifstream fin(name);
  if (!fin)
    {
      string msg = "Cannot open file:" + name + "\n";
      ferror(msg);
    }

  fin >> pc >> pa >> hz >> alpha >> beta >> w >> cp >> pm >> pp >> pz;
  fin.close();
}

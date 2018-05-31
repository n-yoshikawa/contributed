#include <iostream>
#include <fstream>
#include <algorithm>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>
#include <openbabel/elements.h>

using namespace std;
using namespace OpenBabel;

bool asc(pair<int, pair<OBSmartsPattern*, vector<vector3> > >& left,
         pair<int, pair<OBSmartsPattern*, vector<vector3> > >& right) {
    return left.first  < right.first;
}

int main() {
  ifstream ifs("fragments.txt");
  if(!ifs) {
    cerr << "Falied to open fragments.txt" << endl;
    exit(EXIT_FAILURE);
  }
  vector<pair<int, pair<OBSmartsPattern*, vector<vector3> > > > known_fragments;

  char buffer[BUFF_SIZE];
  vector<string> vs;
  OBSmartsPattern *sp = NULL;
  vector<vector3> coords;
  while (ifs.getline(buffer, BUFF_SIZE)) {
    if (buffer[0] == '#')
      continue;

    tokenize(vs, buffer);

    if (vs.size() == 1) { // SMARTS pattern
      if (sp != NULL)
        known_fragments.push_back(pair<int, pair<OBSmartsPattern*, vector<vector3> > > (sp->NumAtoms(), pair<OBSmartsPattern*, vector<vector3> > (sp, coords)));

      coords.clear();
      sp = new OBSmartsPattern;
      if (!sp->Init(vs[0])) {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);
      }
    } else if (vs.size() == 4) { // XYZ coordinates
      vector3 coord(atof(vs[1].c_str()), atof(vs[2].c_str()), atof(vs[3].c_str()));
      coords.push_back(coord);
    } else {
      cerr << "Unexpected input" << endl;
      exit(EXIT_FAILURE);
    }
  }

  known_fragments.push_back(pair<int, pair<OBSmartsPattern*, vector<vector3> > > (sp->NumAtoms(), pair<OBSmartsPattern*, vector<vector3> > (sp, coords)));
  sort(known_fragments.rbegin(), known_fragments.rend(), asc);
  for(size_t i=0; i<known_fragments.size(); ++i) {
    cout << known_fragments[i].second.first->GetSMARTS() << endl;
    for(size_t j=0; j<known_fragments[i].second.second.size(); ++j) {
      cout << known_fragments[i].second.second[j].GetX() << " " 
           << known_fragments[i].second.second[j].GetY() << " " 
           << known_fragments[i].second.second[j].GetZ() << endl;
    }
  }
}


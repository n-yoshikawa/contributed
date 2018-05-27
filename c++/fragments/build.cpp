#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <queue>
#include <map>
#include <vector>
#include <string>

using namespace std;
using namespace OpenBabel;

vector<pair<string, OBMol> > LoadFragments() {
    OBConversion conv;
    conv.SetInFormat("sdf");
    conv.SetOutFormat("can");
    OBMol mol;
    vector<pair<string, OBMol> > fragments;

    bool notatend = conv.ReadFile(&mol, "fragments.sdf");
    while(notatend) {
        string smiles = conv.WriteString(&mol, true);
        fragments.push_back(make_pair(smiles, mol));
        mol.Clear();
        notatend = conv.Read(&mol);
    }
    return fragments;
}

int main(int argc,char *argv[]) {
    std::ios::sync_with_stdio(false);


    if (argc < 1) {
        cout << "Usage: " << argv[0] << " <file>" << endl;
        return -1;
    }

    OBConversion conv;
    OBFormat *inFormat, *canFormat;
    canFormat = conv.FindFormat("can"); // Canonical SMILES format
    conv.SetOutFormat(canFormat);
    //conv.SetOptions("O", conv.OUTOPTIONS);

    OBMol mol;
    ifstream ifs;

    vector<pair<string, OBMol> > known_fragments = LoadFragments();

    for (int i = 1; i < argc; i++) {
        cerr << "Reading file " << argv[i] << endl;

        inFormat = conv.FormatFromExt(argv[i]);
        if(inFormat==NULL || !conv.SetInFormat(inFormat)) {
            cerr << " Cannot read file format for " << argv[i] << endl;
            continue; // try next file
        }

        ifs.open(argv[i]);

        if (!ifs) {
            cerr << "Cannot read input file: " << argv[i] << endl;
            continue;
        }


        while(ifs.peek() != EOF && ifs.good()) {
            conv.Read(&mol, &ifs);

            int size = mol.NumAtoms();
            OBBitVec atomsToCopy(size);
            for (unsigned int i = 1; i <= size; ++i) {
                OBAtom* atom = mol.GetAtom(i);
                atomsToCopy.SetBitOn(atom->GetIdx());
            }

            size = mol.NumBonds();
            OBBitVec bondsToExclude(size);
            for (unsigned int i = 0; i < size; ++i) {
                OBBond* bond = mol.GetBond(i);
                if (bond->IsRotor()) {
                    bondsToExclude.SetBitOn(bond->GetIdx());
                }
            }


            OBMol mol_copy;
            mol.CopySubstructure(mol_copy, &atomsToCopy, &bondsToExclude);
            vector<OBMol> fragments = mol_copy.Separate(); // Copies each disconnected fragment as a separate
            OBBitVec in_frag(mol.NumAtoms()); // Is the atoms in fragment?
            cout << "Original: " << conv.WriteString(&mol, true) << endl;
            for (unsigned int i = 0; i < fragments.size(); ++i) {
                string fragment_smiles = conv.WriteString(&fragments[i], true);

                // need better search method
                vector<pair<string, OBMol> >::iterator f;
                for(f = known_fragments.begin(); f != known_fragments.end(); f++) {
                    if((*f).first == fragment_smiles) {
                        cout << "Found: " << fragment_smiles << endl;

                        // This search may be eliminated?
                        OBSmartsPattern sp;
                        sp.Init(fragment_smiles);
                        sp.Match(mol);

                        vector<vector<int> > mlist; // match list for fragments
                        mlist = sp.GetUMapList();

                        vector<vector<int> >::iterator j;
                        for (j = mlist.begin();j != mlist.end();++j) { // for all matches
                            vector<int>::iterator k;
                            for(k = j->begin(); k != j->end(); ++k) {
                                in_frag.SetBitOn(*k);
                                cout << *k << " ";
                            }
                            cout << endl;
                        }
                    }
                }
            }
            for(size_t i = 1; i < mol.NumAtoms()+1; ++i) {
                cout << (in_frag.BitIsSet(i) ? 1 : 0) << endl;
            }
        }
    }
}

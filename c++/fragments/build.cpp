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

vector<string> LoadFragments() {
    ifstream ifs("fragments.txt");
    string line;
    vector<string> fragments;

    while(ifs && getline(ifs, line)) {
        fragments.push_back(line);
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

    vector<string> known_fragments = LoadFragments();

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

            cout << "1" << endl;

            OBMol mol_copy;
            mol.CopySubstructure(mol_copy, &atomsToCopy, &bondsToExclude);
            vector<OBMol> fragments = mol_copy.Separate(); // Copies each disconnected fragment as a separate
            cout << "2" << endl;
            for (unsigned int i = 0; i < fragments.size(); ++i) {
                cout << "3" << endl;
                string fragment_smiles = conv.WriteString(&fragments[i], true);
                cout << "4" << endl;
                cout << "Search: " << fragment_smiles << endl;
                if (find(known_fragments.begin(), known_fragments.end(), 
                         fragment_smiles) == known_fragments.end()) {
                    cout << "5" << endl;
                    cout << "Not found: " << fragment_smiles << endl;
                } else {
                    cout << "6" << endl;
                    cout << "Found: " << fragment_smiles << endl;
                }
            }
        }
    }
}



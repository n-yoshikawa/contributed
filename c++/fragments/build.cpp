#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>
#include <openbabel/elements.h>

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

vector3 GetCorrectedBondVector(OBAtom *atom1, OBAtom *atom2, int bondOrder = 1)
  {
    double bondLength = 0.0;

    // We create an estimate of the bond length based on the two atoms
    // Scaling is performed by the bond order corrections below
    //  .. so we will use the straight covalent radii
    bondLength += OBElements::GetCovalentRad(atom1->GetAtomicNum());
    bondLength += OBElements::GetCovalentRad(atom2->GetAtomicNum());

    if (bondLength < 1.0)
      bondLength = 1.0;

    // These are based on OBBond::GetEquibLength
    // Numbers come from averaged values of Pyykko and Atsumi
    if (bondOrder == -1) // aromatic
      bondLength *= 0.9475;   // 0.9475 = average of 1.0 and 0.8950
    else if (bondOrder == 2)
      bondLength *= 0.8950;   // 0.8950
    else if (bondOrder == 3)
      bondLength *= 0.8578;   // 0.8578

    return OBBuilder::GetNewBondVector(atom1, bondLength);
  }

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


    while (ifs.peek() != EOF && ifs.good()) {
      conv.Read(&mol, &ifs);

      // Copy to a working molecule
      OBMol workMol = mol;

      // Treat a 2D structure like a 0D one
      // what is this mean?
      if (workMol.GetDimension() == 2) {
        workMol.SetDimension(0);
      }

      // Delete all bonds
      while (workMol.NumBonds()) {
        workMol.DeleteBond(workMol.GetBond(0));
      }

      // Deleting the bonds unsets HybridizationPerceived. To prevent our
      // perceived values being reperceived (incorrectly), we must set
      // this flag again.
      workMol.SetHybridizationPerceived();

      // Get fragments using Copy Substructure
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
      // Copies each disconnected fragment as a separate
      vector<OBMol> fragments = mol_copy.Separate(); 

      OBBitVec vfrag(mol.NumAtoms()); // Is the atoms in fragment?
      cout << "Original: " << conv.WriteString(&mol, true) << endl;

      // Need to implement skip or something here

      // Loop through the fragments and assign the coordinates
      for (unsigned int i = 0; i < fragments.size(); ++i) {
        string fragment_smiles = conv.WriteString(&fragments[i], true);

        // need better search method
        vector<pair<string, OBMol> >::iterator f;
        for(f = known_fragments.begin(); f != known_fragments.end(); f++) {
          if((*f).first == fragment_smiles) {
            cout << "Found: " << fragment_smiles << endl;

            // This search may be omitted?
            OBSmartsPattern sp;
            sp.Init(fragment_smiles);
            sp.Match(mol);

            vector<vector<int> > mlist; // match list for fragments
            mlist = sp.GetUMapList();

            vector<vector<int> >::iterator j;
            for (j = mlist.begin();j != mlist.end();++j) { // for all matches
              bool alreadydone = false;
              vector<int>::iterator k, k2;
              for (k = j->begin(); k != j->end(); ++k) {
                if (vfrag.BitIsSet(*k)) {
                  alreadydone = true;
                  break;
                }
              }
              if (alreadydone) continue;

              int counter;
              for (k = j->begin(), counter=1; k != j->end(); ++k, ++counter) {
                OBAtom *atom = workMol.GetAtom(*k);
                OBAtom *atom2 = (*f).second.GetAtom(counter);
                cout << "atom: " << atom->GetAtomicNum() << " " << atom2->GetAtomicNum() << endl;
                atom->SetVector(atom2->GetVector()); // is this ok?
                vfrag.SetBitOn(*k); // flag for done
              }
              int index2;
              for (k = j->begin(); k != j->end(); ++k) {
                OBAtom *atom1 = mol.GetAtom(*k);

                for (k2 = j->begin(); k2 != j->end(); ++k2) {
                  index2 = *k2;
                  OBAtom *atom2 = mol.GetAtom(index2);
                  OBBond *bond = atom1->GetBond(atom2);

                  if (bond != NULL) {
                    workMol.AddBond(*bond);
                  }
                }
              }
            }
          }
        }
      }
      OBBitVec vdone;
      FOR_DFS_OF_MOL (a, mol) {
        if (vdone.BitIsSet(a->GetIdx())) // continue if the atom is already added
          continue;
        
        // find an atom connected to the current atom that is already added
        OBAtom *prev = NULL;
        FOR_NBORS_OF_ATOM (nbr, &*a) {
          if (vdone.BitIsSet(nbr->GetIdx()))
            prev = &*nbr;
        }

        if (vfrag.BitIsSet(a->GetIdx())) { // Is this atom part of a fragment?
          if (prev != NULL) { // if we have a previous atom, translate/rotate the fragment and connect it
            OBBuilder::Connect(workMol, prev->GetIdx(), a->GetIdx(), mol.GetBond(prev, &*a)->GetBondOrder());
            // set the correct bond order
            int bondOrder = mol.GetBond(prev->GetIdx(), a->GetIdx())->GetBondOrder();
            workMol.GetBond(prev->GetIdx(), a->GetIdx())->SetBondOrder(bondOrder);
          }
          OBBitVec fragment = OBBuilder::GetFragment(workMol.GetAtom(a->GetIdx()));
          vdone |= fragment; // mark this fragment as done

          continue;
        }

        // If this atom is not a part of a fragment
        // get the position for the new atom, this is done with GetNewBondVector
        vector3 molvec, moldir;
        if (prev != NULL) {
          int bondType = a->GetBond(prev)->GetBO();
          if (a->GetBond(prev)->IsAromatic())
            bondType = -1;

          molvec = GetCorrectedBondVector(workMol.GetAtom(prev->GetIdx()),
              workMol.GetAtom(a->GetIdx()),
              bondType);
          moldir = molvec - workMol.GetAtom(prev->GetIdx())->GetVector();
        } else {
          // We don't want to plant all base atoms at exactly the same spot.
          // (or in exactly the same direction)
          // So we'll add a slight tweak -- fixes problem reported by Kasper Thofte
          vector3 randomOffset;
          randomOffset.randomUnitVector();
          molvec = VX + 0.1 * randomOffset;
          moldir = VX + 0.01 * randomOffset;
        }

        vdone.SetBitOn(a->GetIdx());
        // place the atom
        OBAtom *atom = workMol.GetAtom(a->GetIdx());
        atom->SetVector(molvec);

        // add bond between previous part and added atom
        if (prev != NULL) {
          OBBond *bond = a->GetBond(prev); // from mol
          workMol.AddBond(*bond);
        }

      }

      // Make sure we keep the bond indexes the same
      // so we'll delete the bonds again and copy them
      // Fixes PR#3448379 (and likely other topology issues)
      while (workMol.NumBonds())
        workMol.DeleteBond(workMol.GetBond(0));

      int beginIdx, endIdx;
      FOR_BONDS_OF_MOL(b, mol) {
        beginIdx = b->GetBeginAtomIdx();
        endIdx = b->GetEndAtomIdx();
        workMol.AddBond(beginIdx, endIdx, b->GetBO(), b->GetFlags());
      }

      // correct the chirality
      bool success = OBBuilder::CorrectStereoBonds(workMol);
      // we only succeed if we corrected all stereochemistry
      bool stereoWarnings = false;
      success = success && OBBuilder::CorrectStereoAtoms(workMol, stereoWarnings);

      /*
      // if the stereo failed, we should use distance geometry instead
      OBDistanceGeometry dg;
      dg.Setup(workMol);
      dg.GetGeometry(workMol); // ensured to have correct stereo
      */

      mol = workMol;
      mol.SetChiralityPerceived();
      mol.SetDimension(3);

      bool isNanExist = false;
      FOR_ATOMS_OF_MOL(a, mol) {
        vector3 v = a->GetVector();
        if(IsNan(v.x()) || IsNan(v.y()) || IsNan(v.z())) {
          isNanExist = true;
          break;
        }
      }
      if(isNanExist)
        obErrorLog.ThrowError(__FUNCTION__, "There exists NaN in calculated coordinates.", obWarning);
    OBConversion conv2;
    conv2.SetOutFormat("sdf");
    cout << conv2.WriteString(&mol, true) << endl;
    }

  }
}

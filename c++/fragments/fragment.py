import sys

import pybel
ob = pybel.ob

atomsToCopy = ob.OBBitVec()
bondsToExclude = ob.OBBitVec()

def CopyRingSystem(orig, copy):
    atomsToCopy.Clear()
    for atom in ob.OBMolAtomIter(orig):
        if atom.IsInRing():
            atomsToCopy.SetBitOn(atom.GetIdx())
    bondsToExclude.Clear()
    for bond in ob.OBMolBondIter(orig):
        if not bond.IsInRing():
            bondsToExclude.SetBitOn(bond.GetIdx())
    return orig.CopySubstructure(copy, atomsToCopy, bondsToExclude)

if __name__ == "__main__":
    fragments = {}
    copied = ob.OBMol()
    for mol in pybel.readfile("sdf", sys.argv[1]):
        copied.Clear()
        ok = CopyRingSystem(mol.OBMol, copied)
        assert ok
        if copied.NumAtoms() > 0:
            smi = pybel.Molecule(copied).write("smi").rstrip()
            for f in smi.split("."):
                if f in fragments:
                    fragments[f] += 1
                else:
                    fragments[f] = 1
    for f in sorted(fragments.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f)
                


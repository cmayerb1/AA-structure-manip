# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 12:21:13 2021.

@author: umbci
"""
"""
Calculates conformational energy of each molecule, and returns a .sdf file with the input molecules sorted by increasing energy of their lowest energy conformer.  Uses MMFF94 forcefield implementation in RDKit.

Usage:
...>python conf_energy_screen.py <input filename as .sdf> -out <output filename as .sdf>
"""
"""
Module imports
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
import argparse

"""
Function definitions
"""

def minConfEnergy(mol):
    """Find lowest energy (kcal/mol) conformation of the input molecule.
    
    From an input molecule, generates 50 conformers, minimizes and calculates their conformational energies using MMFF94 forcefield, then returns the lowest energy.
    

    Parameters
    ----------
    mol : RDKit Mol object

    Returns
    -------
    Conformational energy (kcal/mol) (float)

    """
    mol_h = Chem.AddHs(mol, addCoords=True) #add hydrogens, create coordinates
    param = rdDistGeom.ETKDGv2()
    cids = rdDistGeom.EmbedMultipleConfs(mol_h, 50, param)
    mp = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant='MMFF94')
    AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=0, mmffVariant='MMFF94') #optimize all conformations
    res = [] # list of conformer energies
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append(e)
    return min(res)

"""
Main Program
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates conformer energies from a user-defined input .sdf file.  Returns those same molecules sorted by lowest conformer energy.')
    parser.add_argument("in_file", help="Input file name (.sdf format) (required)")
    parser.add_argument("-out", help="Name of output file (.sdf format) (required)")
    args = parser.parse_args()
    # load in data file
    data = Chem.SDMolSupplier(args.in_file, removeHs=False)
    smiles = [Chem.MolToSmiles(mol) for mol in data]
    energies = []
    for x in range(len(data)):
        print(f"Processing molecule {x+1} of {len(data)}", end='\r', flush=True)
        energies.append(minConfEnergy(data[x]))
    energies_sorted = [e for e, s in sorted(zip(energies, smiles))]
    mols_sorted = [Chem.MolFromSmiles(s) for e, s in sorted(zip(energies, smiles))]
    # create SD file writer, add energy property to each structure
    writer = Chem.SDWriter(args.out)
    for x in range(len(mols_sorted)):
        mols_sorted[x].SetDoubleProp("Energy", energies_sorted[x])
        writer.write(mols_sorted[x])
    writer.close()
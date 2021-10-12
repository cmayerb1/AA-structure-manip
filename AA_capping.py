# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:27:28 2021.

@author: umbci
"""
"""
Short script for adding amide bond caps to alpha-amino acids.  Concept described in Ilardo et al. (2015) (https://doi.org/10.1038/srep09414) using MOLGEN; current approach using RDKit used in Mayer-Bacon and Freeland (2021) (https://doi.org/10.1016/j.jtbi.2021.110661).
"""
"""
Module imports
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions
import argparse
import re
"""
Function definitions
"""
def aaCapping(mol):
    """Caps the N- and C-termini of an alpha-amino acid (add an acetyl group and N-methyl group, respectively).
    
    Intended to mimic an amino acid in the polypeptide context, i.e. with amide bonds instead of free amine and carboxyl groups.
    
    Parameters
    ----------
    mol : RDKit Mol object
    

    Returns
    -------
    RDKit Mol object with capped C- and N-termini if input is an alpha-amino acid; returns input molecule if not an alpha-amino acid

    """
    # set patterns corresponding to amino acid skeletons
    patt_1a = Chem.MolFromSmarts('[NX3:1][CX4H]([*:3])[CX3:2](=[OX1])[O]') # monosubstituted alpha-AA
    patt_gly = Chem.MolFromSmarts('[NX3:1][CX4H2][CX3:2](=[OX1])[O]') # glycine
    patt_sec_gly = Chem.MolFromSmarts('[*:3][NX3:1][CX4H2][CX3:2](=[OX1])[O]') # N-substituted glycines
    patt_sec_amine = Chem.MolFromSmarts('[*:3][NX3:1][CX4H]([*:4])[CX3:2](=[OX1])[O]') # N-substituted, monosubstituted
    # recalculate implicit valence and ring info before transformation
    mol.UpdatePropertyCache()
    if mol.HasSubstructMatch(patt_1a):
        cap = rdChemReactions.ReactionFromSmarts('[NX3:1][CX4H]([*:3])[CX3:2](=[OX1])[O]>>[NX3:1](C(=O)C)[CX4H]([*:3])[CX3:2](=[OX1])NC')
    elif mol.HasSubstructMatch(patt_gly):
        cap = rdChemReactions.ReactionFromSmarts('[NX3:1][CX4H2][CX3:2](=[OX1])[O]>>[NX3:1](C(=O)C)[CX4H2][CX3:2](=[OX1])NC')
    elif mol.HasSubstructMatch(patt_sec_gly):
        cap = rdChemReactions.ReactionFromSmarts('[*:3][NX3:1][CX4H2][CX3:2](=[OX1])[O]>>[*:3][NX3:1](C(=O)C)[CX4H2][CX3:2](=[OX1])NC')
    elif mol.HasSubstructMatch(patt_sec_amine):
        cap = rdChemReactions.ReactionFromSmarts('[*:3][NX3:1][CX4H]([*:4])[CX3:2](=[OX1])[O]>>[*:3][NX3:1](C(=O)C)[CX4H]([*:4])[CX3:2](=[OX1])NC')
    else:
        cap = None
    if cap != None:
        return cap.RunReactant(mol,0)[0][0]
    else:
        return mol

def bracketRemove(smi):
    """Remove brackets from string.
    
    Removes brackets that surround single characters in a string.  Written to remove brackets from around single atoms in some SMILES strings where they are not useful to the program.

    Parameters
    ----------
    smi : string

    Returns
    -------
    Input string without brackets surrounding single letters

    """
    m = re.sub(r"\[(\w)]", r'\1', smi)
    return m
"""
Main program
"""
if __name__ == "__main__":
    # create argument parser, parse command line arguments
    parser = argparse.ArgumentParser(description='Caps alpha-amino acids from a user-defined input .sdf file.')
    parser.add_argument("in_file", help="Input file name (.sdf format) (required)")
    args = parser.parse_args()
    # load in data file
    data = Chem.SDMolSupplier(args.in_file, removeHs=False)
    # create file writer
    w = Chem.SDWriter("args.in_file"+"_capped.sdf")
    # loop through each input structure, cap and write to file
    for mol in data:
        mol2 = Chem.MolFromSmiles(bracketRemove(Chem.MolToSmiles(mol)))
        mol_cap = aaCapping(mol2)
        w.write(mol_cap)
    w.close()
    

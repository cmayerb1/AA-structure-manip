# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 14:02:49 2021.

@author: umbci
"""
"""
Short script to replace a pseudoatom (currently P) with another, constant substructure.  Currently configured to substitute a trivalent P atom with an alpha-amino backbone plus beta-carbon (e.g. an alanine moiety).  Future work could involve implementing user-specified pseudoatoms or replacement substructures.
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
def pseudoReplace(mol):
    """Replace pseudoatom with alpha-amino acid backbone.
    
    Replaces trivalent phosphorous atom with an alpha-amino acid backbone and beta-carbon.

    Parameters
    ----------
    mol : RDKit Mol object
        Mol object with trivalent P pseudo atom

    Returns
    -------
    RDKit Mol object with P pseudo atom replaced.  Returns input Mol if P atom is not present.

    """
    # SMARTS patterns for different numbers of P heavy-atom neighbors
    three_neighbor = Chem.MolFromSmarts("[PD3,pD3]")
    two_neighbor = Chem.MolFromSmarts("[PD2,pD2]")
    one_neighbor = Chem.MolFromSmarts("[PD1,pD1]")
    three_h = Chem.MolFromSmarts("[PX3&H3]")
    # first match # of P non-H neighbors, then bonding patterns among those neighbors; return input molecule if no matches
    if mol.HasSubstructMatch(three_neighbor): # P with 0 hydrogens
        rxn = rdChemReactions.ReactionFromSmarts("[PD3]([*:2])([*:3])[*:1]>>NC(C(=O)O)C([*:2])([*:3])[*:1]")
        return rxn.RunReactant(mol,0)[0][0]
    elif mol.HasSubstructMatch(two_neighbor): # P with 2 neighbors
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[PD2,pD2]=[*]")):
            rxn = rdChemReactions.ReactionFromSmarts("[PD2,pD2]([*:2])=[*:1]>>NC(C(=O)O)C([*:2])=[*:1]")
            return rxn.RunReactant(mol,0)[0][0]
        else:
            rxn = rdChemReactions.ReactionFromSmarts("[PD2,pD2]([*:2])[*:1]>>NC(C(=O)O)C([*:2])[*:1]")
            return rxn.RunReactant(mol,0)[0][0]
    elif mol.HasSubstructMatch(one_neighbor): # P with 1 neighbor
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[PD1]#[*]")):
            rxn = rdChemReactions.ReactionFromSmarts("[PD1]#[*:1]>>NC(C(=O)O)C#[*:1]")
            return rxn.RunReactant(mol,0)[0][0]
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[PD1]=[*]")):
            rxn = rdChemReactions.ReactionFromSmarts("[PD1]=[*:1]>>NC(C(=O)O)C=[*:1]")
            return rxn.RunReactant(mol,0)[0][0]
        else:
            rxn = rdChemReactions.ReactionFromSmarts("[PD1][*:1]>>NC(C(=O)O)C[*:1]")
            return rxn.RunReactant(mol,0)[0][0]
    elif mol.HasSubstructMatch(three_h): # P with 3 hydrogens, replace with Ala
        return Chem.MolFromSmiles("NC(C(=O)O)C")
    else: # if P pseudo atom isn't present in input molecule
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
    w = Chem.SDWriter(args.in_file+"_pseudo-replace.sdf")
    for mol in data:
        mol2 = Chem.MolFromSmiles(bracketRemove(Chem.MolToSmiles(mol)))
        new_mol = pseudoReplace(mol2)
        w.write(new_mol)
    w.close()
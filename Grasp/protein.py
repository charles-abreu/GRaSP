'''
    Protein module
    Store protein features and data
    Author: Charles Abreu Santana
    Last update: 2020-04-29
'''
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB import is_aa
from Bio.PDB import NeighborSearch

class Protein:
    def __init__(self, file_name):
        self.file_name = file_name
        # get pdb_id in file name
        self.pdb_id = os.path.basename(self.file_name).split('.')[0]
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure(self.pdb_id, self.file_name)
        self.init_matrix()
        self.neighborhood = Neighborhood().run(self, 6)

    ##########################################################
    # Gets the fasta sequence from structure file
    # In: .pdb file
    # Out: respective .fasta file for each chain in .pdb file
    ##########################################################
    def get_fasta(self, out_dir):
        self.fasta_files = []
        # For each chain crate a .fasta file
        for record in SeqIO.parse(self.file_name, "pdb-atom"):
            record.id = self.pdb_id[:4] + record.annotations['chain']
            fasta_output = out_dir + os.sep + record.id  + '.fasta'
            # lista com respectivos arquivos fasta
            self.fasta_files.append(fasta_output)
            SeqIO.write(record, fasta_output , "fasta")

    def set_templates(self, templates):
        self.templates = templates

    def init_matrix(self):
        # Properties
        #properties = ['DON','HPB','ACP','NEG','ARM','POS']
        #list all residues
        residue_list = list(self.structure[0].get_residues())
        self.matrix = pd.DataFrame()#columns=properties
        self.matrix.index.name = 'res_name'

        for residue in residue_list:
            if is_aa(residue):
                # Residue id: pdbID_chain_resName_resNumber
                residue_name = get_residue_index(residue)
                self.matrix.loc[residue_name,:] = 0.0

class Neighborhood:
    def run(self, protein, distance):
        self.adj_list = {}
        # list all atons to compute neighbor residues
        self.atomlist = list(protein.structure.get_atoms())
        # list all residues
        self.residue_list = list(protein.structure.get_residues())
        # object to compute neighbor residues
        self.n = NeighborSearch(self.atomlist)

        for residue in self.residue_list:
            if is_aa(residue):
                # pdbID_residueChain_residueName_residueNumber
                residues_name = get_residue_index(residue)
                # list of neighbor residues
                neighbor_list = []
                for atom in residue.get_atoms():
                    # compute neighbors and add in the list
                    neighbor_list += self.n.search(atom.get_coord(), distance, level='R')

                # save neighbor residues names and delete copies
                neighbor_set = set()
                for neighbor in neighbor_list:
                    if is_aa(neighbor):
                        # pdbID_residueChain_residueName_residueNumber
                        neighbor_name = get_residue_index(neighbor)
                        # Evita que adicione o proprio residuo
                        if neighbor_name != residues_name:
                            neighbor_set.add(neighbor_name)

                self.adj_list[residues_name] = neighbor_set
        return self.adj_list

#######################################################
# Toolbox Functions
#######################################################
# return residue index on format pdb_chain_name_number
def get_residue_index(residue):
    pdb_id = residue.get_full_id()[0]
    chain = residue.get_parent().get_id()
    name = residue.get_resname()
    num = str(residue.get_id()[1])

    return pdb_id + '_' + chain + '_' + name + '_' + num

if __name__ == '__main__':
    p = Protein('test/test.pdb')
    col = ['acc_rel', 'acc_side', 'acc_main', 'acc_apolar', 'acc_polar']
    for c in col:
        p.matrix[c] = np.nan
    print(p.matrix)

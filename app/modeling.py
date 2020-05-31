'''
    Modeling module
    Perform graph representation and compute Interactions
    and properties
    Author: Charles Abreu Santana
    Last update: 2020-04-29
'''
import numpy as np
from Bio.PDB import is_aa
import pandas as pd
import networkx as nx
import os

class Properties:
    def __init__(self, data_dir):
        self.types = read_atom_types(data_dir)

    def compute_properties(self, protein):
        # Insert new columns in protein matrix
        columns = ['DON','HPB','ACP','NEG','ARM','POS']
        for c in columns: protein.matrix[c] = 0
        # for each residue in pdb
        for residue_name in protein.matrix.index:
            # Split residue name
            pdb_id, res_chain, res_name, res_number = residue_name.split('_')
            res_number = int(res_number)
            # Verify if the residue exixts on protein data structure
            if res_number in protein.structure[0][res_chain]:
                residue = protein.structure[0][res_chain][res_number]
                # Counting atom properties to current residue
                for atom in residue:
                    atom_name = residue.get_resname() + '_' + atom.get_name()
                    # Verify if the atom has properties
                    if atom_name in self.types:
                        # Donor
                        if 'DON' in self.types[atom_name]:
                            protein.matrix.loc[residue_name,'DON'] += 1
                        # Hydrophobic
                        if 'HPB' in self.types[atom_name]:
                            protein.matrix.loc[residue_name,'HPB'] += 1
                        # Acceptor
                        if 'ACP' in self.types[atom_name]:
                            protein.matrix.loc[residue_name,'ACP'] += 1
                        # Negative
                        if 'NEG' in self.types[atom_name]:
                            protein.matrix.loc[residue_name,'NEG'] += 1
                        # Aroatic
                        if 'ARM' in self.types[atom_name]:
                            protein.matrix.loc[residue_name,'ARM'] += 1
                        # Positive
                        if 'POS' in self.types[atom_name]:
                            protein.matrix.loc[residue_name,'POS'] += 1

class Interactions:
    def __init__(self, data_dir):
        self.atom_types = read_atom_types(data_dir)

    def compute_interactions(self, protein):
        # Insert new columns in protein matrix
        columns = ['aromatic_stacking','dissulfide_bridge','hydrogen_bond',
                        'hydrophobic', 'repulsive', 'salt_bridge']
        for c in columns: protein.matrix[c] = 0

        # Residue Graph
        protein.residue_graph = nx.Graph()

        # Data frame auxiliar para verificar vertices visitados
        residue_info_df = pd.DataFrame(columns = ['visited'], index = protein.matrix.index)
        residue_info_df['visited'] = False

        for residue_name in protein.matrix.index:
            pdb_id, res_chain, res_name, res_number = residue_name.split('_')
            res_number = int(res_number)

            # Verify if the residue exixts on protein data structure
            if res_number in protein.structure[0][res_chain]:
                residue = protein.structure[0][res_chain][res_number]
                # Computing interactions
                for adj_res_name in protein.neighborhood[residue_name]:
                    # Verify is residue exists in data frame
                    if adj_res_name in protein.matrix.index:
                        # Verify if the residues was visitated
                        if not residue_info_df.loc[adj_res_name, 'visited']:
                            # split residue name
                            pdb_id, adj_chain, adj_name , adj_number = adj_res_name.split('_')
                            adj_number = int(adj_number)
                            # Verify if the neighbor residue exixts on protein data structure
                            if adj_number in protein.structure[0][adj_chain]:
                                adj_residue = protein.structure[0][adj_chain][adj_number]
                                # Access residue atoms
                                for atom_1 in residue:
                                    atom1_name = res_name + '_' + atom_1.get_name().strip()
                                    if atom1_name in self.atom_types:
                                        # Access neighbor residue atoms
                                        for atom_2 in adj_residue:
                                            atom2_name = adj_name + '_' + atom_2.get_name().strip()
                                            if atom2_name in self.atom_types:
                                                # Distance between atons
                                                distance = atom_1 - atom_2
                                                #testing interaction types
                                                if ('ARM' in self.atom_types[atom1_name]) and ('ARM' in self.atom_types[atom2_name]) and (distance >= 1.5) and (distance <= 3.5):
                                                    protein.matrix.loc[residue_name,'aromatic_stacking'] += 1
                                                    protein.matrix.loc[adj_res_name,'aromatic_stacking'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('ACP' in self.atom_types[atom1_name]) and ('DON' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.0):
                                                    protein.matrix.loc[residue_name,'hydrogen_bond'] += 1
                                                    protein.matrix.loc[adj_res_name,'hydrogen_bond'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('DON' in self.atom_types[atom1_name]) and ('ACP' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.0):
                                                    protein.matrix.loc[residue_name,'hydrogen_bond'] += 1
                                                    protein.matrix.loc[adj_res_name,'hydrogen_bond'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('HPB' in self.atom_types[atom1_name]) and ('HPB' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.8):
                                                    protein.matrix.loc[residue_name,'hydrophobic'] += 1
                                                    protein.matrix.loc[adj_res_name,'hydrophobic'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('POS' in self.atom_types[atom1_name]) and ('POS' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    protein.matrix.loc[residue_name,'repulsive'] += 1
                                                    protein.matrix.loc[adj_res_name,'repulsive'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('NEG' in self.atom_types[atom1_name]) and ('NEG' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    protein.matrix.loc[residue_name,'repulsive'] += 1
                                                    protein.matrix.loc[adj_res_name,'repulsive'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('NEG' in self.atom_types[atom1_name]) and ('POS' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    protein.matrix.loc[residue_name,'salt_bridge'] += 1
                                                    protein.matrix.loc[adj_res_name,'salt_bridge'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('POS' in self.atom_types[atom1_name]) and ('NEG' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    protein.matrix.loc[residue_name,'salt_bridge'] += 1
                                                    protein.matrix.loc[adj_res_name,'salt_bridge'] += 1
                                                    protein.residue_graph.add_edge(residue_name,adj_res_name)
            residue_info_df.loc[residue_name,'visited'] = True

class Layer:
    #########################################################
    # Generate the data to respective layers
    #########################################################
    def compute_layer(self, protein, layers):
        # list of data frames with layers information
        layer_df = []
        for l in range(layers):
            layer_df.append(self.sum_properties(protein, l+1 , True))
        for df in layer_df:
            protein.matrix = pd.merge(protein.matrix, df,  left_index=True, right_index=True)

        protein.matrix = protein.matrix.loc[:,:'N2_POS']

    ##################################################
    # Return residues from specific layer
    ##################################################
    def compute_neighbor_layer(self,residue_graph, residue_node, layer):
        # base
        if layer == 1:
            return set(residue_graph.neighbors(residue_node))
        # Recursion
        else:
            neighbors = set()
            # compute the same for each neighbor
            for n in residue_graph.neighbors(residue_node):
                neighbors = neighbors | self.compute_neighbor_layer(residue_graph, n, layer-1)
            return list(neighbors)

    #########################################################
    # Sum the properties from a specific residue's layer
    #########################################################
    def sum_properties(self, protein, layer, weighted):
        layer_name = 'N' + str(layer) + '_'
        # Create new df to store sum information
        new_residue_df = pd.DataFrame(0.0, index = protein.matrix.index,
            columns = [layer_name + c for c in protein.matrix.columns])
        new_residue_df.index.name = 'res_name'
        # Dataframe to control who is neighbor
        is_neighbor = pd.DataFrame(0,index = protein.matrix.index, columns=['is_neighbor'])

        for residue in protein.matrix.index:
            if not residue in protein.residue_graph:
                continue
            # Compute neighbors from especific layer
            neighbors = self.compute_neighbor_layer(protein.residue_graph, residue, layer)
            is_neighbor.loc[neighbors,'is_neighbor'] = 1
            # Sum neighbors properties
            if weighted:
                sum_neighbors = (protein.matrix[is_neighbor['is_neighbor'] == 1]).sum()/len(neighbors)
            else:
                sum_neighbors = (protein.matrix[is_neighbor['is_neighbor'] == 1]).sum()

            for col in sum_neighbors.index:
                new_name = layer_name + col # new name to colum sum
                new_residue_df.loc[residue, new_name] = sum_neighbors[col]

            is_neighbor['is_neighbor'] = 0
        return new_residue_df

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

###################################################################
# Ler arquivos com as propriedades de cada atomo para cada residuo
# Retorna um dicionario com o nome do residuo e suas propriedades
###################################################################
def read_atom_types(data_dir):
    types = {}
    with open(os.path.join(data_dir,"atom_types.csv")) as in_file:
        for line in in_file:
            record = line.strip().split(',')
            atomName = record[0] + '_' + record[1]
            # Atomos sem propriedades
            if len(record) < 3:
                continue
            else:
                types[atomName] = record[2:]
    return types

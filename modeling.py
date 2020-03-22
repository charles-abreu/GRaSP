from Bio.PDB import *
import numpy as np
import networkx as nx
import pandas as pd
import os
import sys
import random
from glob import glob
import getopt

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

########################################################################
# Compute neighborhood
########################################################################
# Get the residue neighborhood
def get_neighborhood(protein_struct, distance):

    adj_list = {}

    # list all atons to compute neighbor residues
    atomlist = list(protein_struct.get_atoms())
    # list all residues
    residue_list = list(protein_struct.get_residues())
    # object to compute neighbor residues
    n = NeighborSearch(atomlist)

    for residue in residue_list:
        if is_aa(residue):
            # pdbID_residueChain_residueName_residueNumber
            residues_name = get_residue_index(residue)
            # list of neighbor residues
            neighbor_list = []
            for atom in residue.get_atoms():
                # compute neighbors and add in the list
                neighbor_list += n.search(atom.get_coord(), distance, level='R')

            # save neighbor residues names and delete copies
            neighbor_set = set()
            for neighbor in neighbor_list:
                if is_aa(neighbor):
                    # pdbID_residueChain_residueName_residueNumber
                    neighbor_name = get_residue_index(neighbor)
                    # Evita que adicione o proprio residuo
                    if neighbor_name != residues_name:
                        neighbor_set.add(neighbor_name)

            adj_list[residues_name] = neighbor_set

    return adj_list

################################
# Iniciando Matriz de resíduos
################################
def init_residue_matrix(protein_struct):
    # Properties
    properties = ['DON','HPB','ACP','NEG','ARM','POS']

    #list all residues
    residue_list = list(protein_struct[0].get_residues())
    residue_matrix = pd.DataFrame(columns=properties)

    for residue in residue_list:
        if is_aa(residue):
            # Residue id: pdbID_chain_resName_resNumber
            residue_name = get_residue_index(residue)
            residue_matrix.loc[residue_name,:] = 0.0

    return residue_matrix

########################################################################
# Ler arquivos com as propriedades de cada atomo para cada residuo
# Retorna um dicionario com o nome do residuo e suas propriedades
########################################################################
def read_atom_types():
    types = {}

    in_file = open('atom_types.csv', 'r')

    for line in in_file:
        record = line.strip().split(',')
        atomName = record[0] + '_' + record[1]
        # Atomos sem propriedades
        if len(record) < 3:
            continue
        else:
            types[atomName] = record[2:]

    return types

########################################################
# Compute all properties for each residue from pdb file
########################################################
def compute_properties(protein_struct):
    # dta frame with residues properties
    residue_df = init_residue_matrix(protein_struct)
    # atom types
    atom_types = read_atom_types()

    # for each residue in pdb
    for residue_name in residue_df.index:
        # Split residue name
        pdb_id, res_chain, res_name, res_number = residue_name.split('_')
        res_number = int(res_number)
        # Verify if the residue exixts on protein data structure
        if res_number in protein_struct[0][res_chain]:
            residue = protein_struct[0][res_chain][res_number]

            # Counting atom properties to current residue
            for atom in residue:
                atom_name = residue.get_resname() + '_' + atom.get_name()
                # Verify if the atom has properties
                if atom_name in atom_types:
                    # Donor
                    if 'DON' in atom_types[atom_name]:
                        residue_df.loc[residue_name,'DON'] += 1
                    # Hydrophobic
                    if 'HPB' in atom_types[atom_name]:
                        residue_df.loc[residue_name,'HPB'] += 1
                    # Acceptor
                    if 'ACP' in atom_types[atom_name]:
                        residue_df.loc[residue_name,'ACP'] += 1
                    # Negative
                    if 'NEG' in atom_types[atom_name]:
                        residue_df.loc[residue_name,'NEG'] += 1
                    # Aroatic
                    if 'ARM' in atom_types[atom_name]:
                        residue_df.loc[residue_name,'ARM'] += 1
                    # Positive
                    if 'POS' in atom_types[atom_name]:
                        residue_df.loc[residue_name,'POS'] += 1

    residue_df.index.name = 'res_name'
    return residue_df

################################################
# Compute ASA imformation
################################################
def compute_rsa(pdb_id, rsa_file):
    # lendo informacao da area acessivel
    # arquivo no formato .rsa
    in_file = open(rsa_file)

    # Data fame com a informacao da acesssibilidade dos residuos
    residue_df = pd.DataFrame(columns=['acc_rel', 'acc_side', 'acc_main', 'acc_apolar', 'acc_polar'])
    residue_df.index.name = 'res_name'

    for line in in_file:
        # Linhas com informacao sobre a area acessivel
        if line.startswith('RES') or line.startswith('HEM'):
            # All-atoms
            area = float(line[22:28].strip())

            # Ignora registros com area negativa
            if area < 0:
                continue
            # Total-Side
            area_side = float(line[36:41].strip())
            # Main-Chain
            area_main_chain = float(line[49:54].strip())
            # Non-polar
            area_apolar = float(line[62:67].strip())
            # All polar
            area_polar = float(line[75:80].strip())

            res_name = line[4:7].strip()
            res_chain = line[8]
            res_number = line[9:13].strip()
            residue_id = pdb_id + '_' + res_chain + '_' + res_name + '_' + res_number
            # Add informacao ao dataframe
            residue_df = residue_df.append(pd.Series({'acc_rel':area,
                                                     'acc_side':area_side,
                                                     'acc_main':area_main_chain,
                                                     'acc_apolar':area_apolar,
                                                     'acc_polar':area_polar}).rename(residue_id))
    in_file.close()

    return residue_df

################################################
# Compute Exposure imformation
################################################
def compute_exposure(protein_struct):
    # Data frame with exposure information
    residue_df = pd.DataFrame(columns=['hseu', 'hsed'])
    residue_df.index.name = 'res_name'
    # Compute residues exposure
    exp_cb = HSExposure.HSExposureCB(protein_struct[0])

    for key in exp_cb.keys():
        res_chain = key[0]
        res_number = key[1][1]
        try:
            res_name = protein_struct[0][res_chain][res_number].get_resname()
        except KeyError:
            continue

        residue_id = protein_struct.get_id() + '_' + res_chain + '_' + res_name + '_' + str(res_number)

        residue_df.loc[residue_id,'hseu'] = float(exp_cb[key][0])
        residue_df.loc[residue_id,'hsed'] = float(exp_cb[key][1])

    return residue_df

#############################
# Computing Interactions between residues
#############################
def compute_interactions(protein_struct, residue_df):

    # step 1 : compute neighborhood
    residue_adj = get_neighborhood(protein_struct, 6.0)

    # Interactions
    interactions = ['aromatic_stacking','dissulfide_bridge','hydrogen_bond',
                    'hydrophobic', 'repulsive', 'salt_bridge']

    # Add new columns to data frame with interactions
    for inter in interactions:
        residue_df[inter] = 0.0

    # atom types
    atom_types = read_atom_types()

    # Residue Graph
    residue_graph = nx.Graph()

    # Data frame auxiliar para verificar vertices visitados
    residue_info_df = pd.DataFrame(columns = ['visited'], index = residue_df.index)
    residue_info_df['visited'] = False

    for residue_name in residue_df.index:
        pdb_id, res_chain, res_name, res_number = residue_name.split('_')
        res_number = int(res_number)

        # Verify if the residue exixts on protein data structure
        if res_number in protein_struct[0][res_chain]:
            residue = protein_struct[0][res_chain][res_number]
            # Computing interactions
            for adj_res_name in residue_adj[residue_name]:
                # Verify is residue exists in data frame
                if adj_res_name in residue_df.index:
                    # Verify if the residues was visitated
                    if not residue_info_df.loc[adj_res_name, 'visited']:
                        # split residue name
                        pdb_id, adj_chain, adj_name , adj_number = adj_res_name.split('_')
                        adj_number = int(adj_number)
                        # Verify if the neighbor residue exixts on protein data structure
                        if adj_number in protein_struct[0][adj_chain]:
                            adj_residue = protein_struct[0][adj_chain][adj_number]
                            # Access residue atoms
                            for atom_1 in residue:
                                atom1_name = res_name + '_' + atom_1.get_name().strip()
                                if atom1_name in atom_types:
                                    # Access neighbor residue atoms
                                    for atom_2 in adj_residue:
                                        atom2_name = adj_name + '_' + atom_2.get_name().strip()
                                        if atom2_name in atom_types:
                                            # Distance between atons
                                            distance = atom_1 - atom_2
                                            #testing interaction types
                                            if ('ARM' in atom_types[atom1_name]) and ('ARM' in atom_types[atom2_name]) and (distance >= 1.5) and (distance <= 3.5):
                                                residue_df.loc[residue_name,'aromatic_stacking'] += 1
                                                residue_df.loc[adj_res_name,'aromatic_stacking'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('ACP' in atom_types[atom1_name]) and ('DON' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.0):
                                                residue_df.loc[residue_name,'hydrogen_bond'] += 1
                                                residue_df.loc[adj_res_name,'hydrogen_bond'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('DON' in atom_types[atom1_name]) and ('ACP' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.0):
                                                residue_df.loc[residue_name,'hydrogen_bond'] += 1
                                                residue_df.loc[adj_res_name,'hydrogen_bond'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('HPB' in atom_types[atom1_name]) and ('HPB' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.8):
                                                residue_df.loc[residue_name,'hydrophobic'] += 1
                                                residue_df.loc[adj_res_name,'hydrophobic'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('POS' in atom_types[atom1_name]) and ('POS' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                residue_df.loc[residue_name,'repulsive'] += 1
                                                residue_df.loc[adj_res_name,'repulsive'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('NEG' in atom_types[atom1_name]) and ('NEG' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                residue_df.loc[residue_name,'repulsive'] += 1
                                                residue_df.loc[adj_res_name,'repulsive'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('NEG' in atom_types[atom1_name]) and ('POS' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                residue_df.loc[residue_name,'salt_bridge'] += 1
                                                residue_df.loc[adj_res_name,'salt_bridge'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            if ('POS' in atom_types[atom1_name]) and ('NEG' in atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                residue_df.loc[residue_name,'salt_bridge'] += 1
                                                residue_df.loc[adj_res_name,'salt_bridge'] += 1
                                                residue_graph.add_edge(residue_name,adj_res_name)
                                            #if ('SSB' in atom_types[atom1_name]) and ('SSB' in atom_types[atom2_name]) and (distance >= 1.4) and (distance <= 2.3):
                                            #    residue_df.loc[residue_name,'dissulfide_bridge'] += 1
                                            #    residue_df.loc[adj_res_name,'dissulfide_bridge'] += 1

        residue_info_df.loc[residue_name,'visited'] = True
    # Graph with real interactions
    return residue_graph

##################################################
# Return residues from specific layer
##################################################
def compute_neighbor_layer(residue_graph, residue_node, layer):
    # base
    if layer == 1:
        return set(residue_graph.neighbors(residue_node))
    # Recursion
    else:
        neighbors = set()
        # compute the same for each neighbor
        for n in residue_graph.neighbors(residue_node):
            neighbors = neighbors | compute_neighbor_layer(residue_graph, n, layer-1)

        return list(neighbors)

#########################################################
# Sum the properties from a specific residue's layer
#########################################################
def sum_properties(residue_df, residue_graph, layer, weighted):
    columns = residue_df.columns
    layer_name = 'N' + str(layer) + '_'

    # Create new df to store sum information
    new_residue_df = pd.DataFrame(0.0, index = residue_df.index, columns = [layer_name + c for c in residue_df.columns])
    new_residue_df.index.name = 'res_name'

    # Dataframe to control who is neighbor
    is_neighbor = pd.DataFrame(0,index = residue_df.index, columns=['is_neighbor'])

    for residue in residue_df.index:
        if not residue in residue_graph:
            continue
        # Compute neighbors from especific layer
        neighbors = compute_neighbor_layer(residue_graph, residue, layer)

        for n in neighbors:
            is_neighbor.loc[n,'is_neighbor'] = 1

        # Sum neighbors properties
        if weighted:
            sum_neighbors = (residue_df[is_neighbor['is_neighbor'] == 1]).sum()/len(neighbors)
        else:
            sum_neighbors = (residue_df[is_neighbor['is_neighbor'] == 1]).sum()

        for col in sum_neighbors.index:
            # new name to colum sum
            new_name = layer_name + col
            new_residue_df.loc[residue, new_name] = sum_neighbors[col]

        is_neighbor['is_neighbor'] = 0

    return new_residue_df

#########################################################
# Generate the data to respective layers
#########################################################
def compute_layer(residue_df, residue_graph, layers):
    # list of data frames with layers information
    layer_df = []

    for l in range(layers):
        layer_df.append(sum_properties(residue_df,residue_graph, l+1 , True))

    complete_df = layer_df[0]

    for df in layer_df[1:]:
        complete_df = pd.merge(complete_df, df,  left_index=True, right_index=True)

    return complete_df

#########################################
# 				MAIN 					#
#########################################
def compute_matrix(pdb_file, rsa_file, out_dir):
    pdb_id = os.path.basename(pdb_file).split('.')[0]
    # Creating protein struct from Biopython
    parser = PDBParser(QUIET=True)
    protein_struct = parser.get_structure(pdb_id, pdb_file)
    # STEP 1: Compute properties
    print('Computing residue properties to ' + pdb_id)
    residue_properties = compute_properties(protein_struct)
    # STEP 2:  Compute asa
    print('Computing asa information to ' + pdb_id)
    residue_asa = compute_rsa(protein_struct.get_id(), rsa_file)
    # STEP 3: Compute exposure
    #residue_exo = compute_exposure(protein_struct)
    # STEP 4: Compute Interactions
    print('Computing Interactions to ' + pdb_id)
    residue_graph = compute_interactions(protein_struct, residue_properties)
    # STEP 5 : merge matrix
    residue_data = pd.merge(residue_asa, residue_properties, left_index=True, right_index=True)
    #residue_data = pd.merge(residue_exo, residue_data, on = 'res_name')
    # STEP 6: Compute layers
    print('Computing layers properties to ' + pdb_id)
    residue_data = pd.merge(residue_data, compute_layer(residue_data, residue_graph, 2), left_index=True, right_index=True)
    residue_data = residue_data.loc[:,:'N2_POS']
    # STEP 8: Save .csv file
    print('Saving matrix to ' + pdb_id)
    residue_data.round(3).to_csv(out_dir + os.sep + pdb_id + '.csv')
    #nx.write_graphml(residue_graph, out_dir + os.sep + pdb_id + '.graphml')
    return out_dir + os.sep + pdb_id + '.csv'

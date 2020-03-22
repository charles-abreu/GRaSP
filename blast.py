import os
import warnings
from glob import glob
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

##########################################################
# Gets the fasta sequence from structure file
# In: .pdb file
# Out: respective .fasta file for each chain in .pdb file
##########################################################
def get_fasta(file_name, out_dir):
    # get pdb_id in file name
    pdb_id = os.path.basename(file_name).split('.')[0]
    structure_len = 0
    # For each chain crate a .fasta file
    for record in SeqIO.parse(file_name, "pdb-atom"):
        record.id = pdb_id[:4] + record.annotations['chain']
        structure_len += len(record)
        SeqIO.write(record, out_dir + os.sep + record.id  + '.fasta', "fasta")
    return structure_len
###################################################
# Run blast for all .fasta data in project path
###################################################
def run_blast(pdb_id, fasta_path, template_path):
    fasta_list = glob(fasta_path + os.sep + pdb_id + '*.fasta')

    num_templates = 0
    e_value = 0.001

    while num_templates == 0:
        e_value *= 10

        if e_value > 100:
            break

        for fasta_file in fasta_list:
            out_file = fasta_file.replace('.fasta', '.out')
            # BLAST parametres
            blast = './blast/ncbi-blast-2.9.0+/bin/blastp'
            db = ' -db biolip/biolip_blast/biolip_organic.fasta' #
            ev = ' -evalue ' + str(e_value)
            query = ' -query ' + fasta_file
            outfmt = ' -outfmt "6 delim=,"'
            result = ' -out ' + out_file
            # Run blast
            os.system( blast + db + ev + query + outfmt + result)

        template_list = get_blast_result(pdb_id, fasta_path, template_path)
        num_templates = len(template_list)
    return template_list
##########################################
# GET CANDIDATES
##########################################
#Fields: query acc., subject acc., evalue, q. start, q. end, s.start, s. end
def get_cadidate_list(result_file, pdb_id):
    out_file = open(result_file)

    pdb_list = set()
    # Adicionando o pdb sem dividir por cadeia
    for line in out_file:
        pdb_chain = line.split(',')[1]
        #percent = float(line.split(',')[2])

        if pdb_chain[:4] != pdb_id[:4]:
            #if percent > 30.0:
            pdb_list.add(pdb_chain)

    #if len(pdb_list) > 100:
    #    return list(pdb_list)[:100]
    #else:
    return pdb_list

##########################################################
# Gets the blast results of the project and save in the .dat file
# In: Path with .out files
# Out: .dat file containing the structure templates
##########################################################
def get_blast_result(pdb_id, data_path, template_path):
    out_list = glob(data_path + os.sep + pdb_id + '*.out')
    # pdb list in string format
    template_set = set()
    for blast_result in out_list:
        template_set = template_set | get_cadidate_list(blast_result, pdb_id)

    #out_name = template_path + os.sep + pdb_id + '.tpl'
    #with open(out_name, 'w') as out_file:
    #    out_file.write(','.join(template_set))

    return template_set
##########################################################
# Get the num of templates for a specific pdb_id
# In: pdb_id, Path with .out files
# Out: len of template list
##########################################################
def get_num_templates(pdb_id, out_path, template_path):
    out_list = glob(out_path + os.sep + pdb_id + '*.out')
    # pdb list in string format
    template_set = set()
    for blast_result in out_list:
        template_set = template_set | get_cadidate_list(blast_result, pdb_id)

    return len(template_set)

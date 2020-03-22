# My modules
import naccess
import blast
import modeling
import prediction
#
import sys
import os
import time
from glob import glob

def main(pdb_dir, project_name):
    project_path = 'projects' + os.sep + project_name + os.sep
    pdb_list = glob(pdb_dir + os.sep + '*.pdb')
    log_file = open(project_path + 'log.csv', 'w')
    log_file.write('pdb_id,seq_len,time\n')
    for pdb_file in pdb_list:
        ####################
        inicio = time.time()
        ####################
        pdb_id = os.path.basename(pdb_file).split('.')[0]
        # Compute acesible surface area
        rsa_file = naccess.run_naccess(pdb_file, project_path + 'rsa_files/')
        len_seq = blast.get_fasta(pdb_file, project_path + 'fasta_files/')
        template_list = blast.run_blast(pdb_id, project_path + 'fasta_files/', project_path + 'templates/')
        if template_list:
            #template_file = blast.get_blast_result(pdb_id, project_path + 'fasta_files/', project_path + 'templates/')
            csv_file = modeling.compute_matrix(pdb_file, rsa_file, project_path + 'matrix/')
            prediction.class_experiment(csv_file, template_list, '../experiments/biolip/data/', project_path + 'predictions/')
        else:
            print(pdb_id + ': Templates not found!')
        ###################
        fim = time.time()
        ###################
        print('Done!')
        tempo = round(fim - inicio, 2)
        log_file.write(pdb_id + ',' + str(len_seq) + ',' + str(tempo) + '\n')

    log_file.close()
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

from glob import glob
import os
######################################
# Executando o NACCESS
######################################
def run_naccess(file_name, out_dir):
    pdb_id = os.path.basename(file_name).split('.')[0]
    os.system('./../naccess/naccess -h ' + file_name)
    os.system('mv ' + pdb_id + '.rsa ' + out_dir)
    os.system('rm ' + pdb_id + '.log')
    os.system('rm ' + pdb_id + '.asa')
    # Return the name of .rsa file
    return out_dir + pdb_id + '.rsa'

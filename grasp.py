import os, sys, getopt, time, warnings
from glob import glob
from app.exposure import Exposure
from app.blast import BLAST
from app.protein import Protein
from app.modeling import Properties, Interactions, Layer
from app.prediction import BalancedPrediction
from Bio import BiopythonWarning

basedir = os.path.dirname(__file__)
config = {
    'BLAST_DIR' : os.path.join(basedir,'app/data/'),
    'DATA_DIR' : os.path.join(basedir,'app/data/'),
    'TEMPLATES_DIR': os.path.join(basedir,'app/data/templates/')
}

warnings.simplefilter('ignore', BiopythonWarning)

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hp:o:n:", ['help', 'protein=',
                                    'output=','naccess='])
    except getopt.GetoptError:
        print ('grasp.py -p <protein_data> -o <out_dir>')
        sys.exit(2)

    protein_file = ''
    out_dir = ''
    naccess_dir = False

    for opt, arg in opts:
        if opt == '-h':
            print ('grasp.py -d <dataset_dir> -o <out_dir>')
            sys.exit()
        elif opt in ("-p", "--protein"):
            protein_file = arg
        elif opt in ("-o", "--output"):
            out_dir = arg
        elif opt in ("-n", "--naccess"):
            naccess_dir = arg

    # Suport both file and directory with pdb files
    if not os.path.isdir(protein_file):
        pdb_list = [protein_file]
    else:
        pdb_list = glob(protein_file + os.sep + '*.pdb')

    for pdb in pdb_list:
        ####################
        inicio = time.time()
        ####################
        this_protein = Protein(pdb)
        print('--------------------------------------')
        print('Start prediction to protein ' + this_protein.pdb_id)
        print('--------------------------------------')
        # STEP 1:  Compute asa or exposure
        if naccess_dir:
            if not os.path.exists(naccess_dir):
                print( 'WARNING: naccess file not found!')
            print('Computing asa information to ' + this_protein.pdb_id)
            Exposure().run_naccess(this_protein, naccess_dir)
        else:
            print('Computing exposure information to ' + this_protein.pdb_id)
            Exposure().compute_hse(this_protein)

        # STEP 2: Get templates using BLAST
        this_protein.get_fasta(out_dir)
        templates = BLAST(config['BLAST_DIR']).run(this_protein)
        this_protein.set_templates(templates)

        if templates:
            # STEP 3: Compute properties
            print('Computing residue properties for ' + this_protein.pdb_id)
            p = Properties(config['DATA_DIR'])
            p.compute_properties(this_protein)
            # STEP 4: Compute Interactions
            print('Computing Interactions for ' +  this_protein.pdb_id)
            i = Interactions(config['DATA_DIR'])
            i.compute_interactions(this_protein)
            # STEP 5: Compute layers
            print('Computing layers properties for ' +  this_protein.pdb_id)
            l = Layer()
            l.compute_layer(this_protein, 2)
            # STEP 6: Prediction
            b = BalancedPrediction(config['TEMPLATES_DIR'])
            if naccess_dir:
                b.balanced_prediction(this_protein, out_dir, True)
            else:
                b.balanced_prediction(this_protein, out_dir, False)
        else:
            print( 'WARNING: Templates not found for ' + this_protein.pdb_id)

        ###################
        fim = time.time()
        ###################
        print('Done in {} seconds.'.format(round(fim - inicio, 2)))

if __name__ == '__main__':
    main(sys.argv[1:])

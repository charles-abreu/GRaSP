'''
    BLAST module
    Search templates from biolip using BLAST
    Author: Charles Abreu Santana
    Last update: 2020-04-29
'''
import os

class BLAST:
    def __init__(self, blast_dir):
        self.blast_dir = blast_dir
    ###################################################
    # Run blast for all .fasta data in protein object
    # In: Protein instance
    # Out: list with pdb_id from templates
    ###################################################
    def run(self, query):
        # increase evalue to get mor templates
        num_templates = 0
        e_value = 0.001

        while num_templates == 0:
            e_value *= 10

            if e_value > 100:
                break

            for fasta in query.fasta_files:
                # File with blast results
                out_file = fasta.replace('.fasta', '.out')
                # BLAST parametres
                #blast_dir = './blast/ncbi-blast-2.9.0+/bin/'
                blast_type = 'blastp'
                db = ' -db Grasp/data/biolip_blast/biolip_organic.fasta'
                ev = ' -evalue ' + str(e_value)
                q = ' -query ' + fasta
                outfmt = ' -outfmt "6 delim=,"'
                result = ' -out ' + out_file
                # Run blast
                os.system( './' + self.blast_dir + blast_type + db + ev + q + outfmt + result)

            num_templates = len(BlastResult(out_file))
        return self.get_templates(query)

    ##########################################################
    # Gets the blast results of the project and save in the .dat file
    # In: Path with .out files
    # Out: .dat file containing the structure templates
    ##########################################################
    def get_templates(self, query):
        out_list = [file.replace('.fasta', '.out') for file in query.fasta_files]
        # pdb list in string format
        template_set = set()
        for blast_result in out_list:
            template_set = template_set | BlastResult(blast_result).get_result()
        return template_set

class BlastResult:
    '''
    Fields:
    query acc., subject acc., evalue,
    q. start, q. end, s.start, s. end
    '''
    def __init__(self, file_name):
        self.pdb_set= set()
        with open(file_name) as file:
            for line in file:
                # pdb_id + chain
                self.pdb_set.add(line.split(',')[1])

    def get_result(self):
        return self.pdb_set

    def __len__(self):
        return len(self.pdb_set)

if __name__ == '__main__':
    print('')

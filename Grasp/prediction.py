'''
    Prediction module
    Use machine learning to predict biding site residues
    Author: Charles Abreu Santana
    Last update: 2020-04-29
'''
from sklearn.ensemble import ExtraTreesClassifier
import pandas as pd
import numpy as np
import os

class BalancedPrediction:
    def __init__(self, template_path):
        self.template_path = template_path

    def get_model(self, train_data):
        #class_weight ='balanced' or 'balanced_subsample' max_features = 10
        clf = ExtraTreesClassifier(n_estimators=100, max_depth=None,
            min_samples_split=2, bootstrap = True, random_state=0, n_jobs=2)
        # matrix data (except class), class column
        clf = clf.fit(train_data.iloc[:, :-1], train_data.iloc[:,-1].astype('int'))
        return clf

    def set_train_data(self, protein, naccess):
        df_list = []
        for pdb in protein.templates:
            file_name = self.template_path + os.sep + pdb + '.csv.zip'
            if os.path.exists(file_name):
                df_list.append(pd.read_csv(file_name, index_col = 'res_name'))

        self.train_set = pd.DataFrame()

        if df_list:
            self.train_set = pd.concat(df_list, axis = 0)

            if not naccess:
                # drop rsa columns from templates
                drop_list = ['acc_rel','acc_side','acc_main','acc_apolar','acc_polar',
                'N1_acc_rel','N1_acc_side','N1_acc_main','N1_acc_apolar','N1_acc_polar',
                'N2_acc_rel','N2_acc_side','N2_acc_main','N2_acc_apolar','N2_acc_polar']
                self.train_set = self.train_set.drop(drop_list, axis=1)
            else:
                # drop exposure columns from templates
                drop_list = ['hseu','hsed','N1_hseu','N1_hsed','N2_hseu','N2_hsed']
                self.train_set = self.train_set.drop(drop_list, axis=1)

            # Removendo residuos enterrados
            '''
            if 'acc_rel' in self.train_set.columns:
                buried = self.train_set[self.train_set['acc_rel'] < 0.1]
                buried_list = buried.index
                self.train_set = self.train_set.drop(buried_list, axis=0)
            '''
    def split_data(self, train_data):
        # rows with positive class
        pos_data = train_data[train_data['class'] == 1]
        # rows with negative class
        neg_data = train_data[train_data['class'] == 0]
        # The number of partitions is set according the positive class length
        proportion = int(len(neg_data)/len(pos_data))
        if proportion < 4:
            n_parts = proportion
        else:
            n_parts = int((len(neg_data)/len(pos_data))/4)

        return pos_data, neg_data, n_parts

    # Faz a predição a partir do vetor de probabilidades de um classificador
    def predict_prob(self,prob_vector, threshold):
        pred_vector = []
        for prob in prob_vector:
            if prob < threshold:
                pred_vector.append(0)
            else:
                pred_vector.append(1)
        return pred_vector
    # voting using calssification matrix
    def voting(self, binary_matrix, treshold):
        binary_matrix = binary_matrix.T
        confidence = []
        for i in range(binary_matrix.shape[0]):
            # computa a porcentagem de classificadores na votacao
            confidence.append(np.sum(binary_matrix[i])/float(binary_matrix[i].shape[0]))#
        return self.predict_prob(confidence, treshold)

    def balanced_prediction(self, protein, out_dir, naccess):
        # search and store templates matricies
        self.set_train_data(protein, naccess)
        pos_data, neg_data, N_PARTITIONS = self.split_data(self.train_set)
        # shuffle negative index
        permuted_indices = np.random.permutation(len(neg_data))
        # matrix with class values foa all predictions
        # each line is a prediction of one ensenble
        class_matrix = []
        for i in range(N_PARTITIONS):
            # Concat positive and a fragment os negative instances
            final_matrix = pd.concat([pos_data, neg_data.iloc[permuted_indices[i::N_PARTITIONS]]])
            class_model = self.get_model(final_matrix)
            # probability_vector
            probs = class_model.predict_proba(protein.matrix)[:,1]
            # voting probabilities
            class_matrix.append(self.predict_prob(probs, 0.5))
            # cleaning memory
            del class_model
        # Fzendo predicao final usando a combinacao dos ensembles
        vector_pred = self.voting(np.array(class_matrix), 0.5)
        protein.matrix['prediction'] = vector_pred
        protein.matrix.to_csv(out_dir + os.path.basename(protein.file_name).replace('.pdb', '.csv'), columns=['prediction'])

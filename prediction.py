# -*- coding: utf-8 -*-
import numpy as np
import os
import pandas as pd
from sklearn import preprocessing
from sklearn import metrics
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_validate
from sklearn.neural_network import MLPClassifier
# Oversampling
from collections import Counter
import sys
import math
from sklearn.utils import shuffle
# saving model
from glob import glob
#from sklearn.externals import joblib

def run_model(data):
    #class_weight = 'balanced'
    #class_weight='balanced_subsample'
    #max_features = 10
	clf = ExtraTreesClassifier(n_estimators=100, max_depth=None, min_samples_split=2, bootstrap = True, random_state=0, n_jobs=2)
	# matrix data (except column), class column
	clf = clf.fit(data.iloc[:, :-1], data.iloc[:,-1].astype('int'))
	return clf

def balanced_prediction(data , query_data, out_dir):

    # rows with positive class
    pos_data = data[data['class'] == 1]
    # rows with negative class
    neg_data = data[data['class'] == 0]
    # The number of partitions is set according the positive class length
    proportion = int(len(neg_data)/len(pos_data))
    if proportion < 4:
        N_PARTITIONS = proportion
    else:
        N_PARTITIONS = int((len(neg_data)/len(pos_data))/4)

    neg_data = shuffle(neg_data)
    # The number of partitions is set according the positive class length
    permuted_indices = np.random.permutation(len(neg_data))
    # Negative partitions
    neg_data_list = []

    for i in range(N_PARTITIONS):
       neg_data_list.append(neg_data.iloc[permuted_indices[i::N_PARTITIONS]])

    classification_matriz = []

    for i, neg in enumerate(neg_data_list):
        # Concatenate positive class with negative partition
        final_matrix = pd.concat([pos_data, neg])
        #final_matrix = shuffle(final_matrix)

        classification_model = run_model(final_matrix)
        classification_matriz.append(classification_model.predict_proba(query_data)[:,1])
        del classification_model

    # Preenchendo matriz com as predições de cada ensemble
    aux_matrix = []
    for row in classification_matriz:
        aux_matrix.append(prediction(row, 0.5))

    aux_matrix = np.array(aux_matrix)
    aux_matrix = aux_matrix.T

    # Fzendo predicao final usando a combinacao dos ensembles
    vector_pred = []

    for i in range(aux_matrix.shape[0]):
        # computa a porcentagem de classificadores na votacao
        confidence = np.sum(aux_matrix[i])/float(aux_matrix[i].shape[0])#
        if confidence < 0.5:
            vector_pred.append(0)
        else:
            vector_pred.append(1)

    query_data['prediction'] = vector_pred
    query_data.to_csv(out_dir, columns=['prediction'])
    #result_data = pd.merge(query_data['class'],query_data['prediction'], left_index=True, right_index=True)

    ###################################################
    # Usind average probability
    ###################################################
    #query_data['average'] = np.array(classification_matriz).mean(axis = 0)
    #result_data = pd.merge(result_data,query_data['average'], left_index=True, right_index=True)
    #result_data['average_pred'] = prediction(result_data['average'], 0.5)

    ###################################################
    # Usind max probability
    ###################################################
    #query_data['max'] = np.array(classification_matriz).max(axis = 0)
    #result_data = pd.merge(result_data,query_data['max'], left_index=True, right_index=True)
    #result_data['max_pred'] = prediction(result_data['max'], 0.5)

    # Save matrix data
    #result_data.to_csv(out_dir)
    #return compute_result(query_data['class'],query_data['prediction'])

# Faz a predição a partir do vetor de probabilidades de um classificador
def prediction(probability_vector, threshold):
	pred_vector = []

	for prob in probability_vector:
		if prob < threshold:
			pred_vector.append(0)
		else:
			pred_vector.append(1)

	return pred_vector



#########################################
# Obtem lista com pdbids que serao sados como Templates
# para gerar os modelo de predicao
#########################################
def get_train_list(template_file):
    with open(template_file) as in_file:
        return in_file.read().split(',')

#######################################
# Get Templates and Build the train set
#######################################
def get_train_set(pdb_list, template_path):
    df_list = []
    for pdb in pdb_list:
        file_name = template_path + os.sep + pdb + '.csv'
        if os.path.exists(file_name):
            df_list.append(pd.read_csv(file_name, index_col = 'res_name'))

    new_db = pd.DataFrame()

    if df_list:
        new_db = pd.concat(df_list, axis = 0)

        # Removendo residuos enterrados
        buried = new_db[new_db['acc_rel'] < 0.1]
        buried_list = buried.index

        new_db = new_db.drop(buried_list, axis=0)

    return new_db, len(df_list)

def class_experiment(csv_file, template_list, template_path, out_dir):
  # listca com os templates que serao utilizados para construcao do modelo preditivo
  pdb_id = os.path.basename(csv_file).split('.')[0]
  #template_list = get_train_list(template_file)
  if len(template_list) > 0:
    # obtendo dados de entrada:
    train_data, num_templates = get_train_set(template_list, template_path)

    if train_data.empty:
      print('No templates found')
      return
    # Conjunto de teste
    test_data = pd.read_csv(csv_file, index_col = 'res_name')

    #columns_to_remove = ['acc_side','acc_main','acc_apolar','acc_polar',
    #                    'N1_acc_side','N1_acc_main','N1_acc_apolar','N1_acc_polar',
    #                    'N2_acc_side','N2_acc_main','N2_acc_apolar','N2_acc_polar']

    #train_data = train_data.drop(columns_to_remove, axis=1)
    #test_data = test_data.drop(columns_to_remove, axis=1)
    out_file = out_dir + os.sep + pdb_id + '.csv'
    balanced_prediction(train_data , test_data, out_file)

import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, roc_curve, auc

from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neural_network import MLPClassifier

import sys

#from sklearn.model_selection import RandomizedSearchCV
#from sklearn.model_selection import GridSearchCV
#from sklearn.model_selection import train_test_split


my_dat = sys.argv[1]


def run_acc(clf, method_str):
	clf.fit(X_train, y_train)
	pred = clf.predict(X_valid)
	scores = clf.predict_proba(X_valid)[:,1]
	fpr, tpr, thresholds = roc_curve(y_valid, scores)
	acc = accuracy_score(y_valid, pred)
	# make a data frame from all of this
	auc_df = pd.DataFrame(data={"threshold": thresholds, \
		"tpr": tpr, "fpr": fpr})
	# write it out
	auc_df['acc']= acc
	auc_df['auc'] = auc(fpr, tpr)
	auc_df['method'] = method_str
	print(f'{method_str} has an accuracy of {acc}')
	return auc_df

X_train = np.loadtxt("data/05_train_df/%s_x_train.csv" %my_dat, skiprows=1, delimiter=",")
Y_train = np.loadtxt("data/05_train_df/%s_y_train.csv" %my_dat, skiprows=1, delimiter=",")

# do a different split to try these
X = X_train
y = Y_train
#X_train, X_test, y_train, y_test = train_test_split(X, y, \
# test_size=0.2, stratify = y, random_state=28)

fold_brk = pd.read_csv("data/05_train_df/%s_folds.csv" %my_dat)

for fold in range(1,7):
	train_idx = [x!=fold for x in fold_brk['partition'].tolist()]
	valid_idx = [x==fold for x in fold_brk['partition'].tolist()]

	X_train = X[train_idx, :]
	y_train = y[train_idx]
	X_valid = X[valid_idx,:]
	y_valid = y[valid_idx]

	# RF
	rf_acc = run_acc(RandomForestClassifier(), "RandomForest")

	# SVM - RBF
	svm_acc = run_acc(svm.SVC(probability=True), "SVC-RBF")

	# NB
	nb_acc = run_acc(GaussianNB(), "Naive Bayes")

	# KNN
	knn_acc = run_acc(KNeighborsClassifier(), "KNN")


	# GLMNET
	lasso_acc = run_acc(LogisticRegression(penalty='l1', \
	 tol=0.01, solver='saga'), "Lasso")
	rr_acc = run_acc(LogisticRegression(penalty='l2', tol=0.01, \
		solver='saga'), "Ridge")
	en_acc = run_acc(LogisticRegression(penalty='elasticnet', \
	 tol=0.01, l1_ratio=0.5, solver='saga'), "ElasticNet")

	# DT
	dt_acc = run_acc(DecisionTreeClassifier(), "DecisionTree")

	# LDA
	lda_acc = run_acc(LinearDiscriminantAnalysis(), "LDA")

	# Boosting
	adab_acc = run_acc(AdaBoostClassifier(), "AdaBoost")

	# MLP - slow
	mlp_acc = run_acc(MLPClassifier(), "Multi-Layer Perceptron")


	frames = [rf_acc, svm_acc, nb_acc, knn_acc, lasso_acc, rr_acc, \
	en_acc, dt_acc, lda_acc, adab_acc, mlp_acc]
	df_result = pd.concat(frames)
	df_result.to_csv("%s_alt_class_%s.csv" %(my_dat, fold))


# hyperparameter tuning
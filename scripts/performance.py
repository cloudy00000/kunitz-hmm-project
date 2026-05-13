#!/usr/bin/env python3

'''
This script is for the optimization
This script evaluates binary classification performance for a given E-value threshold.
Usage: python preformance.py <predictions_file> <threshold> 

Input file format:
protein_id    evalue    true_label
Where true_label is 1 for positive and 0 for negative.
#For example, you can run the following command on terminal to evaluate the performance
#for i in $(seq 1 10); do python3 preformance.py blast_predictions_$i.txt 0.001; done
#for i in $(seq 1 10); do python3 preformance.py blast_prediction.txt "1e-$i" ; done
#This will run the script for each of the prediction files and print the accuracy and MCC for each threshold. You can adjust the threshold as needed to find the best performance.
'''

import sys
import numpy as np

def get_predictions(fname):
    '''Reads the predictions from the given file and returns a list of tuples (id, score).
    Expects the file to have three columns: id, score, and label.
    true_label is 1 for positive and 0 for negative.'''    
    preds = []
    fh = open(fname)
    for line in fh:
        v = line.strip().split()
        if len(v) != 3:
            continue
        protein_id = v[0]
        evalue = float(v[1])
        true_label = int(v[2])
        preds.append((protein_id, evalue, true_label))
    return preds  

def get_confusion_matrix(preds, threshold):
    '''Returns the confusion matrix for the given predictions and threshold.
    if evalue <= threshold, predict positive (1), else predict negative (0).'''    
    confusion_matrix = np.zeros((2, 2), dtype=int)
    for protein_id, score, true_label in preds:
        predicted_label = 0
        if score <= threshold:
            predicted_label = 1
        confusion_matrix[true_label][predicted_label] += 1
    return confusion_matrix

def get_accuracy(confusion_matrix):
    '''Returns the accuracy for the given confusion matrix.'''
    return (confusion_matrix[0][0] + confusion_matrix[1][1]) / np.sum(confusion_matrix)

def get_mcc(confusion_matrix):
    '''Returns the Matthews correlation coefficient for the given confusion matrix.'''
    tp = confusion_matrix[1][1]
    tn = confusion_matrix[0][0]
    fp = confusion_matrix[0][1]
    fn = confusion_matrix[1][0]

    numerator = tp * tn - fp * fn
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    return numerator / denominator if denominator != 0 else 0            

def get_metrics(confusion_matrix):
    '''Returns TP, FP, FN, TN  the accuracy and MCC for the given confusion matrix.'''
    tp = confusion_matrix[1][1]
    tn = confusion_matrix[0][0]
    fp = confusion_matrix[0][1]
    fn = confusion_matrix[1][0] 

    Sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    Specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1_score = 2 * (precision * Sensitivity) / (precision + Sensitivity) if (precision + Sensitivity) > 0 else 0

    return  Sensitivity, Specificity, fpr, precision, f1_score

if __name__ == "__main__":
    fname = sys.argv[1]
    threshold = float(sys.argv[2])

    preds = get_predictions(fname)

    confusion_matrix = get_confusion_matrix(preds, threshold)

    accuracy = get_accuracy(confusion_matrix)
    mcc = get_mcc(confusion_matrix)
    sensitivity, specificity, fpr, precision, f1_score = get_metrics(confusion_matrix)

    print(f"Threshold: {threshold}, Accuracy: {accuracy:}, MCC: {mcc:}")
    print(f"Sensitivity: {sensitivity:}, Specificity: {specificity:}, FPR: {fpr:}, Precision: {precision:}, F1 Score: {f1_score:}")
    print(f"Confusion Matrix:\n{confusion_matrix}")
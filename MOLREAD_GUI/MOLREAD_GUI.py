# -*- coding: utf-8 -*-
"""
"""
import tkinter as tk
from GUI_class import *
from rdkit import Chem
from tkinter.filedialog import askopenfilename
import pandas as pd
import os
import numpy as np
from keras.models import load_model
import time
from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys, DataStructs

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = ""

# LOAD MODEL
model = load_model("models/v1.h5") # Load the model here.

def convert_time(second):
    day = second/86400
    hour = (day - int(day))*24
    minute = (hour - int(hour))*60
    second = round((minute - int(minute))*60,4)
    return(str(int(day)) + ' DAYS: '+ str(int(hour)) + ' HOURS: '+ str(int(minute)) + ' MINUTES: ' + str(second) + ' SECONDS')

def generate_MACCS(smiles):
    header = ['bit' + str(i) for i in range(167)]
    data = []
    for i in range(len(smiles)):
        mol = Chem.MolFromSmiles(smiles[i])
        ds = list(MACCSkeys.GenMACCSKeys(mol).ToBitString())
        data.append(ds)
    return data,header

def decode_message(data,header):
    df = pd.DataFrame(data,columns=header)
    message_predictions = model.predict(df.values)
    message_predictions_list = []
    for i in range(len(message_predictions)):
        message_predictions_list.append(np.argmax(message_predictions[i]))
    return message_predictions_list

def label_switching_decoder(key_smiles, bit_list, nmol_df):
    '''
    :param key_smiles: key molecules
    :param bit_list: model predictions
    :param df: df where to pick key molecule and the 'neighbor' molecules
    :return: list; ACSII code
    '''
    bit_list = list(map(int, bit_list)) #conver string to integers

    # build a list from 0 to 127
    orig_label = [i for i in range(128)]

    key_mol = Chem.MolFromSmiles(key_smiles)
    key_fp = MACCSkeys.GenMACCSKeys(key_mol)
    # rebuild root_seed and rotor_seed based on MW and number of atoms of key_mol
    root_seed = int(Chem.Descriptors.ExactMolWt(key_mol))
    rotor_seed = key_mol.GetNumAtoms()

    #pick 128 neighbor molecules
    # Pick the 128 reference molecules
    np.random.seed(root_seed)
    ref_smiles = np.random.choice(nmol_df.SMILES, size=128, replace=False)
    #compute the distance
    dist = []
    for i in range(len(ref_smiles)):
        mol = Chem.MolFromSmiles(ref_smiles[i])
        fp = MACCSkeys.GenMACCSKeys(mol)
        dist.append(DataStructs.FingerprintSimilarity(key_fp, fp))

    decoded_message = []
    for index, bit in enumerate(bit_list):
        SEED = root_seed + index * rotor_seed
        # Pick the 128 reference molecules
        np.random.seed(SEED)
        step_dist = np.random.choice(dist, size=len(dist), replace=False)
        # Base on the distance, swap the original cluster labels

        # get the index of ordered distances
        dict_rank = [0] * len(dist)
        for i, x in enumerate(sorted(range(len(step_dist)), key=lambda y: step_dist[y])):
            dict_rank[x] = i
        swaper_dict = dict(zip(orig_label, dict_rank))
        # print(swaper_dict)
        decoded_message.append(swaper_dict.get(bit))
        output = ''.join([chr(i) for i in decoded_message])
    return output

def placeholder_call_back(gui):

    neighbor_mol = pd.read_csv('Libraries/V1/nmol.csv', index_col='ID')

    #Get input data from the user
    encrypted_msg = gui.send_data()

    molecular_key = gui.get_molecular_key()
    smiles = encrypted_msg.rstrip().split('.')


    #Decode and update the output message box
    data,header = generate_MACCS(smiles)
    print('Decoding...')
    start_time = time.time()
    decrypted_msg = decode_message(data,header)


    decoded_msg = label_switching_decoder(molecular_key, decrypted_msg, neighbor_mol)
    duration = convert_time(time.time() - start_time)
    todisplay = 'Time Elapsed: ' + str(duration)
    print(todisplay)

    gui.update_decoded(decoded_msg)
    gui.display_run_info(todisplay)


def main():

    root = tk.Tk()
    gui = GUI(root)
    gui.start(lambda: placeholder_call_back(gui))   # call input message function here
    root.mainloop()


if __name__ == '__main__':
    main()

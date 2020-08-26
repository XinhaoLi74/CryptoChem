import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox as msg

import pandas as pd
import numpy as np
import random
import math
import time
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs


############ Functions #######################
def eudis(v1, v2):
    dist = [(a - b)**2 for a, b in zip(v1, v2)]
    dist = math.sqrt(sum(dist))
    return dist

def convert_time(second):
    day = second/86400
    hour = (day - int(day))*24
    minute = (hour - int(hour))*60
    second = round((minute - int(minute))*60,4)
    return(str(int(day)) + ' DAYS: '+ str(int(hour)) + ' HOURS: '+ str(int(minute)) + ' MINUTES: ' + str(second) + ' SECONDS')

def text_to_bits(text):
    Bits = []
    for letter in text:
        # translate text to ASCII code
        Bits.append(str(ord(letter)))
    return Bits

def choose_molkey(nmol_df, random_pick=True):
    if random_pick == True:
        key_smiles = random.choice(nmol_df.SMILES)

    return key_smiles

def getCen(sampler_Fname):
    g = open(sampler_Fname,'r')
    cens = {e[0]:[float(r) for r in e[1:]] for e in [i.rstrip().split(',') for i in g.readlines()]}
    g.close()
    return cens

def ret_anaF(sample_frame, target_fp):
    Bdist = 10
    Fscore = {}
    for i in list(sample_frame.values()):
        dst = eudis(i,target_fp)
        if dst < Bdist:
            Bdist = dst
        Fscore[list(sample_frame.keys())[list(sample_frame.values()).index(i)]] = Bdist
    Fname = list(Fscore.keys())[list(Fscore.values()).index(min(list(Fscore.values())))]
    return Fname

def label_switching_encoder(key_smiles, bits, df, nmol_df):
    '''

    :param bits: a list of ACSII code
    :param df: df where to pick key molecule and the 'pad' molecules
    :return: key molecule and chemical messages
    '''
    # molecular key
    key_mol = Chem.MolFromSmiles(key_smiles)
    key_fp = MACCSkeys.GenMACCSKeys(key_mol)
    # build root_seed and rotor_seed based on MW and number of atoms of key_mol
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

    # build a list from 0 to 127
    orig_label = [i for i in range(128)]
    message_mol_list = []

    for index, bit in enumerate(bits):
        SEED = root_seed + index * rotor_seed
        np.random.seed(SEED)
        # Base on the random seed, swap the original distance.
        step_dist = np.random.choice(dist, size=len(dist), replace=False)

        # get the index of ordered distances
        dict_rank = [0] * len(step_dist)
        for i, x in enumerate(sorted(range(len(step_dist)), key=lambda y: step_dist[y])):
            dict_rank[x] = i
        swaper_dict = dict(zip(dict_rank,orig_label))
        # pick mol from df
        # fix the problem that the original text has the ACSII code larger than 127
        if int(bit) < 128:
            rand_mol = random.choice(df[df.clusters == swaper_dict.get(int(bit))]['smiles'])
        else:
            rand_mol = random.choice(df[df.clusters == swaper_dict.get(int(42))]['smiles']) # use * as replacement
        message_mol_list.append(rand_mol)
    return message_mol_list

################## Molwriter class #########################################
class molwriter(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title('CryptoChem V1.0')
        self.geometry('800x800')

        self.notebook = ttk.Notebook(self)

        #tabs
        molwrite_tab = tk.Frame(self.notebook)

        #################################molwrite_tab###################################

        #init the output dataframe (cipher text)
        self.msg_mol_list = None
        self.key_mol = None

        #upload files for encryption
        self.text_from_file = None
        self.input_file_name = tk.StringVar(value='') # show the file name and path on the screen

        self.input_file_frame = tk.Frame(molwrite_tab)
        self.input_file_frame.pack(padx=15, pady=(15, 0), anchor='w')
        self.input_file_button = tk.Button(self.input_file_frame, text='Select file(s)', command=self.file_upload, anchor='w')
        self.input_file_button.pack(side='left')

        self.file_label = tk.Label(self.input_file_frame, textvariable=self.input_file_name, anchor='w')
        self.file_label.pack(side='right', padx=15)

        # Input plain text
        tk.Label(molwrite_tab, text="Plain Text:").pack(padx=15, pady=(15, 0), anchor='w')
        self.input_plain_frame = tk.Frame(molwrite_tab, borderwidth=1, relief='sunken')
        self.input_plain_frame.pack(padx=15, pady=5)

        # make sure the scroll bar works properly
        self.input_plain_frame.grid_rowconfigure(0, weight=1)
        self.input_plain_frame.grid_columnconfigure(0, weight=1)

        self.input_plain_yscrollbar = tk.Scrollbar(self.input_plain_frame)
        self.input_plain_yscrollbar.grid(row=0, column=1, sticky=tk.N+tk.S)

        self.text_from_keyboard = tk.Text(self.input_plain_frame, width=60, height=8, highlightbackground='#ffffff',
                                          highlightcolor="#7baedc", bg='#ffffff', wrap=tk.WORD, font=("Arial", 14),
                                          yscrollcommand=self.input_plain_yscrollbar.set)
        self.text_from_keyboard.focus_set()
        self.text_from_keyboard.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)

        self.input_plain_yscrollbar.config(command=self.text_from_keyboard.yview)

        # Encryption checkboxes and Buttons
        self.encryption_button_frame = tk.Frame(molwrite_tab)
        self.encryption_button_frame.pack(padx=15, pady=(0, 15), anchor='e')

        self.encryption_button = tk.Button(self.encryption_button_frame, text='Encryption', command=self.click_encryption)
        self.encryption_button.pack(side='right')

        self.save_plain_button = tk.Button(self.encryption_button_frame, text='Save', command=self.click_save_plain)
        self.save_plain_button.pack(side='right', padx=15)

        # show molecular key
        tk.Label(molwrite_tab, text="*******************************************************************************************************************").pack(padx=15, pady=(15, 0), anchor='w')
        tk.Label(molwrite_tab, text="Molecular Key:").pack(padx=15, pady=(5, 0), anchor='w')
        self.molkey_frame = tk.Frame(molwrite_tab, borderwidth=1, relief='sunken')
        self.molkey_frame.pack(padx=15, pady=5)

        self.molkey_onscreen = tk.StringVar(value='Selected Molecular Key')

        self.molkey_label = tk.Label(self.molkey_frame, textvariable=self.molkey_onscreen, anchor='w')
        self.molkey_label.pack(side='right', padx=15)


        # show cipher text (a pandas dataframe)
        tk.Label(molwrite_tab, text="Cipher Text:").pack(padx=15, pady=(15, 0), anchor='w')
        self.cipher_show_frame = tk.Frame(molwrite_tab)
        self.cipher_show_frame.pack(padx=15, pady=5)

        #make sure the scroll bar works properly
        self.cipher_show_frame.grid_rowconfigure(0, weight=1)
        self.cipher_show_frame.grid_columnconfigure(0, weight=1)

        self.cipher_text_xscrollbar = tk.Scrollbar(self.cipher_show_frame, orient=tk.HORIZONTAL)
        self.cipher_text_xscrollbar.grid(row=1, column=0, sticky=tk.E+tk.W)
        self.cipher_text_yscrollbar = tk.Scrollbar(self.cipher_show_frame)
        self.cipher_text_yscrollbar.grid(row=0, column=1, sticky=tk.N+tk.S)


        self.cipher_text_show = tk.Text(self.cipher_show_frame, width=80, height=15, bg="lightgrey", fg="black", pady=10,
                                         wrap = tk.WORD, state=tk.DISABLED,
                                        xscrollcommand=self.cipher_text_xscrollbar.set,
                                        yscrollcommand=self.cipher_text_yscrollbar.set)
        self.cipher_text_show.focus_set()
        self.cipher_text_show.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)

        self.cipher_text_xscrollbar.config(command=self.cipher_text_show.xview)
        self.cipher_text_yscrollbar.config(command=self.cipher_text_show.yview)

        # save button
        self.save_cipher_button_frame = tk.Frame(molwrite_tab)
        self.save_cipher_button_frame.pack(padx=15, pady=(0, 15), anchor='e')

        self.save_cipher_button = tk.Button(self.save_cipher_button_frame, text='Export', command=self.click_export)
        self.save_cipher_button.pack(side='right')


        # # pack all the tabs
        self.notebook.add(molwrite_tab, text="MOlWRITE")
        self.notebook.pack(fill=tk.BOTH, expand=1)

    def file_upload(self):
        selected_files = filedialog.askopenfilename(filetypes=[("Text files","*.txt")])
        try:
            # os.path.split is used to show the file name without file path
            self.input_file_name.set(os.path.split(selected_files)[1])
            with open(selected_files, 'r') as file:
                self.text_from_file = file.read()
        except:
            print("No file exists")

    def molkey_upload(self):
        selected_files = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
        try:
            # os.path.split is used to show the file name without file path
            self.input_file_name.set(os.path.split(selected_files)[1])
            with open(selected_files, 'r') as file:
                self.text_from_file = file.read()
        except:
            print("No file exists")


    def click_encryption(self):
        pd.set_option("display.max_colwidth", 1000) #avoid the df_to_string function turncate the long smiles

        if self.text_from_file == None and len(self.text_from_keyboard.get(1.0, 'end-1c')) == 0:
            msg.showerror("No input", "Please select a way to input")

        elif self.text_from_file != None and len(self.text_from_keyboard.get(1.0, 'end-1c')) == 0:
            try:
                start_time = time.time()
                #molwrtie functions from cryptochem
                bit_list = text_to_bits(self.text_from_file)
                key_mol = choose_molkey(neighbor_mol, random_pick=True)
                self.molkey_onscreen.set(key_mol)
                #label_swapper
                mol_list = label_switching_encoder(key_mol, bit_list, ref_frame, neighbor_mol)

                # AR
                sample_frame = getCen(Fname)
                msg_mols = []

                for i in mol_list:
                    target = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(i))
                    anaF = ret_anaF(sample_frame, target)
                    mol = random.choice(open(directory + anaF, 'r').read().splitlines())
                    msg_mols.append(mol)
                duration = convert_time(time.time() - start_time)
                # save the dataframe
                self.msg_mol_list = '.'.join(msg_mols)
                self.key_mol = key_mol
                #show the dataframe on screen
                self.cipher_text_show.config(state=tk.NORMAL) #enable texting
                self.cipher_text_show.delete(1.0, tk.END) # delete all the text first
                self.cipher_text_show.insert(tk.INSERT, self.msg_mol_list)
                self.cipher_text_show.config(state=tk.DISABLED)

                msg.showinfo("Encryption Successful", 'Running Time: ' + str(duration)) # pop-up window shows running time
            except Exception as e:
                msg.showerror("Translation Failed", str(e))

        elif self.text_from_file == None and len(self.text_from_keyboard.get(1.0, 'end-1c')) > 0:
            try:
                start_time = time.time()
                # molwrtie functions from cryptochem
                bit_list = text_to_bits(self.text_from_keyboard.get(1.0, 'end-1c'))
                key_mol = choose_molkey(neighbor_mol, random_pick=True)
                self.molkey_onscreen.set(key_mol)
                # label_swapper
                mol_list = label_switching_encoder(key_mol, bit_list, ref_frame, neighbor_mol)

                # AR
                sample_frame = getCen(Fname)
                msg_mols = []

                for i in mol_list:
                    target = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(i))
                    anaF = ret_anaF(sample_frame, target)
                    mol = random.choice(open(directory + anaF, 'r').read().splitlines())
                    msg_mols.append(mol)
                duration = convert_time(time.time() - start_time)
                # save the dataframe
                self.msg_mol_list = '.'.join(msg_mols)
                self.key_mol = key_mol
                #show the dataframe on screen
                self.cipher_text_show.config(state=tk.NORMAL) #enable texting
                self.cipher_text_show.delete(1.0, tk.END) # delete all the text first
                self.cipher_text_show.insert(tk.INSERT, self.msg_mol_list)
                self.cipher_text_show.config(state=tk.DISABLED)

                msg.showinfo("Encryption Successful", 'Running Time: ' + str(duration))  # pop-up window shows running time
            except Exception as e:
                msg.showerror("Translation Failed", str(e))
        else:
            msg.showerror("Duplicated Input", "PLease choice only one way to input plain text")

    def click_save_plain(self):
        if len(self.text_from_keyboard.get(1.0, 'end-1c')) > 0:
            try:
                export_file_path = filedialog.asksaveasfilename(defaultextension='.txt',
                                                            filetypes=(('Text files', '*.txt'), ('All files', '*.*')))
                if export_file_path is None:
                    return
                contents = self.text_from_keyboard.get(1.0, tk.END)
                with open(export_file_path, "w") as text_file:
                    text_file.write(contents)
                msg.showinfo("Save Successful", os.path.split(export_file_path)[1] + ' is saved')
            except Exception as e:
                msg.showerror("SavE Failed", str(e))
        else:
            msg.showerror("Empty Input", "Empty Input")

    def click_export(self):
        export_msg_file_path = filedialog.asksaveasfilename(defaultextension='.txt')
        with open(export_msg_file_path, 'w') as msg_file:
            msg_file.write(self.msg_mol_list)
        with open((os.path.split(export_msg_file_path)[0] + '/key_' + os.path.split(export_msg_file_path)[1]), 'w') as key_file:
            key_file.write(self.key_mol)
        msg.showinfo("Save Successful", os.path.split(export_msg_file_path)[1] + ' and its key are saved')

if __name__ == "__main__":
    # Load datasets
    directory = 'Libraries/V1/External/'
    Fname = 'Libraries/V1/internal_centroids'
    ref_frame = pd.read_csv('Libraries/V1/refined_Reference.csv',index_col='IDs')
    neighbor_mol = pd.read_csv('Libraries/V1/nmol.csv', index_col='ID')

    app = molwriter()
    app.resizable(width=False, height=True) #make the window unchangeable
    app.mainloop()

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 01:53:12 2018

"""
#from tkinter import ttk
import tkinter as tk
import tkinter.filedialog
#import PIL
#from PIL import ImageTk

#--------------------------------------------
def append_dicts(base,to_add):
    new_dic = {**base,**to_add}

    return new_dic

def _toint_helper(str_input):
    '''assumes input is a string'''
    output = 0
    if len(str_input) < 1:
        output = 0
    else:
        try:
            output = int(str_input)
        except:
            raise ValueError('''Input "{}" cannot be converted into an integer'''.format(str_input))

    return output

def toint(str_input):
    '''It will check if input is string or tuple or int and return an int of values provided '''

    _type = type(str_input)
    output = None
    if _type == str:
        output = _toint_helper(str_input)
    elif _type == tuple:
        output = tuple(_toint_helper(i) for i in str_input)
    elif _type == int:
        output = str_input
    else:
        print("Values entered must be integers")
        raise ValueError

    return output
#-==============================================================================
class GUI(object):
    def __init__(self,root,window_width=1300,window_height=500):
        self.root = root
        self.mainFrame = tk.Frame(self.root)
        win_size = str(window_width) + 'x' + str(window_height)
        self.root.geometry(win_size)
        self.root.title("MOLREAD")
        self.encrypted_msg  = None
        self.decrypted_msg  = None
        self.MOLECULAR_KEY  = None           #LINE ADDED HERE
        self.USE_TEXTBOX    = 1              #flag to tell if program should use data directly from textbox or from internal varible. 1 for textbox, 0 for internal variable

    def _make_mainframe(self, call_back):
        self._make_inputArea()
        self._make_outputArea(call_back)

        self.mainFrame.pack(padx=10, fill='x')
        self.inputArea.pack(side='left')
        self.outputArea.pack(side='right')

    def _make_inputArea(self):
        self.inputArea   = tk.Frame(self.mainFrame,bd=1,relief='sunken')
        self.inputTitle  = tk.Label(self.inputArea,text="Message")

        self.buttonFrame = tk.Frame(self.inputArea)
        self.keyFrame    = tk.Frame(self.inputArea)         #LINE ADDED HERE
        #===================================
        self.inputBox       = tk.Text(self.inputArea)
        self.orLabel        = tk.Label(self.inputArea,text="*"*120)
        self.keyButton      = tk.Button(self.keyFrame,text="Import Key",command=self._open_key)     #LINE ADDED HERE
        self.varKeyLabel    = tk.StringVar(value='')      #LINE ADDED HERE
        self.keyEntry       = tk.Entry(self.keyFrame,textvariable=self.varKeyLabel)
        #self.keyImportLabel = tk.Label(self.keyFrame,textvariable=self.varKeyLabel)      #LINE ADDED HERE
        self.uploadButton   = tk.Button(self.buttonFrame,text='Upload File',command=self._open_file)
        self.varUploadLabel = tk.StringVar(value="                          ")
        self.uploadEntry    = tk.Entry(self.buttonFrame,textvariable=self.varUploadLabel)
        self.uploadLabel    = tk.Label(self.buttonFrame,textvariable=self.varUploadLabel)

        #==================================+
        ##Actually drawing the items on the screen
        #---Input text area
        self.inputTitle.pack()
        self.inputBox.pack()

        #---File Packing area
        self.buttonFrame.pack()
        #tk.Label(self.buttonFrame,text=' ').pack()
        self.uploadButton.pack(side='left')
        #self.uploadEntry.pack(side='right')
        self.uploadLabel.pack(side='right')

        #---Star area area
        self.orLabel.pack()

        #---Import key area
        self.keyFrame.pack()
        self.keyButton.pack(side='left')
        self.keyEntry.pack(side='right')

    def _make_outputArea(self,call_back):
        self.outputArea   = tk.Frame(self.mainFrame,bd=1,relief='sunken')
        self.outputTitle  = tk.Label(self.outputArea, text="Decoded")

        self.output_buttonFrame = tk.Frame(self.outputArea)
        #=====================================
        self.outputBox    = tk.Text(self.outputArea)
        self.decodeButton = tk.Button(self.output_buttonFrame,text="Decode Message",command=call_back)
        self.saveButton   = tk.Button(self.output_buttonFrame,text="Save Ouptut", command=self._save_file)

        self.varInfoLabel = tk.StringVar(value="")
        self.infoLabel    = tk.Label(self.output_buttonFrame,textvariable=self.varInfoLabel)

        self.outputTitle.pack()
        self.outputBox.pack()
        tk.Label(self.outputArea,text=" ").pack()
        self.infoLabel.pack()
        #tk.Label(self.outputArea,text=" ").pack()
        self.output_buttonFrame.pack()
        self.decodeButton.pack(side='left')
        self.saveButton.pack(side='right')

    def start(self, call_back):
        self._make_mainframe(call_back)

    def _call_decode(self):
        #get the encrypted data and send it to the outputArea class
        self.encrypted_msg = self.__encrypted_msg()

        #Here call the module that decrypts the message:
        decrypted = self.encrypted_msg

        #Finally, put the result into the text field of the output box
        self.outputBox.insert(tk.INSERT,decrypted)

    ############################################################################
    ##======================== API Interface Below =============================
    def display_run_info(self, msg):
        """Display a message into the GUI"""
        #self.varInfoLabel = msg
        self.varInfoLabel.set(str(msg))

    def get_molecular_key(self):
        return self.MOLECULAR_KEY
##        return self.keyEntry.get("1.0",'end-1c')

    def send_data(self):
        """Returns the contents of the text box with the message to be decoded"""
        self.encrypted_msg = self.__encrypted_msg()      #collect input message into variable to send_data
        #data = self.inputBox.get()
        return self.encrypted_msg

    def update_decoded(self,decrypted_msg):

        self.outputBox.delete(1.0,tk.END)
        self.outputBox.insert(tk.INSERT,decrypted_msg)

    ##==========================================================================
    def _open_key(self):
        """Ask user to select a file containing a molecular key"""
        filename = tk.filedialog.askopenfilename(initialdir='.',title='import Key',filetypes = (("text files","*.txt"),("all files","*.*")))
        #self.varKeyLabel.set(filename)

        with open(filename,'r') as f:
            self.MOLECULAR_KEY = f.read()
            self.keyEntry.configure(width=len(self.MOLECULAR_KEY))
            self.keyEntry.delete(0,tk.END)
            self.keyEntry.insert(0,self.MOLECULAR_KEY)

    def _open_file(self):
        """Ask user to select a file to load as input"""
        filename = tk.filedialog.askopenfilename(initialdir='.',title='Open File',filetypes = (("text files","*.txt"),("all files","*.*")))
        #self.varUploadLabel.set(filename)       #updates label with selected file name
        self.uploadEntry.delete(0,tk.END)
        self.uploadEntry.insert(0,filename)

        with open(filename,'r') as f:
            data = f.read()
            self.inputBox.delete(1.0,tk.END)
            if len(data) > 100:
                self.inputBox.insert(tk.INSERT, data[:1600] + "\n-----file continues")
            else:
                self.inputBox.insert(tk.INSERT, data)
            self.encrypted_msg = data
            self.USE_TEXTBOX   = 0

    def _save_file(self):
        """Create a file where to save the decrypted information"""
        filename = tk.filedialog.asksaveasfilename(initialdir = ".",title = "Save Filee",filetypes = (("text files","*.txt"),("all files","*.*")))

        with open(filename,'w') as f:
            f.write(self.__decrypted_msg())

    def __encrypted_msg(self):
        #return self.inputBox.get("1.0",'end-1c')
        if self.USE_TEXTBOX:
            return self.inputBox.get("1.0",'end-1c')
        else:
            self.USE_TEXTBOX = 1            #to allow for reseting
            return self.encrypted_msg
    def __decrypted_msg(self):
        return self.outputBox.get("1.0",'end-1c')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:12:09 2019

@author: akshat
"""
import numpy as np

def unique_chars_iterator(smile):
     """
     Iterate over the characters of a smile string. 
     Note that 'Cl' & 'Br' are considered as one character
     """
     atoms = []
     for i in range(len(smile)):
         atoms.append(smile[i])
     return atoms

def smile_to_hot(smile, largest_smile_len, alphabet, type_of_encoding):
    """
    Go from a single smile string to a one-hot encoding.
    """
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
#    print('ENCODING: ', char_to_int)
    # integer encode input smile
    if type_of_encoding==0:
        for _ in range(largest_smile_len-len(smile)):
            smile+=' ' 
    elif type_of_encoding==1: 
        for _ in range(largest_smile_len-len(smile)):
            smile+=' '    
    elif type_of_encoding==2: 
        for _ in range(largest_smile_len-len(smile)):
            smile+='A'        
        
    integer_encoded = [char_to_int[char] for char in unique_chars_iterator(smile)]

        
    # one hot-encode input smile
    onehot_encoded = list()
    for value in integer_encoded:
    	letter = [0 for _ in range(len(alphabet))]
    	letter[value] = 1
    	onehot_encoded.append(letter)
    return integer_encoded, np.array(onehot_encoded)
    

def multiple_smile_to_hot(smiles_list, largest_smile_len, alphabet, type_of_encoding):
    """
    Convert a list of smile strings to a one-hot encoding
    
    Returned shape (num_smiles x len_of_largest_smile x len_smile_encoding)
    """
    hot_list = []
    for smile in smiles_list:
        _, onehot_encoded = smile_to_hot(smile, largest_smile_len, alphabet, type_of_encoding)
        hot_list.append(onehot_encoded)
    return np.array(hot_list)
        

def hot_to_smile(onehot_encoded, alphabet):
    """
    Go from one-hot encoding to smile string
    """    
    # From one-hot to integer encoding
    integer_encoded = onehot_encoded.argmax(1)
#    print('integer_encoded ', integer_encoded)
    
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
#    print('DECODING: ', int_to_char)
    # integer encoding to smile
    regen_smile = "".join(int_to_char[x] for x in integer_encoded)
    regen_smile = regen_smile.strip()
    return regen_smile


def check_conversion_bijection(smiles_list, largest_smile_len, alphabet):
    """
    This function should be called to check successful conversion to and from 
    one-hot on a data set.
    """
    for i, smile in enumerate(smiles_list):
        _, onehot_encoded = smile_to_hot(smile, largest_smile_len, alphabet, type_of_encoding=0)
        regen_smile = hot_to_smile(onehot_encoded, alphabet)
        # print('Original: ', smile, ' shape: ', len(smile))
        # print('REcon: ', regen_smile , ' shape: ', len(regen_smile))
        if smile != regen_smile:
            print('Filed conversion for: ', smile, ' @index: ', i)
            raise Exception('FAILEDDDD!!!')
    print('All conditions passed!')



#with open('smiles_qm9.txt') as f:
#    content = f.readlines()
#content = content[1:]
#content = [x.strip() for x in content] 
#A = [x.split(',')[1] for x in content]
#
#alphabets = ['N', '1', '(', '#', 'C', '3', '5', 'O', '2', 'F', '=', '4', ')', ' ']
#                         
#data = multiple_smile_to_hot(A, len(max(A, key=len)), alphabets, 0)
                         
#check_conversion_bijection(smiles_list=A, largest_smile_len=len(max(A, key=len)), alphabet=alphabets)
                         
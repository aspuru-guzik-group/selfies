"""
This file is meant to go from various representations to 1HOT and back
"""
import numpy as np
from GrammarVAE_codes import to_one_hot, prods_to_eq
import GrammarVAE_grammar as zinc_grammar

def unique_chars_iterator(smile):
     """
     """
     atoms = []
     for i in range(len(smile)):
         atoms.append(smile[i])
     return atoms



def grammar_one_hot_to_smile(one_hot_ls):
    _grammar = zinc_grammar
    _productions = _grammar.GCFG.productions()
    
    # This is the generated grammar sequence
    grammar_seq = [[_productions[one_hot_ls[index,t].argmax()] 
                        for t in range(one_hot_ls.shape[1])] 
                        for index in range(one_hot_ls.shape[0])]
    #print(grammar_seq)
    smile = [prods_to_eq(prods) for prods in grammar_seq]
    
    return grammar_seq, smile


def smile_to_hot(smile, largest_smile_len, alphabet, type_of_encoding):
    """
    Go from a single smile string to a one-hot encoding.
    """
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
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
        

def hot_to_smile(onehot_encoded,alphabet):
    """
    Go from one-hot encoding to smile string
    """    
    # From one-hot to integer encoding
    integer_encoded = onehot_encoded.argmax(1)
    
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    
    # integer encoding to smile
    regen_smile = "".join(int_to_char[x] for x in integer_encoded)
    regen_smile = regen_smile.strip()
    return regen_smile


def check_conversion_bijection(smiles_list, largest_smile_len):
    """
    This function should be called to check successful conversion to and from 
    one-hot on a data set.
    """
    for i, smile in enumerate(smiles_list):
        _, onehot_encoded = smile_to_hot(smile, largest_smile_len)
        regen_smile = hot_to_smile(onehot_encoded)
#        print('Original: ', smile, ' shape: ', len(smile))
#        print('REcon: ', regen_smile , ' shape: ', len(regen_smile))
#        return
        if smile != regen_smile:
            print('Filed conversion for: ', smile, ' @index: ', i)
            break
    print('All conditions passed!')


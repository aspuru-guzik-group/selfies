import nltk
import pdb
#import zinc_grammar
import numpy as np
import h5py
#import molecule_vae

import GrammarVAE_grammar as zinc_grammar

#f = open('data/250k_rndm_zinc_drugs_clean.smi','r')
#L = []    
#
#count = -1
#for line in f:
#    line = line.strip()
#    L.append(line)      # The zinc data set
#f.close()
#
NCHARS = len(zinc_grammar.GCFG.productions())



def prods_to_eq(prods):
    seq = [prods[0].lhs()]
    for prod in prods:
        if str(prod.lhs()) == 'Nothing':
            break
        for ix, s in enumerate(seq):
            if s == prod.lhs():
                seq = seq[:ix] + list(prod.rhs()) + seq[ix+1:]
                break
    try:
        return ''.join(seq)
    except:
        return ''


def get_zinc_tokenizer(cfg):
    long_tokens = [a for a in list(cfg._lexical_index.keys()) if len(a) > 1]
    replacements = ['$','%','^'] # ,'&']
    assert len(long_tokens) == len(replacements)
    for token in replacements: 
        assert token not in cfg._lexical_index
    
    def tokenize(smiles):
        for i, token in enumerate(long_tokens):
            smiles = smiles.replace(token, replacements[i])
        tokens = []
        for token in smiles:
            try:
                ix = replacements.index(token)
                tokens.append(long_tokens[ix])
            except:
                tokens.append(token)
        return tokens
    
    return tokenize


def to_one_hot(smiles, MaxNumSymbols, check=True):
    """ Encode a list of smiles strings to one-hot vectors """
    assert type(smiles) == list
    prod_map = {}
    for ix, prod in enumerate(zinc_grammar.GCFG.productions()):
        prod_map[prod] = ix
    tokenize = get_zinc_tokenizer(zinc_grammar.GCFG)
    tokens = list(map(tokenize, smiles))
    parser = nltk.ChartParser(zinc_grammar.GCFG)
    parse_trees = [next(parser.parse(t)) for t in tokens]
    productions_seq = [tree.productions() for tree in parse_trees]
    
    #if check:
    #    print(productions_seq)
        
    indices = [np.array([prod_map[prod] for prod in entry], dtype=int) for entry in productions_seq]
    one_hot = np.zeros((len(indices), MaxNumSymbols, NCHARS), dtype=np.float32)
    for i in range(len(indices)):
        num_productions = len(indices[i])
        one_hot[i][np.arange(num_productions),indices[i]] = 1.
        one_hot[i][np.arange(num_productions, MaxNumSymbols),-1] = 1.
    return one_hot



def SizeOneHot(smiles, check=True):
    """ Encode a list of smiles strings to one-hot vectors """
    assert type(smiles) == list
    prod_map = {}
    for ix, prod in enumerate(zinc_grammar.GCFG.productions()):
        prod_map[prod] = ix
    tokenize = get_zinc_tokenizer(zinc_grammar.GCFG)
    tokens = list(map(tokenize, smiles))
    parser = nltk.ChartParser(zinc_grammar.GCFG)
    parse_trees = [next(parser.parse(t)) for t in tokens]
    productions_seq = [tree.productions() for tree in parse_trees]
    
    indices = [np.array([prod_map[prod] for prod in entry], dtype=int) for entry in productions_seq]
    return len(indices[0])


# SINGLE EXAMPLE
#smile = [L[0]]
##smile = ['C']
#one_hot_single =  to_one_hot(smile, )
#print(one_hot_single.shape)
#print(one_hot_single)


# GOING THROUGH ALL OF ZINC....

#OH = np.zeros((len(L),MAX_LEN,NCHARS))
#for i in range(0, len(L), 100):
#    print('Processing: i=[' + str(i) + ':' + str(i+100) + ']')
#    onehot = to_one_hot(L[i:i+100], False)
#    OH[i:i+100,:,:] = onehot
#
#h5f = h5py.File('zinc_grammar_dataset.h5','w')
#h5f.create_dataset('data', data=OH)
#h5f.close()

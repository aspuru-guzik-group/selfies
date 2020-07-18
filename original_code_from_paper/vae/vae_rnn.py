#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 SELFIES: a robust representation of semantically constrained graphs with an
           example application in chemistry (https://arxiv.org/abs/1905.13741)
           by Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich, Alan Aspuru-Guzik


           Variational Auto Encoder (VAE) for chemistry
                  comparing SMILES, DeepSMILES, GrammarVAE and SELFIES
                  representation using reconstruction quality, diversity and
                  latent space validity as metrics of interest
                  v0.1.0 -- 04. August 2019

 information:
     This is the original code used to generate the data in our paper.
     It used a hand-written SELFIES encoding (Table2 of paper), and cannot
     easily be adapted to other situations. If you want to use the VAE, please
     see ..\WithFullEncoder\chemistryVAE.py
     That file is connected with the selfies.encoder/selfies.decoder, and can
     be applied on general datasets. For more documentation, please look there.


 settings*.yml
     these four files contain the settings for values for the best model described in the paper



For comments, bug reports or feature ideas, please send an email to
mario.krenn@utoronto.ca and alan@aspuru.com

"""
import os, sys
sys.path.append('VAE_dependencies')
print(sys.argv)

import numpy as np
import yaml
import torch
from torch import nn
from random import shuffle
from data_loader import multiple_smile_to_hot, grammar_one_hot_to_smile
import pandas as pd


from GPlus2S import GrammarPlusToSMILES, IncludeRingsForSMILES
import deepsmiles
converter = deepsmiles.Converter(rings=True, branches=True) # Coverter object, described by authors

from rdkit.Chem import MolFromSmiles
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

import time

def _make_dir(directory):
    os.makedirs(directory)

def save_models(encoder, decoder, epoch):
    out_dir = './saved_models/{}'.format(epoch)
    _make_dir(out_dir)
    torch.save(encoder, '{}/E'.format(out_dir))
    torch.save(decoder, '{}/D'.format(out_dir))



class VAE_encode(nn.Module):

    def __init__(self, layer_1d, layer_2d, layer_3d, latent_dimension):
        """
        Fully Connected layers for the RNN.
        """
        super(VAE_encode, self).__init__()

        # Reduce dimension upto second last layer of Encoder
        self.encode_4d = nn.Sequential(
            nn.Linear(len_max_molec1Hot, layer_1d),
            nn.ReLU(),
            nn.Linear(layer_1d, layer_2d),
            nn.ReLU(),
            nn.Linear(layer_2d, layer_3d),

			nn.ReLU(),
        )

        # Latent space mean
        self.encode_mu = nn.Linear(layer_3d, latent_dimension)

        # Latent space variance
        self.encode_log_var = nn.Linear(layer_3d, latent_dimension)


    def reparameterize(self, mu, log_var):
        """
        This trick is explained well here:
            https://stats.stackexchange.com/a/16338
        """
        #print('reparameterize(self, mu, log_var)')
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return eps.mul(std).add_(mu)


    def forward(self, x):
        """
        Pass throught the Encoder
        """
        # Go down to dim-4
        h1 = self.encode_4d(x)

        # Go down to dim-2 & produce the mean & variance
        mu = self.encode_mu(h1)
        log_var = self.encode_log_var(h1)

        # Reparameterize
        z = self.reparameterize(mu, log_var)
        return z, mu, log_var



class VAE_decode(nn.Module):

    def __init__(self, latent_dimension, gru_stack_size, gru_neurons_num):
        """
        Through Decoder
        """
        super(VAE_decode, self).__init__()
        self.gru_stack_size = gru_stack_size
        self.gru_neurons_num = gru_neurons_num

        # Simple Decoder
        self.decode_RNN  = nn.GRU(
                input_size  = latent_dimension,
                hidden_size = gru_neurons_num,
                num_layers  = gru_stack_size,
                batch_first = False)

        self.decode_FC = nn.Sequential(
            nn.Linear(gru_neurons_num, len_alphabet),
        )


    def init_hidden(self, batch_size = 1):
        weight = next(self.parameters())
        return weight.new_zeros(self.gru_stack_size, batch_size, self.gru_neurons_num)


    def forward(self, z, hidden):
        """
        A forward pass throught the entire model.
        """
        # Decode
        l1, hidden = self.decode_RNN(z, hidden)
        decoded = self.decode_FC(l1)        # fully connected layer

        return decoded, hidden


def IsCorrectSMILES(smiles):
    try:
        resMol=MolFromSmiles(smiles, sanitize=True)
    except Exception:
        resMol=None

    if resMol==None:
        return 0
    else:
        return 1



def sample_latent_space(latent_dimension):
    model_encode.eval()
    model_decode.eval()

    fancy_latent_point=torch.normal(torch.zeros(latent_dimension),torch.ones(latent_dimension))
    #print('fancy_latent_point: ',fancy_latent_point)
    #print('fancy_latent_point length: ',len(fancy_latent_point))
    hidden = model_decode.init_hidden()
    gathered_atoms = []
    for ii in range(len_max_molec):                 # runs over letters from SMILES (len=size of largest molecule)
        fancy_latent_point = fancy_latent_point.reshape(1, 1, latent_dimension)
        fancy_latent_point=fancy_latent_point.to(device)
        decoded_one_hot, hidden = model_decode(fancy_latent_point, hidden)

        decoded_one_hot = decoded_one_hot.flatten()
        decoded_one_hot = decoded_one_hot.detach()

        soft = nn.Softmax(0)
        decoded_one_hot = soft(decoded_one_hot)
        #print('decoded_one_hot: ', decoded_one_hot)
        #print('decoded_one_hot length: ', len(decoded_one_hot))

        _,MaxIdx=decoded_one_hot.max(0)
        gathered_atoms.append(MaxIdx.data.cpu().numpy().tolist())

    model_encode.train()
    model_decode.train()

    return gathered_atoms


def list_to_molecule_string(mollist,current_alphabet):
    molecule=''
    #print('mollist: ',mollist)
    #print('current_alphabet: ',current_alphabet)
    for ii in mollist:
        molecule+=current_alphabet[ii]

    molecule=molecule.replace(' ','')

    return(molecule)


def latent_space_quality(latent_dimension,sample_num=100):
    total_samples=0
    total_correct=0
    all_correct_molecules=[];
    print('latent_space_quality: sample_num: ',sample_num)
    while total_samples<=sample_num:
        Molecule=''
        while len(Molecule)==0:
            is_decoding_error=0
            if type_of_encoding==0: # SMILES
                Molecule=list_to_molecule_string(sample_latent_space(latent_dimension),' FONC()=#12345')

            if type_of_encoding==1: # DeepSMILES
                Molecule=list_to_molecule_string(sample_latent_space(latent_dimension),' FONC)=#3456789')
                #print(Molecule)
                try:
                    Molecule=converter.decode(Molecule)
                except Exception:
                    is_decoding_error=1
                    Molecule='err'

            if type_of_encoding==2: # GRIN
                Molecule=list_to_molecule_string(sample_latent_space(latent_dimension),'ABCDEFGHIJKLMN')
                Molecule=IncludeRingsForSMILES(GrammarPlusToSMILES(Molecule,'X0'))


            if type_of_encoding==3: # GrammarVAE
                RuleIdx=sample_latent_space(latent_dimension)
                Rule1Hot=np.zeros((1,len_max_molec,len_alphabet ))
                for jj in range(len(RuleIdx)):
                    Rule1Hot[0,jj,RuleIdx[jj]]=1
                gramm, MoleculeL=grammar_one_hot_to_smile(Rule1Hot)
                Molecule=MoleculeL[0]
                #print('grammar: ',gramm)

        total_samples+=1
        if is_decoding_error==0:
            isItCorrect=IsCorrectSMILES(Molecule)
        else:
            isItCorrect=0

        if isItCorrect==1:
            #print('correct: ',Molecule)
            total_correct+=1
            SameMol=0
            for jj in range(len(all_correct_molecules)):
                if Molecule==all_correct_molecules[jj]:
                    SameMol=1
                    break

            if SameMol==0:
                all_correct_molecules.append(Molecule)
                #print('#', len(all_correct_molecules), ': ', Molecule)
        #else:
            #print('wrong OrigMole: ', Molecule)
            #print(str(total_correct)+': wrong molecule: ', Molecule)


    return total_correct, len(all_correct_molecules)





def quality_in_validation_set(data_valid):
    x = [i for i in range(len(data_valid))]  # random shuffle input
    shuffle(x)
    data_valid = data_valid[x]

    quality_list=[]
    for batch_iteration in range(min(25,num_batches_valid)):  # batch iterator

        current_smiles_start, current_smiles_stop = batch_iteration * batch_size, (batch_iteration + 1) * batch_size
        inp_smile_hot = data_valid[current_smiles_start : current_smiles_stop]

        inp_smile_encode = inp_smile_hot.reshape(inp_smile_hot.shape[0], inp_smile_hot.shape[1] * inp_smile_hot.shape[2])
        latent_points, mus, log_vars = model_encode(inp_smile_encode)
        latent_points = latent_points.reshape(1, batch_size, latent_points.shape[1])

        hidden = model_decode.init_hidden(batch_size = batch_size)
        decoded_one_hot = torch.zeros(batch_size, inp_smile_hot.shape[1], inp_smile_hot.shape[2]).to(device)
        for seq_index in range(inp_smile_hot.shape[1]):
            decoded_one_hot_line, hidden  = model_decode(latent_points, hidden)
            decoded_one_hot[:, seq_index, :] = decoded_one_hot_line[0]

        decoded_one_hot = decoded_one_hot.reshape(batch_size * inp_smile_hot.shape[1], inp_smile_hot.shape[2])
        _, label_atoms  = inp_smile_hot.max(2)
        label_atoms     = label_atoms.reshape(batch_size * inp_smile_hot.shape[1])

        # assess reconstruction quality
        _, decoded_max_indices = decoded_one_hot.max(1)
        _, input_max_indices   = inp_smile_hot.reshape(batch_size * inp_smile_hot.shape[1], inp_smile_hot.shape[2]).max(1)

        differences = 1. - torch.abs(decoded_max_indices - input_max_indices)
        differences = torch.clamp(differences, min = 0., max = 1.).double()
        quality     = 100. * torch.mean(differences)
        quality     = quality.detach().cpu().numpy()
        quality_list.append(quality)

    return(np.mean(quality_list))






def train_model(data_train, data_valid, num_epochs, latent_dimension, tensorBoard_graphing, checkpoint, lr_enc, lr_dec, KLD_alpha, sample_num):
    """
    Train the Variational Auto-Encoder
    """

    print('num_epochs: ',num_epochs)

    if tensorBoard_graphing:
        from tensorboardX import SummaryWriter
        writer = SummaryWriter()

    # initialize an instance of the model
    optimizer_encoder = torch.optim.Adam(model_encode.parameters(), lr=lr_enc)
    optimizer_decoder = torch.optim.Adam(model_decode.parameters(), lr=lr_dec)

    data_train=torch.tensor(data_train, dtype=torch.float)
    data_train=data_train.to(device)

    #print(data)
    quality_valid_list=[0,0,0,0];
    for epoch in range(num_epochs):
        x = [i for i in range(len(data_train))]  # random shuffle input
        shuffle(x)

        #B = [data[iai] for ii in x]         # Shuffled inputs (TODO: unnecesary variable)
        data_train  = data_train[x]
        start = time.time()
        for batch_iteration in range(num_batches_train):  # batch iterator

            loss, recon_loss, kld = 0., 0., 0.

            current_smiles_start, current_smiles_stop = batch_iteration * batch_size, (batch_iteration + 1) * batch_size
            inp_smile_hot = data_train[current_smiles_start : current_smiles_stop]

#            print('data', data.shape)
#            print('INP_SMILE_HOT', inp_smile_hot.shape)
            inp_smile_encode = inp_smile_hot.reshape(inp_smile_hot.shape[0], inp_smile_hot.shape[1] * inp_smile_hot.shape[2])
            latent_points, mus, log_vars = model_encode(inp_smile_encode)
            latent_points = latent_points.reshape(1, batch_size, latent_points.shape[1])
#            print('LATENT_POINTS', latent_points.shape)

            kld += -0.5 * torch.mean(1. + log_vars - mus.pow(2) - log_vars.exp())

            hidden = model_decode.init_hidden(batch_size = batch_size)
            decoded_one_hot = torch.zeros(batch_size, inp_smile_hot.shape[1], inp_smile_hot.shape[2]).to(device)
            for seq_index in range(inp_smile_hot.shape[1]):
                decoded_one_hot_line, hidden  = model_decode(latent_points, hidden)
#                print('DECODED_ONE_HOT_LINE', decoded_one_hot_line.shape)
#                print('DECODED_ONE_HOT', decoded_one_hot.shape)
                decoded_one_hot[:, seq_index, :] = decoded_one_hot_line[0]

            decoded_one_hot = decoded_one_hot.reshape(batch_size * inp_smile_hot.shape[1], inp_smile_hot.shape[2])
            _, label_atoms  = inp_smile_hot.max(2)
            label_atoms     = label_atoms.reshape(batch_size * inp_smile_hot.shape[1])

            criterion   = torch.nn.CrossEntropyLoss()
            recon_loss += criterion(decoded_one_hot, label_atoms)

            loss += recon_loss + KLD_alpha * kld
#            loss = loss
            if tensorBoard_graphing:
                writer.add_scalar('Batch Loss', loss, epoch*(num_batches_train) + batch_iteration)

            # perform back propogation
            optimizer_encoder.zero_grad()
            optimizer_decoder.zero_grad()
            loss.backward(retain_graph=True)
            nn.utils.clip_grad_norm_(model_decode.parameters(), 0.5)
            optimizer_encoder.step()
            optimizer_decoder.step()

            if batch_iteration % 30 == 0:
                end = time.time()

                # assess reconstruction quality
                _, decoded_max_indices = decoded_one_hot.max(1)
                _, input_max_indices   = inp_smile_hot.reshape(batch_size * inp_smile_hot.shape[1], inp_smile_hot.shape[2]).max(1)

                differences = 1. - torch.abs(decoded_max_indices - input_max_indices)
                differences = torch.clamp(differences, min = 0., max = 1.).double()
                quality     = 100. * torch.mean(differences)
                quality     = quality.detach().cpu().numpy()

                qualityValid=quality_in_validation_set(data_valid)

                new_line = 'Epoch: %d,  Batch: %d / %d,\t(loss: %.4f\t| quality: %.4f | quality_valid: %.4f)\tELAPSED TIME: %.5f' % (epoch, batch_iteration, num_batches_train, loss.item(), quality, qualityValid, end - start)
                print(new_line)
                start = time.time()



        qualityValid = quality_in_validation_set(data_valid)
        quality_valid_list.append(qualityValid)

        # only measure validity of reconstruction improved
        quality_increase = len(quality_valid_list) - np.argmax(quality_valid_list)
        if quality_increase == 1 and quality_valid_list[-1] > 60.:
            corr, unique = latent_space_quality(latent_dimension,sample_num = sample_num)
        else:
            corr, unique = -1., -1.

        new_line = 'Validity: %.5f %% | Diversity: %.5f %% | Reconstruction: %.5f %%' % (corr * 100. / sample_num, unique * 100. / sample_num, qualityValid)

        print(new_line)
        with open('results.dat', 'a') as content:
            content.write(new_line + '\n')

        if quality_valid_list[-1] < 70. and epoch > 200:
            break

        if quality_increase > 20: # increase less than one percent in 3 episodes
            print('Early stopping criteria')
            break

        if quality_increase == 1:
            if checkpoint:
                save_models(model_encode, model_encode, epoch)






if __name__ == '__main__':
    try:
        content = open('logfile.dat', 'w')
        content.close()
        content = open('results.dat', 'w')
        content.close()

        if os.path.exists("settings.yml"):
            user_settings=yaml.load(open("settings.yml", "r"))
            settings = user_settings
        else:
            print("Expected a file settings.yml but didn't find it.")
            print("Create a file with the default settings.")
            print("Please check the file before restarting the code.")
            print()
            exit()

        cuda_device = settings['data']['cuda_device']
        os.environ['CUDA_VISIBLE_DEVICES'] = str(cuda_device)

        type_of_encoding = settings['data']['type_of_encoding']

        if type_of_encoding == 0:
            settings['data']['smiles_file']='VAE_dependencies/Datasets/QM9/0SelectedSMILES_QM9.txt'
            encoding_alphabet=' FONC()=#12345'
        elif type_of_encoding == 1:
            settings['data']['smiles_file']='VAE_dependencies/Datasets/QM9/1DeepSMILES_QM9.txt'
            encoding_alphabet=' FONC)=#3456789'
        elif type_of_encoding == 2:
            settings['data']['smiles_file']='VAE_dependencies/Datasets/QM9/2RGSMILES_QM9.txt'
            encoding_alphabet='ABCDEFGHIJKLMN'
        elif type_of_encoding == 3:
            settings['data']['smiles_file']='VAE_dependencies/Datasets/QM9/3GrammarVAE_QM9.txt'
            encoding_alphabet='ABCDEFGHIJKLMNOPQRSTUVQXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%[]&*()-_='

        data_parameters = settings['data']
        batch_size = data_parameters['batch_size']
        file_name = data_parameters['smiles_file']

        df = pd.read_csv(file_name)
        smiles_list = np.asanyarray(df.smiles)
        largest_smile_len = len(max(smiles_list, key=len))
        print('Acquiring data...')
        data = multiple_smile_to_hot(smiles_list, largest_smile_len, encoding_alphabet, type_of_encoding)
        print('Data Acquired.')


        len_max_molec = data.shape[1]
        len_alphabet = data.shape[2]
        len_max_molec1Hot = len_max_molec * len_alphabet


        encoder_parameter = settings['encoder']
        decoder_parameter = settings['decoder']
        training_parameters = settings['training']

        model_encode = VAE_encode(**encoder_parameter)
        model_decode = VAE_decode(**decoder_parameter)

        model_encode.train()
        model_decode.train()

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print('*'*15, ': -->', device)

        data = torch.tensor(data, dtype=torch.float).to(device)

        train_valid_test_size=[0.5, 0.5, 0.0]
        x = [i for i in range(len(data))]  # random shuffle input
        shuffle(x)
        data = data[x]
        idx_traintest=int(len(data)*train_valid_test_size[0])
        idx_trainvalid=idx_traintest+int(len(data)*train_valid_test_size[1])
        data_train=data[0:idx_traintest]
        data_valid=data[idx_traintest:idx_trainvalid]
        data_test=data[idx_trainvalid:]

        num_batches_train = int(len(data_train) / batch_size)
        num_batches_valid = int(len(data_valid) / batch_size)

        model_encode = VAE_encode(**encoder_parameter).to(device)
        model_decode = VAE_decode(**decoder_parameter).to(device)
        print("start training")
        train_model(data_train=data_train, data_valid=data_valid, **training_parameters)

        with open('COMPLETED', 'w') as content:
            content.write('exit code: 0')

#    except Exception as e:
    except AttributeError:
        _, error_message,_ = sys.exc_info()
        print(error_message)

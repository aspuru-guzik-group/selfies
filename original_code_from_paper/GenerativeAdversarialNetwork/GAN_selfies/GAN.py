"""
@author: Akshat
"""
import torch
import random
import torch.utils.data
from torch import nn, optim
from torch.autograd.variable import Variable
import os
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem import Draw
from one_hot_converter import multiple_smile_to_hot, hot_to_smile, check_conversion_bijection
from GPlus2S import GrammarPlusToSMILES, IncludeRingsForSMILES
import numpy as np
import matplotlib.pyplot as plt
from tensorboardX import SummaryWriter
writer = SummaryWriter()

random.seed(4001)

def load_data(cut_off=None):
    '''
    Ensuring correct Bijection:
        check_conversion_bijection(smiles_list=A, largest_smile_len=len(max(A, key=len)), alphabet=alphabets)         
    '''
    with open('2RGSMILES_QM9.txt') as f:
        content = f.readlines()
    content = content[1:]
    content = [x.strip() for x in content] 
    A = [x.split(',')[1] for x in content]
    if cut_off is not None:
        A = A[0:cut_off]
    return A

alphabets = ['I', 'N', 'D', 'K', 'L', 'B', 'C', 'J', 'G', 'A', 'H', 'E', 'F', ' ']
            
                         
# Read in the QM9 dataset
A = load_data(cut_off=None)


data = multiple_smile_to_hot(A, len(max(A, key=len)), alphabets, 0)

print('Data shape: ', data.shape)
one_hot_len_comb = data.shape[1]*data.shape[2]
data = data.reshape(( data.shape[0], one_hot_len_comb))
data = torch.tensor(data, dtype=torch.float)

data_loader = torch.utils.data.DataLoader(data, batch_size=1024, shuffle=True)
num_batches = len(data_loader)
print('DATA Acquired!')

def get_canon_smi_ls(smiles_ls):
    '''
    return a list of canonical smiles in smiles_ls
    '''    
    canon_ls = [Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=False, canonical=True) for smi in smiles_ls]
    return canon_ls

def _make_dir(directory):
    os.makedirs(directory)
        

def save_models(generator, discriminator, epoch, dir_name):
    out_dir = './{}/saved_models/{}'.format(dir_name, epoch)
    _make_dir(out_dir)
    torch.save(generator,     '{}/G'.format(out_dir))
    torch.save(discriminator, '{}/D'.format(out_dir))
        

def display_status(epoch, num_epochs, n_batch, num_batches, d_error, g_error, d_pred_real, d_pred_fake):
    
    # var_class = torch.autograd.variable.Variable
    if isinstance(d_error, torch.autograd.Variable):
        d_error = d_error.data.cpu().numpy()
    if isinstance(g_error, torch.autograd.Variable):
        g_error = g_error.data.cpu().numpy()
    if isinstance(d_pred_real, torch.autograd.Variable):
        d_pred_real = d_pred_real.data
    if isinstance(d_pred_fake, torch.autograd.Variable):
        d_pred_fake = d_pred_fake.data
    
    
    print('Epoch: [{}/{}], Batch Num: [{}/{}]'.format(
        epoch,num_epochs, n_batch, num_batches)
         )
    print('Discriminator Loss: {:.4f}, Generator Loss: {:.4f}'.format(d_error, g_error))
    print('D(x): {:.4f}, D(G(z)): {:.4f}'.format(d_pred_real.mean(), d_pred_fake.mean()))
    writer.add_scalar('D(x)', d_pred_real.mean(), epoch)
    writer.add_scalar('D(G(z)', d_pred_fake.mean(), epoch)




class DiscriminatorNet(torch.nn.Module):
    """
    A three hidden-layer discriminative neural network
    """
    def __init__(self, drop_rate, layer_2_size):
        super(DiscriminatorNet, self).__init__()
        n_features = one_hot_len_comb
        n_out = 1
        
        self.hidden0 = nn.Sequential( 
            nn.Linear(n_features, layer_2_size),
            nn.LeakyReLU(0.2),
            nn.Dropout(drop_rate)
        )

        self.out = nn.Sequential(
            torch.nn.Linear(layer_2_size, n_out),
            nn.Sigmoid()        
        )

    def forward(self, x):
        x = self.hidden0(x)
        x = self.out(x)
        return x



class GeneratorNet(torch.nn.Module):
    """
    A three hidden-layer generative neural network
    """
    def __init__(self, prior_lv_size, layer_interm_size):
        super(GeneratorNet, self).__init__()
        n_out = one_hot_len_comb  
        
        self.hidden0 = nn.Sequential(
            nn.Linear(prior_lv_size, layer_interm_size),
            nn.LeakyReLU(0.2)
        )
        
        self.out = nn.Sequential(
            nn.Linear(layer_interm_size, n_out),
            nn.Sigmoid() 
        )

    def forward(self, x):
        x = self.hidden0(x)
        x = self.out(x)
        return x
    

def noise(size, G_start_layer_size):
    '''
    Standard nois, which acts as input to the GAN generator .
    '''
    n = Variable(torch.randn(size, G_start_layer_size))
    if torch.cuda.is_available(): return n.cuda() 
    return n


def train_discriminator(optimizer, real_data, fake_data, discriminator, criterion):
    optimizer.zero_grad()
    
    # 1.1 Train on Real Data
    prediction_real = discriminator(real_data)
    y_real = Variable(torch.ones(prediction_real.shape[0], 1))
    if torch.cuda.is_available(): 
        D_real_loss = criterion(prediction_real, y_real.cuda())
    else: 
        D_real_loss = criterion(prediction_real, y_real)

    # 1.2 Train on Fake Data
    prediction_fake = discriminator(fake_data)
    y_fake = Variable(torch.zeros(prediction_fake.shape[0], 1))
    if torch.cuda.is_available(): 
        D_fake_loss = criterion(prediction_fake, y_fake.cuda())
    else: 
        D_fake_loss = criterion(prediction_fake, y_fake)
    
    D_loss = D_real_loss + D_fake_loss
    D_loss.backward()
    optimizer.step()
    
    return D_real_loss + D_fake_loss, prediction_real, prediction_fake, discriminator


def train_generator(optimizer, fake_data, criterion, discriminator):
    optimizer.zero_grad()
    prediction = discriminator(fake_data)
    y = Variable(torch.ones(prediction.shape[0], 1))
    if torch.cuda.is_available(): 
        G_loss = criterion(prediction, y.cuda(0))
    else: 
        G_loss = criterion(prediction, y)
    G_loss.backward()

    optimizer.step()
    return G_loss.data.item(), discriminator


def train_gan(lr_disc, lr_genr, prior_lv_size, layer_interm_size_G, discr_dropout, discr_layer_2_size, num_epochs, dir_name, num_unique):
    '''
    All the hyper parameters are to be added as parameters to this function!!!
    '''    
    criterion = nn.BCELoss() 
    discriminator = DiscriminatorNet(drop_rate=discr_dropout, layer_2_size=discr_layer_2_size)
    generator     = GeneratorNet(prior_lv_size=prior_lv_size, layer_interm_size=layer_interm_size_G)
    criterion = nn.BCELoss() 
    if torch.cuda.is_available():
        discriminator.cuda()
        generator.cuda()
        
    # Optimizers (notice the use of 'discriminator'<-Object class)
    d_optimizer = optim.Adam(discriminator.parameters(), lr=lr_disc)
    g_optimizer = optim.Adam(generator.parameters(), lr=lr_genr) # 1: 3e-4; 2: 1e-5
    
    
    for epoch in range(num_epochs+1): 
        print('Epoch: ', epoch)
        #   batch_num   real_data
        for n_batch, real_batch in enumerate(data_loader):
    
            # 1. Train Discriminator
            real_data = Variable(real_batch)
            
            if torch.cuda.is_available(): 
                real_data = real_data.cuda()
                
            # Generate fake data
            fake_data = generator(noise(real_data.size(0), prior_lv_size)).detach()        
            
            # Train D
            d_error, d_pred_real, d_pred_fake, discriminator = train_discriminator(d_optimizer, real_data, fake_data, discriminator, criterion)
    
    
            # 2. Train Generator
            # Generate fake data
            fake_data = generator(noise(real_batch.size(0), prior_lv_size))
            
            # Train G
            g_error, discriminator = train_generator(g_optimizer, fake_data, criterion, discriminator)
    
                    
            # Log onto a graph
            writer.add_scalar('D_error', d_error, epoch * num_batches + n_batch)
            writer.add_scalar('G_error', g_error, epoch * num_batches + n_batch)
            
        if epoch % 10 == 0 and epoch > 0:
            generator = generator.eval()
            # Display complete training stats
            display_status(
                epoch, num_epochs, n_batch, num_batches,
                d_error, g_error, d_pred_real, d_pred_fake
            )
            
            smiles_ls = []
            print('Sampling....')
            for _ in range(10000):
                
                T = generator(noise(1, prior_lv_size)).cpu().detach().numpy().flatten() # Just chose one molecule
                T = T.reshape((21, 14))
                
                smile = hot_to_smile(T, alphabets)
                
                if ' ' in smile:
                    smile = smile[:smile.find(' ')]
                    
                # Convert SELFIE to SMILE here
                smile = IncludeRingsForSMILES(GrammarPlusToSMILES(smile,'X0'))
                print(smile)
                mol = Chem.MolFromSmiles(smile)
                if mol is not None:
                    smiles_ls.append(smile)
            print('unique molecules: ', len(set(get_canon_smi_ls(smiles_ls))))
            writer.add_scalar('Num Smiles', len(set(get_canon_smi_ls(smiles_ls))), epoch)
            
            # Write TensorBoard curv on to a text file
            f = open('{}/smiles_curve.txt'.format(dir_name), 'a+')
            f.write(str(len(set(get_canon_smi_ls(smiles_ls))))+ '\n')
            f.close()            
            
            # Save the TensorBoard models
            save_models(generator, discriminator, epoch, dir_name)
                
            num_unique.append(len(set(get_canon_smi_ls(smiles_ls))))
            
            
            #  Early stopping criteria for the model to stop training
            if epoch > 100:
                A = np.array(num_unique)
                stopping_criteria = A.max() - ( 6 *   ( np.sqrt(A.max())))
                if stopping_criteria < 0:
                    stopping_criteria = 0
                print('Stopping Criteria: ', stopping_criteria)
                print('Array A: ', A)
                if num_unique[-1] <= stopping_criteria or max(A) < 50:
                    # Write the stopping epoch onto a text file
                    f = open('{}/stoping_epoch.txt'.format(dir_name), 'a+')
                    f.write('Early stopping Epoch: ' +str(epoch)+ '\n')
                    f.close()            
                    return 

            print(set(smiles_ls))
            generator = generator.train()
            
            
            
if __name__ == '__main__':
    
    num_model = 1 # Number of GAN models to be run (each will have different - randomly initialized hyperparms)
    for model_iter in range(100):
        dir_name = str(model_iter)    # Directory for saving all the data
                          
        os.makedirs(dir_name)
        num_unique = []
    
        ## HYPERPARAMETER SELECTION!

        num_epochs = 1500
    
        # Let me select the learning rates:
        lr_disc = 10 ** random.uniform(-7, -4)
        lr_genr = 10 ** random.uniform(-7, -4)
        
        # Generator
        prior_lv_size       = random.randint(50,  300)
        layer_interm_size_G = random.randint(300, 3000)
        
        # Discriminator
        discr_dropout       = random.uniform(0, 0.8)
        discr_layer_2_size  = random.randint(5, 100)

        
        # Save all the selected hyperparamters: 
        f = open('{}/hyperparams.txt'.format(dir_name), 'a+')
        f.write('lr Discr: ' + str(lr_disc) + '\n')
        f.write('lr Gener: ' + str(lr_genr)+ '\n')
        f.write('Gener Sampling layer size: ' + str(prior_lv_size)+ '\n')
        f.write('Gener middle layer size: ' + str(layer_interm_size_G)+ '\n')
        f.write('Discr dropout: ' + str(discr_dropout)+ '\n')
        f.write('Discr middle layer size: ' + str(discr_layer_2_size)+ '\n')
        f.close()
            
        train_gan(lr_disc, lr_genr, prior_lv_size, layer_interm_size_G, discr_dropout, discr_layer_2_size, num_epochs, dir_name, num_unique)
    
    
    
    
    

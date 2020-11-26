# SELFIES Example: Variational Autoencoder (VAE) for Chemistry

An implementation of a variational autoencoder that runs on both SMILES and 
SELFIES. Included is code that compares the SMILES and SELFIES representations
for a VAE using reconstruction quality, diversity, and latent space validity 
as metrics of interest. 
 
## Dependencies 
Dependencies are ``pytorch``, ``rdkit``, and ``pyyaml``, which can be installed 
using Conda. 
      
## Files 

 * ``chemistry_vae.py``: the main file; contains the model definitions, 
    the data processing, and the training.
 * ``settings.yml``: a file containing the hyperparameters of the 
    model and the training. Also configures the VAE to run on either SMILES 
    or SELFIES. 
 * ``data_loader.py``: contains helper methods that convert SMILES and SELFIES
    into integer-encoded or one-hot encoded vectors. 
    
### Tested at:
- Python 3.7 
              
CPU and GPU supported

For comments, bug reports or feature ideas, please send an email to
mario.krenn@utoronto.ca and alan@aspuru.com 
 

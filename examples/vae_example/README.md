## SELFIES Example: Variational Auto Encoder (VAE) for chemistry

by Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich,
Alan Aspuru-Guzik
- https://arxiv.org/abs/1905.13741
           
- comparing SMILES and SELFIES representation using reconstruction quality, 
diversity and latent space validity as metrics of interest

v0.1.0 -- 04. August 2019             
                  
### information:
     ML framework: pytorch
     chemistry framework: RDKit     
          
     chemistry_vae.py
             main file, includes the model definitions, the training, the 
             representation conversion
     
     settings.yml
             contains link to data file containing SMILES encoded molecule, and
             hyperparameters of neural network model and training
             
     chemistry_vae.py: get_selfie_and_smiles_encodings_for_dataset
             generate complete encoding (inclusive alphabet) for SMILES and 
             SELFIES given a data file
                          
     chemistry_vae.py: VAEEncoder
             fully connection, 3 layer neural network - encodes a one-hot 
             representation of molecule (in SMILES or SELFIES representation) 
             to latent space
             
     chemistry_vae.py: VAEDecoder
             decodes point in latent space using an RNN
             
     chemistry_vae.py: latent_space_quality
             samples points from latent space, decodes them into molecules,
             calculates chemical validity (using RDKit's MolFromSmiles), 
             calculates diversity
             
     environment.yml
             shows dependencies
             Particularly important: RDKit and SELFIES 
             (via 'pip install selfies')
             
     original\
             contains the original codes and settings of hyperparameters in 
             the paper. Not so easy to adapt to other datasets, so we 
             recommend starting from chemistryVAE.py
    
    
### tested at:
- Python 3.7.1
- Python 3.6.8               
CPU and GPU supported
     

### Robustness:
- semantic validity is only implemented so far for atoms described in Table 2 
  of our paper. This corresponds to (non-ionic) QM9. Other chemical constraints
  might generate additional mistakes. Syntactical constraints are always
  fulfilled
- aromatic Symbols: they have additional semantic constraints, thus to reduce 
  invalidity due to aromatic constraints, one can de-aromatize molecules 
  (aromatic symbols are simplifications in SMILES). Otherwise, one could add 
  the semantic constraints (this could be done in an automated way, but is 
  not implemented yet)
  

For comments, bug reports or feature ideas, please send an email to
mario.krenn@utoronto.ca and alan@aspuru.com 
 

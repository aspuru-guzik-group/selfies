# Package requirements for running the code: 
- Pytorch
- TensorBoardX
- rdkit 
- numpy 
- matplotlib


# File Navigator: 

GAN_selfies: 
- 2RGSMILES_QM9.txt : Dearomatized QM9 dataset (selfies representation, from alphabet described in the main text (symbols shortened for simplicity))
- GAN.py: Code for running the generative adversarial network
- one_hot_converter.py: Code for creating one-hot-encodings of molecular strings
- adjusted_selfies_fcts.py: SMILES to SELFIES conversion file
- GPlus2S.py: SMILES to SELFIES conversion file
- translate.py: General helper functions file

GAN_smiles: 
- GAN.py: Code for running the generative adversarial network
- one_hot_converter.py: Code for creating one-hot-encodings of molecular strings
- smiles_qm9.txt: Dearomatized QM9 dataset (smiles representation)

# How to Run the code:? 
Step1: cd inside either 'GAN_selfies' or 'GAN_smiles' depending on which type of molecular representation you would like to run the GAN

Step2: run ./GAN.
       The code will automatically detect the availability of a GPU on your device, and run multiple models with different hyperparameters

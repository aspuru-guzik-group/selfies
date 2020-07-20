The main source code for ``selfies``. 
 * ``decoder.py`` contains the ``selfies.decoder`` method, which translates
    from SELFIES to SMILES.
 * ``encoder.py`` contains the ``selfies.encoder`` method, which translates
    from SMILES to SELFIES.  
 * ``grammar_rules.py`` and ``state_library.py`` contain various global
    variables and helper methods that specify and enforce the SELFIES grammar. 
 * ``kekulize.py`` contains the code that kekulizes SMILES with aromatic 
    symbols, before they are translated into SELFIES. 
 * ``utils.py`` contains the helper methods of ``selfies``, such as 
    ``selfies.len_selfies``.


Pathminer: A simple biological pathway mining algorithm
========================================================

This package provides a simple Python API for biological pathway mining using the 
Biocyc (+MetaCyc) web services. It is built on the separate package `biocyc` which 
handles the actual connections.

This algorithm has successfully been used in the Pathomx software to identify gene
and metabolite changes within pathways. 

Parameters to the function as follows:

.. code:: python

    data: a list of tuples in the format (biocyc_id, score); multiple datasets as a list of lists
    target: what to mine - pathways, reactions or compartments (not supported yet)
    algorithm: one of 'c', 'u', 'd', 'm', 't' 
    
        'Compound change scores for pathway': 'c'
        'Compound up-regulation scores for pathway': 'u'
        'Compound down-regulation scores for pathway': 'd'
        'Number compounds with data per pathway': 'm'
        'Pathway overall tendency': 't'
    
    no_of_results: number of targets to return
    shared: 'share' scores between pathways, to minimise influence of promiscuous metabolites
    relative: 'scale' scores to the size of pathways, to reduce influence of pathway size


target=TARGET_PATHWAYS, algorithm='c', no_of_results=5, shared=True, relative=False


Simple usage is as follows.

.. code:: python

    from biocyc import biocyc
    biocyc.set_organism('meta')
    
    from pathminer import mining
    
    a = biocyc.get('FRUCTOSE-16-DIPHOSPHATE')
    b = biocyc.get('GLC')
    
    data = [(a,5),(b,10)]
    
    mining(data)
    
    [(glycolysis, 2.5), (gluconeogenesis, 2.5), (glycolysis, 2.5), (gluconeogenesis, 2.5), (glycogenolysis, 0.625)]
    
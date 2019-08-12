'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy as np
import RNA as rna 
import random 
import matplotlib.pyplot as plt
from rnaevol.ea import Parameter, inverse, run
from rnaevol.landscape import Landscape
import pandas


#Main function 
def main() : 
    
    target_structure = "((....)).((....)).((....)).((....))"
    population_size = 100
    number_of_generation = 100
    mut_rate = 0.4
    mut_bp = 0.5
    mut_prob = 1./len(target_structure)
    njobs = 1

    type_ = "ham"
    method = "MED"
    success = [] 
    nucleotides = ['A', 'C', 'G', 'U']
    base_paires = ["AU", "UA", "GC", "CG", "GU", "UG"]

    p_n = [.25,.25,.25,.25]
    p_c = [0.1, 0.1, 0.2,0.2,0.2,0.2]
    
    mut_probs = np.array(rna.ptable(target_structure)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0

    landscape = Landscape(type_, target_structure)
    #print(landscape.__dict__)
    #params = Parameter(nucleotides, base_paires, method,population_size,landscape,mu=mut_probs,mu_bp=mut_bp,generation=number_of_generation,p_n=p_n, p_c=p_c)
    #print(params.__dict__)
    
    params = []
    for i in range(njobs) : 
        param = Parameter(nucleotides, base_paires, method,population_size,landscape,mu=mut_probs,mu_bp=mut_bp,generation=number_of_generation,p_n=p_n, p_c=p_c, log_folder="rnaevol_log/"+str(i))
        params.append(param)
    
    df = run(params, pe=True, bt_log=True)
    print (df)
if __name__== "__main__" : 

    main()

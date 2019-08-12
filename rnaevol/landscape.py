'''
@author: Nono Saha Cyrille Merleau 
@email: csaha@aims.edu.gh/ nonosaha@mis.mpg.de

This program is an implementation of a landscape function. 

The landscape function is given by: 


'''


#Importing the libraries 
import numpy as np 
import RNA
import subprocess 
from dataclasses import dataclass

@dataclass(init=True)
class Landscape(object) : 

        TREE_TYPE = "tree"
        BP_TYPE = "bp" 
        HAM_TYPE = "ham"

        type_ : str 
        target_structure : str
        min_energy : float = 5.0


        
        #Compute the fitness of an RNA Structure
        def eval_fitness(self, rna_structure:str) -> float :  
                if self.type_ == self.BP_TYPE : 
                        return 1./(1.+RNA.bp_distance(self.target_structure,rna_structure))
                elif self.type_ == self.TREE_TYPE : 
                        return tree_edit_fitness(rna_structure)
                else: 
                        return 1/(1. + RNA.hamming_distance(self.target_structure, rna_structure))

        
   

#Compute the tree edit distance between two RNA structures
def tree_edit_fitness(rna_structure1:str, rna_structure2:str) -> float : 
        ref_xstrc = RNA.expand_Full(rna_structure1)
        xstrc = RNA.expand_Full(rna_structure2)
        
        return 1./(1.+RNA.tree_edit_distance(RNA.make_tree(ref_xstrc), RNA.make_tree(xstrc)))


def min_ens_distance(target_structure:str, rna_sequence:str, min_energy=5.0) -> float : 

        rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        rnasubopt_out, rnasubopt_err = rnasubopt.communicate(rna_sequence.encode()) 

        cut_pipe = subprocess.Popen(['cut', '-f', '1', '-d', ' '], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cut_out, stderr = cut_pipe.communicate(rnasubopt_out)
        
        list_subopt =  cut_out.split()[1:] 

        return max([1./(1+RNA.hamming_distance(target_structure,str(subopt))) for subopt in list_subopt])

def ens_defect(target_structure:str, rna_sequence:str) -> float: 
        p = subprocess.Popen("defect", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_input = rna_sequence+"\n"+target_structure+"\n"
        
        defect, error = p.communicate(cmd_input.encode())
        defect = defect.split("\n".encode())
        return 1./float(defect[-3])

def ens_diversity(rna_sequence:str, min_energy = 3.5) -> float: 
        rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        rnasubopt_out, rnasubopt_err = rnasubopt.communicate(rna_sequence) 
        result = np.array([[s.split()[0], s.split()[1]] for s in rnasubopt_out.split("\n")[1:-1]])
        
        kt = 0.612 
        z = sum(np.exp(-np.array(result[:,1], dtype=float)/kt))
        p = np.exp(-np.array(result[:,1], dtype=float)/kt)/z

        return 1-sum(p**2) 
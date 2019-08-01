
'''
    @author: Nono Saha Cyrille Merleau 
    @email: nonosaha@mis.mpg.de 

    This class encode the structure of an RNA individual.
'''

import RNA
from dataclasses import dataclass

@dataclass(init= True)
class Individual(object) : 

    rna_sequence : str 
    rna_structure : str 
    mfe : float 
    fitness : float 
    

    def as_dict(self):
        return self.__dict__

def get_bp_position(structure: str) ->  list: 
    position = RNA.ptable(structure)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
        if position[i] != 0 :
            if (position[i]-1,i) in base_paire_pos : 
                continue;  
            else :
                base_paire_pos.append((i,position[i]-1))

    return base_paire_pos
    
    
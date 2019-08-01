'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''


#Native Lib
import random 
import os

#External Lib
from dataclasses import dataclass
from typing import List
import numpy 
import RNA 
import pandas

#Internal module
from rnaevol.individual import *
from rnaevol import landscape 


@dataclass(init=True)
class Initializer(object) : 


    DEFAULT_CONONICAL_BP = list(["GC","CG","GU","UG","AU","UA"])
    DEFAULT_NUCLEOTIDES = list(["A", "C", "G","U"])

    
    population_size: int
    landscape: landscape.Landscape


    def init_with_const(self,nucleotides=DEFAULT_NUCLEOTIDES,base_paires=DEFAULT_CONONICAL_BP )->List[Individual] : 

        pos = get_bp_position(self.landscape.target_structure)
        pop = []
        i = 0
        while i < self.population_size : 
            if i < 4 : 
                print(self.landscape.target_structure,nucleotides)
                arn = numpy.random.choice(nucleotides[i:i+1],len(self.landscape.target_structure))
            else : 
                arn = numpy.random.choice(nucleotides,len(self.landscape.target_structure))

            for bp_cord in pos : 
                bp = numpy.random.choice(base_paires,1)
                arn[bp_cord[0]] = bp[0][0]
                arn[bp_cord[1]] = bp[0][1]
            pop.append(''.join(arn))
            i = len(set(pop))

        pop = numpy.array(list(set(pop)))

        init_pop = []
        for seq in pop:
            strc,mfe = RNA.fold(seq)
            init_pop.append(Individual(seq, strc, mfe, self.landscape.eval_fitness(strc)))

        return init_pop

        

    def simple_init(self,nucluotides=DEFAULT_NUCLEOTIDES ) -> List[Individual] : 

        population = []
        init_depth = len(self.landscape.target_structure)
        for i in range(self.population_size):
            if i < 4 : 
                arn = numpy.random.choice(nucluotides[i:i+1],init_depth)
                seq = ''.join(arn)
            else :
                arn = numpy.random.choice(nucluotides,init_depth)
                seq = ''.join(arn)
            (strc, mfe) = RNA.fold(seq)

            ind = Individual(seq,strc, mfe,self.landscape.eval_fitness(strc))
            population.append(ind)

        return population



    def initialize_from_csv(self, path:str, sep:str) -> List[Individual] : 
        
        try:
            dataFrame = pandas.read_csv(path,sep=sep)

        except pandas.errors.EmptyDataError:
            print("Error occured while opening a file ", path)
            return
        
        ref_ind = Individual((dataFrame.values)[0,1:][0],(dataFrame.values)[0,1:][1],(dataFrame.values)[0,1:][3],(dataFrame.values)[0,1:][2])
        
        init_pop = []
        for ind in (dataFrame.values)[1:self.population_size+1,1:] : 
            init_pop.append(Individual(ind[0],ind[1],ind[3],ind[2]))

        return init_pop



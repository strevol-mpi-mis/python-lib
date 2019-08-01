'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA 
import pandas
import os
import subprocess

from dataclasses import dataclass
from typing import List

from rnaevol.individual import *
from rnaevol.initializer import Initializer
from rnaevol.logger import Logger
from rnaevol import landscape


@dataclass(init = True)
class Parameter(object) : 

    nucleotides: List[str]
    base_paires: List[str]
    select_meth: str
    population_size: int
    landscape: landscape.Landscape
    mu: List[float] 
    mu_bp:float 
    p_n:List[float]
    p_c:List[float]
    generation: int
    

#Mutation function 
def mutateOne(ind:Individual, params:Parameter) -> Individual :  
    
    assert 0<params.mu_bp<1, "mu_bp should be between 0 and 1"
    assert len(params.p_n) == 4, "The length of p_n should be 4"
    assert len(params.p_c) == 6, "The length of p_n should be 6"

    RNA_seq = numpy.array(list(ind.rna_sequence))
    r = numpy.random.rand(len(ind.rna_sequence))
    mut_pos =RNA_seq[r<params.mu] 
    choices = numpy.random.choice(params.nucleotides, len(mut_pos), p=params.p_n)
    RNA_seq[r<params.mu] = choices 
    pos = get_bp_position(params.landscape.target_structure)

    for bp_cord in pos : 
        r = numpy.random.uniform(0,1)
        if r < params.mu_bp : 
            bp = numpy.random.choice(params.base_paires,1, p=params.p_c)
            RNA_seq[bp_cord[0]] = bp[0][0]
            RNA_seq[bp_cord[1]] = bp[0][1] 
    

    (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
    return Individual(''.join(RNA_seq), RNA_strc, mef,params.landscape.eval_fitness(RNA_strc))



def mutateAll(population: List[Individual], params:Parameter) -> List[Individual]: 
    mutated_pop = [mutateOne(ind,params) for ind in population] 
    return mutated_pop


# Natural selection based on fitness proportionate method
def fitness_proportion_selection(population : List[Individual], size:int) -> List[Individual] : 
    fitnesses = numpy.array([ind.fitness for ind in population])
    selected = numpy.random.choice(population,size=size,p=fitnesses/numpy.sum(fitnesses))
    return selected

def ensDiversity_proportion_selection(population: List[Individual], size:int, params: Parameter)-> List[Individual]: 
    
    ensDiv = numpy.array([landscape.ens_diversity(params.landscape.target_structure,ind.rna_sequence) for ind in population])
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDiv)/sum(ensDiv))
    return selected


def ensDefect_proportion_selection(population: List[Individual], size:int, params: Parameter) : 
    ensDefect = [landscape.ens_defect(params.landscape.target_structure,ind.rna_sequence) for ind in population]
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDefect)/sum(ensDefect))
    return selected

def min_ens_distance_proportion_selection(population: List[Individual], size:int, params: Parameter)->List[Individual] : 
    ensDist = [landscape.min_ens_distance(params.landscape.target_structure,ind.rna_sequence, 5.0) for ind in population]
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDist)/sum(ensDist))
    return selected


def reproduce(population: List[Individual], size:int) : 
    list_fitness = []

    for ind in population : 
        list_fitness.append(ind.fitness)
    list_fitness = sorted(list_fitness, reverse=True) 

    sorted_pop = [ ] 
    for fitness in list_fitness : 
        for ind in population : 
            if ind.fitness == fitness : 
                sorted_pop.append(ind)
            
    return sorted_pop[:size]


##################
'''
    This method is used to cross two given genotypes to form two new offspring for the next generation.
    
    For two given sequences of bit, the crossover consists of selecting n-bit in genotype1 and n-bit in genotype2 
    and exchange them to form two new genotype. 

        INPUT
        =====
                genotype1(Type: string) : a sequence of bit of lenth N. that represent a parent
                genotype2(Type: string) : a sequence of bit of lenth N. that represent a donor
                number_of_bit(Type: int): the default value is 1. It is the number bit to cross(to exchange)
        OUTPUT
        ======
                The function returns two new offsprings.
'''
def crossover(self, ind1, ind2) : 
    vect1 = map(str, ind1.RNA_seq)
    vect2 = map(str, ind2.RNA_seq)
    
    r = numpy.random.randint(0, len(vect1))
    swap = vect1[: r]
    vect1[:r] = vect2[:r]
    vect2[:r] = swap 
    pos = ind1.get_bp_position(self.landscape.target_structure)
    for i in range(r) : 
        for bp_cord in pos:  
            if i in bp_cord :
                pos.remove(bp_cord)
                if vect1[i] == "G" : 
                    n = numpy.random.choice(["U","C"],1)
                    vect1[bp_cord[bp_cord!=i]] = n[0]
                if vect1[i] == "U" : 
                    n = numpy.random.choice(["G","A"],1)
                    vect1[bp_cord[bp_cord!=i]] = n[0]
                if vect1[i] == "A" : 
                    vect1[bp_cord[bp_cord!=i]] = "U"
                if vect1[i] == "C" : 
                    vect1[bp_cord[bp_cord!=i]] = "G"
                    
                if vect2[i] == "G" : 
                    n = numpy.random.choice(["U","C"],1)
                    vect2[bp_cord[bp_cord!=i]] = n[0]
                if vect2[i] == "U" : 
                    n = numpy.random.choice(["G","A"],1)
                    vect2[bp_cord[bp_cord!=i]] = n[0]
                if vect2[i] == "A" : 
                    vect2[bp_cord[bp_cord!=i]] = "U"
                if vect2[i] == "C" : 
                    vect2[bp_cord[bp_cord!=i]] = "G"
            else : 
                n = numpy.random.choice(["A","G","U","C"],1)
                vect1[i] = n[0]
                n = numpy.random.choice(["A","G","C","U"],1)
                vect2[i] = n[0]
                
    (child1_structure, mef_child1) = RNA.fold(''.join(vect1))
    (child2_structure, mef_child2) = RNA.fold(''.join(vect2))
    child1 = Individual(''.join(vect1), child1_structure, mef_child1,self.landscape.fitness(child1_structure))
    child2 = Individual(''.join(vect2), child2_structure, mef_child2,self.landscape.fitness(child2_structure))
    
    return  child1, child2



'''
This function is implementing the simple genetic algorithm
        INPUT
        =====
                population_size(Type: int) : the number of RNA sequences to generate.
                number_of_generation(Type: int) : the number of generation.
                mut_prob(Type: float between [0,1]): the mutation probability.
        OUTPUT
        ======
                The function returns a list of individual of the new generation.
'''

def ea_without_crossover(init_pop: List[Individual], params:Parameter ) -> List[Individual]: 

    print (" Start of evolution ")
    
    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    population_size = len(init_pop)
    n = params.generation
    
    #logger = Logger(str(log_folder),str(self.select_meth))
    #logger.save_population(prev_population,0)
    max_fitness = max([ind.fitness for ind in prev_population])
    newgeneration = numpy.copy(prev_population)
    while (n > 0) and (max_fitness < 1):
    
        if (params.generation - n)%10 == 0 : 
            print ('Generation '+str(params.generation - n), "Max fitness = ", max_fitness)
        newgeneration = []
        newgeneration =  reproduce(prev_population,int(0.1*population_size))
        
        if params.select_meth == "F" : 
            selected_ind =fitness_proportion_selection(prev_population,population_size,params)
        elif  params.select_meth == "ED":
            selected_ind = ensDefect_proportion_selection(prev_population,population_size,params)
        elif params.select_meth == "MED"  : 
            selected_ind = min_ens_distance_proportion_selection(prev_population, population_size,params)
        elif params.select_meth == "EDV"  : 
            selected_ind = ensDiversity_proportion_selection(prev_population, population_size,params)
        
        mutated = mutateAll(selected_ind,params)

        
        newgeneration = numpy.insert(newgeneration, len(newgeneration),mutated)
        prev_population = numpy.copy(newgeneration)
        #logger.save_population(prev_population,number_of_generation-n)
        
        new_max = max([ind.fitness for ind in prev_population])
        if new_max > max_fitness : 
            max_fitness = new_max
        n -=1

    return params.generation - n, prev_population



    
'''
This function is implementing the simple genetic algorithm
        INPUT
        =====
                population_size(Type: int) : the number of RNA sequences to generate.
                number_of_generation(Type: int) : the number of generation.
                mut_prob(Type: float between [0,1]): the mutation probability.
        OUTPUT
        ======
                The function returns a list of individual of the new generation.
'''
"""
def ea_with_crossover(self,number_of_generation, mut_probs, log_folder, mut_bp) : 

    print (" Starting of evolution ")
    prev_population = self.initializer.init() #Initialize the population of RNA
    population_size = self.initializer.population_size
    
    logger = Logger(str(log_folder),str(self.select_meth))
    n = number_of_generation
    maxfitness = max([ind.fitness for ind in prev_population])
    newgeneration = numpy.copy(prev_population)
    logger.save_population(prev_population,0)
    while (n > 0) and (maxfitness<1)   :
        
        if (number_of_generation - n)%20 == 0 : 
            print ('Generation '+str(number_of_generation - n)), "Max fitness = ", maxfitness
        newgeneration =  self.reproduce(prev_population,int(0.1*population_size))
        
        cross_ind_num = int(population_size*0.9)
        
        while cross_ind_num > 0 :
            if self.select_meth == "F" : 
                selected_ind = self.fitness_proportion_selection(prev_population,2)
            elif  self.select_meth == "ED":
                selected_ind = self.ensDefect_proportion_selection(prev_population,2)
            elif self.select_meth == "MED"  : 
                selected_ind = self.min_ens_distance_proportion_selection(prev_population, 2)
            elif self.select_meth == "EDV"  : 
                selected_ind = self.ensDiversity_proportion_selection(prev_population, 2)
            
            new_ind1, new_ind2 = self.crossover(selected_ind[0],selected_ind[1])
            newgeneration= numpy.insert(newgeneration, len(newgeneration),[new_ind1, new_ind2])
            cross_ind_num = cross_ind_num - 2
        
        prev_population = numpy.copy(newgeneration)
        new_max = max([ind.fitness for ind in prev_population])

        if new_max > maxfitness : 
            maxfitness = new_max
        n -=1
        
        logger.save_population(prev_population,number_of_generation-n)
    return prev_population

"""
def inverse(params:Parameter) -> List[Individual] : 
    initializer = Initializer(params.population_size,params.landscape) 
    init_pop = initializer.init_with_const(params.nucleotides, params.base_paires)

    max_gen,best_population = ea_without_crossover(init_pop, params)
    if max_gen < params.generation : 
        print("Solution found after ",max_gen, "generations : ")
        best_solutions = []
        for ind in best_population : 
            if ind.fitness == 1.0 : 
                best_solutions.append([ind.rna_sequence, ind.fitness, ind.mfe,1/landscape.ens_defect(params.landscape.target_structure,ind.rna_sequence)])
        
        df = pandas.DataFrame(best_solutions, columns=["Sequence", "Fitness","MFE","ED"])
        print (df)
        return df
    else: 
        print("Solution not found after ",max_gen, "generations  ")
        return None

def run(self, number_of_generation,mut_probs, mut_bp, nbjobs, log_fold) : 
    

    return None



def run_with_crossover(self, number_of_generation,mut_probs, mut_bp, nbjobs) : 
    
    return None




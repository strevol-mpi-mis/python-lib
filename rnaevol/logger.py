
'''
    @author: Nono Saha Cyrille Merleau  
    @email: nonosaha@mis.mpg.de

'''
import os
import pandas
import time



class Logger(object) : 

    def __init__(self, logfolder, sublogfolder): 
        loc_time = time.localtime()
        date_time= str(loc_time.tm_year)+str(loc_time.tm_mon)+str(loc_time.tm_mday)+str(loc_time.tm_hour)+str(loc_time.tm_min)+str(loc_time.tm_sec)
        self.logfolder = str(os.getcwd())+'/'+logfolder+"/"+date_time
        self.root_path = self.logfolder+"/"+sublogfolder
        try:
            if not os.path.exists(self.root_path):
                os.makedirs(self.root_path)
            
            with open(self.root_path+'/log.txt', 'w'):
                pass
        except OSError :
            print (" Can not initialize the log folder ")

    def save_population(self, population,gen) : 
        data = []
        for ind in population : 
            data.append([ind.rna_sequence, ind.rna_structure,ind.mfe, ind.fitness])
    
        dataFrame = pandas.DataFrame(data)
        dataFrame.to_csv(self.root_path+"/gen"+str(gen)+".csv")

    def bt_save_population(self, prev_pop, population,gen) : 
        data = []
        prev_data = []
        for i in range(len(population)) : 
            data.append([population[i].rna_sequence, population[i].rna_structure,population[i].mfe, population[i].fitness])
            prev_data.append([prev_pop[i].rna_sequence, prev_pop[i].rna_structure,prev_pop[i].mfe, prev_pop[i].fitness])

        
        dataFrame = pandas.DataFrame(data)
        prev_dataFrame = pandas.DataFrame(prev_data)
        prev_dataFrame.to_csv(self.root_path+"/prev_gen"+str(gen)+".csv")
        dataFrame.to_csv(self.root_path+"/gen"+str(gen)+".csv")

    
    def info(self, message:str) : 
        print("[INFO] [" +time.ctime()+"]: "+message)
        log_file = open(self.root_path+'/log.txt', 'a')
        log_file.write("[INFO] [" +time.ctime()+"]: "+message+"\n")
        log_file.close()
        return

    def error(self, message:str) : 
        print("[ERROR]: "+message)
        log_file = open(self.root_path+'/log.txt', 'a')
        log_file.write("[ERROR]: "+message+"\n")
        log_file.close()
        return

    def debug(self, message:str) : 
        print("[DEBUG]: "+message)
        log_file = open(self.root_path+'/log.txt', 'a')
        log_file.write("[DEBUG]: "+message+"\n")
        log_file.close()
        return
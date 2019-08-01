
'''
    @author: Nono Saha Cyrille Merleau 
    @email: nonosaha@mis.mpg.de 

    This class tests if the instantiation of an individual. 
'''

from rnaevol import individual

def test_getposition() : 

    ind = individual.Individual("CGGAAACG", "((....))",1., -1.2)
    pos = individual.get_bp_position("..((.....))..")

    return ind

def main() : 

    print(test_getposition().as_dict())

if __name__ == '__main__' : 

    main()
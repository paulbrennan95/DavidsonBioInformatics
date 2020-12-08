
class Gene:
    def __init__(self, dna_seq):
        ''' Class constructor 
        Args -
        dna_seq - string of sequence of dna bases for gene
        
        returns: self.seq and self.rna (transcribed rna sequence of self.seq)
        '''
        
        self.seq = dna_seq.upper()
        rna = ''
        for base in self.seq:
            if base == 'C':
                rna += 'G'
            elif base == 'G':
                rna += 'C'
            elif base == 'T':
                rna += 'A'
            elif base == 'A':
                rna += 'U'                
        self.rna = rna
        
    def gc(self):
        '''calculates GC content of dna sequence
        
           Returns - floating point value of GC content
        '''
        return (self.seq.count('G') + self.seq.count('C'))/float(len(self.seq))
    
    
    def find_word(self, word):
        ''' 
        The function finds the quanitity and locations of instances of string 
        word in the dna seq
        
        Args:
         word - sequence to search for in dna seq 
                 
        Returns: a tuple of (word quantity, list of locations)
        '''
        
          
        #initialize list of locations and number of instances of word
        word = word.upper() 
        list_loc = []
        num_words = 0
        
        #determine the locations and number of instances of word in dna seq
        for i in range(len(self.seq)):
            if self.seq[i: i + len(word)] == word:
                num_words = num_words + 1
                list_loc.append(i)
                  
        return num_words, list_loc
          
    def __str__(self):
        '''prints self sequence when called
        '''
        return self.seq
    
    def get_rna(self):
        '''prints self rna sequence determined in constructor
        '''
        return self.rna
    
    def __add__(self, other):
        '''enables class Gene objects to be concatenated to make new class Gene
            object
            
            Arg:
            other - string dna sequence to concatenate to original string dna
            
            Returns: new class Gene object 
        '''
        new_gene = self.seq + other.seq
        return Gene(new_gene)
   

def main():   
    # Construct a new gene from its DNA sequence, and print it
    g = Gene("acatgaatacaatcgaatcg")
    print g
    
    # get_rna is a "getter" method that returns an already existing 
    #      rna attribute of the object
    print g.get_rna()
    
    
    # Compute the G/C proportion in the DNA sequence
    print g.gc()
    
    # Compute the number of occurrences, and each location, of the input word 
    print g.find_word('gaa')
    
    # Construct a second gene and then make a third gene by adding the first two
    #    The DNA of the sum is the concatenated DNA of the first and second genes
    g2 = Gene("actgacatgcatgactca")
    newg = g + g2
    print newg.get_rna()
    
    # Confirm that the first gene has not been modified by any of the above methods
    print g


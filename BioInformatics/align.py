'''
Authors: Paul Brennan and Susannah Cooley
10-24-2017
Time:  1 hr 50 min
'''
import random

def match_score(base1, base2):
    """ Computes the match/mismatch score for a single base comparison
        Args:
        base1 - Base in first sequence
        base2 - Base in second sequence
        
        Returns: 1 if bases are equal, and -3 if bases are not equal
    """
    # if bases are same, score is 1
    if base1 == base2:
        score = 1 
        
    # if bases are not same, score is  -3
    else: 
        score = -3
    return score

def align_matrix(seqA,seqB):
    """ Fills a local sequence alignment matrix 
    
        Uses the Smith-Waterman algorithm to compute the local alignment  
        scoring matrix of two input sequences. The maximum entry in this 
        matrix is the local alignment score. Assumes the cost of a gap is -4.
        Calls the function match_score to find the match/mismatch score 
        of aligning two bases.
    
        Args: 
        seqA - first DNA sequence to be aligned
        seqB - second DNA sequence to be aligned
        
        Returns: a list of lists representing the alignment score matrix        
    """
    #function modified from Paul's HW7
    
    #initialize matrix, and add non-base character to beginnings of sequences
    mat = []
    seqA = '-' + seqA
    seqB = '+' + seqB
    
    #look at base in seqA first
    for (i, base1) in enumerate(seqA):
        #Make all values in first row equal zero
        if i == 0:
            row1 = []
            for pos in range(len(seqB)):
                row1.append(0)
            mat.append(row1)
        #on all other rows, compare base in seqA to all bases in seqB
        else:
            row = []
            #building one row per matrix per one value of seqA
            for (j, base2) in enumerate(seqB):
                #make the first column in matrix all zeros
                if j == 0:
                    row.append(0)
                #for a mismatch, subract 3 from previous pair, or subtract 4 
                # from above or from left in matrix.
                #for a match, add the maximum value achieved of (adding 1 from  
                #  previous pair, subract 4 from above or from left in matrix 
                #    (gap), or value of 0)              
                else:
                    score = match_score(base1, base2)
                    row.append(max((mat[i-1][j-1]+score),(mat[i-1][j]-4),\
                                   (row[j-1]-4), 0))
            
            #add row to matrix     
            mat.append(row)
        
    return mat

def maxmat(A):
    """ Determines the maximum value in a matrix 
        Arg:
        A - nested list representing local alignment matrix
        
        Returns: maximum integer in matrix representing local alignment score
    """
    #initialize list of maximum integers in each row of matrix
    maxes = []
    for row in A:
        #add maximum integer in each row to new list called maxes
        maxes.append(max(row))
    #return the maximum integer in the list of maxes
    return max(maxes)


def randseq(n, p_a):
    """ Generates a random DNA sequence of length n.
    
        Calls randbase to generate a single randomly chosen DNA base using
        Chargoff's rules.
        
        Args:
        n - length of the desired sequence
        p_a - the probability of A in the desired sequence
        
        Returns: string of randomly generated dna sequence of length n
    """
    dna_seq = '' 
    for i in range(n):
        dna_seq += randbase(p_a)
    return dna_seq
   
def randbase(p_a):
    ''' Generates random base according to Chargaff's rules using the
        probability of a 
        
        Arg:
        p_a  - probability of a in desired sequence
        
        Returns: string of single base (a, t, g, c)
        
    '''
    #generates random number between 0 and 1
    r = random.random()
    # assigns a, t, g, or c according to value of r
    if r < p_a:
        base = 'a'
    elif r < p_a*2:
        base = 't'
    elif r < 0.5 + p_a:
        base = 'c'
    else:
        base = 'g' 
    return base    
    
def mean(numbers):
    ''' Finds average of a list of numbers
        
        Arg:
        numbers - list of numbers
        
        Returns: average value of list of numbers as floating point value 
    
    '''
    #initialize summation of numbers in list
    sum_num = 0
    #iterate through list to sum all numbers in list
    for num in numbers:
        sum_num += num
    
    #determine average of list
    m = sum_num/float(len(numbers))
    
    return m

    
def stdev(numbers, m):
    """ Computes the standard deviation of a list of numbers.
    
        Args:
        numbers - a list of numbers to find the standard deviation for
        m - the mean of the list of numbers
        
        Returns: the standard deviation as a floating point value
    """
    #initialize summation of list of numbers
    sum_dev = 0
    #iterate through list of numbers to sum differences between number and mean
    for num in numbers:
        sum_dev += abs(m - num)
    
    # find average deviation from the mean
    std_dev = sum_dev/float(len(numbers))
    
    return std_dev


def alignsim(numseq, lenseq, pa):
    """ Performs a simulation of local alignment for a set of random sequences 
        and runs statistics on local alignment scores. 
        
        Args-
        numseq - number of random sequences to generate
        lenseq - length of each random sequence to generate
        pa - probability of a in each random sequence
        
        Returns: average and standard deviation of all alignment scores between
        each sequence comparison as floating point values
    """
    #initialize lists of sequences and alignment scores
    sequences = []
    align_scores = []

    #creates list of all randomly generated sequences
    for i in range(numseq):
        sequences.append(randseq(lenseq, pa))
    
    #performs local sequence alignment between each sequence combination
    # in list sequences, and creates list of all local alignment scores
    for i in range(len(sequences)):
        for j in range(i+1, len(sequences)):
            mat = align_matrix(sequences[i], sequences[j])
            align_scores.append(maxmat(mat))
    
    #determines mean and standard deviation of list of alignment scores
    align_mean = mean(align_scores)
    align_std_dev = stdev(align_scores, align_mean)
    
    return align_mean, align_std_dev
                     

print align_matrix("gactcctgagggctggcttacgtccgccgtgcagag", "acgtg")
print align_matrix("gactcctgagggctggcttacgtccgccgtgcagag", "gcgtg")

if __name__ == "__main__":
    print alignsim(10, 100, 0.2)
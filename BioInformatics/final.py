'''
Bioinformatics Fall 2017 Take-Home Final
Author: Paul Brennan
Time spent: 8 hours

'''
GAP = -4

def pwn_scan(pwmdict, dna):
    ''' Determines presence of motifs in dna sequences given a certain motifs 
        defined by DNA position weight matrix (PWM) dictionary.
        
    Args:
    pwmdict - dictionary with bases as keys containing the PWM value for each 
               base's value at each respective position in the motif
               
    dna - any sequence in which to scan for motif scores
    
    Returns:  list of motif scores of the first n nucleotides (where n is the 
               length of the motif), then nucleotide 2 to n+1, then nucleotide
               3 to n+2, and so on until the end of the dna sequence. 
               
    '''
    #make dna lowercase
    dna = dna.lower()
    
    #find motif length, regardless of capitalization of dict key
    if 'A' in pwmdict:
        motif_len = len(pwmdict['A'])
    elif 'a' in pwmdict:
        motif_len = len(pwmdict['a'])
    
    #initialize list
    pwm_score_list = []
    
    #determine scoring windows depending on motif length and dna seq length
    for x in range(len(dna)-(motif_len-1)):
        #look at one scoring window
        dna_seg = dna[x:x+motif_len]
        #initialize motif score of scoring window
        pwm_score = 0
        
        #iterate through scoring window, scoring bases in their respective 
        #  positions.
        for i in range(len(dna_seg)):
            if dna_seg[i] == 'a':
                pwm_score += pwmdict['A'][i]
            elif dna_seg[i] == 'c':
                pwm_score += pwmdict['C'][i] 
            elif dna_seg[i] == 'g':
                pwm_score += pwmdict['G'][i]                  
            elif dna_seg[i] == 't':
                pwm_score += pwmdict['T'][i]   
        # add score (rounded to two decimals) to output list       
        pwm_score_list.append(round(pwm_score, 2))
            

    return pwm_score_list


def hairpin(seq, n):
    ''' 
        Decides whether the input sequence can form a hairpin, in which the 
        first base pairs with the last, the second pairs with the next to last,
        and so on, leaving no more than n bases in the middle that are 
        (possibly) unable to pair.
            
        Args:
            seq - input DNA sequence 
            n - number of bases ignored in the middle
            
        Returns: a boolean
            True if sequence satisfies the pairing properties
            False otherwise
    '''
    #Base Case: if the seq is only ignorable bases, its a hairpin!
    if len(seq) == n:
        return True
    
    #if the seq has more than ignorable bases:
    if len(seq) > n:
        #if the first and last bases are complementary, check if interior 
        #   of those bases are complementary, and if so all the way until base
        #   case, the sequence is a hairpin
        if comp(seq[0], seq[-1]):
            return hairpin(seq[1:-1], n)
        else:
            return False
            

def comp(a, b):
    ''' Decides if bases a and b are complementary
    
        Args: 
            a - a DNA base (upper or lower case)
            b - another DNA base (upper or lower case)
        
        Returns: a boolean
            True if bases a and b are complementary
            False otherwise
    '''
    a = a.lower()
    b = b.lower()
    
    #make complementary bases dictionary
    comp = {'g':'c','c':'g','a':'t','t':'a'}
    
    # x is not determinable 
    if a == 'x' or b == 'x':
        return False

    #if a and b are complementary, the complement of a should equal b
    if comp[a] == b:
        return True
    else:
        return False


def match_score(base1, base2):
    base1 = base1.lower()
    base2 = base2.lower()
    if base1 == base2:
        match_score = 1
    else:
        match_score = -3
    return match_score

# Fills matrix with scores using Smith-Waterman algorithm.
def align2seq(seqA,seqB):
    """ Fills a global sequence alignment matrix (Smith-Waterman alignment)
    
        Args: 
        seqA - first DNA sequence to be aligned
        seqB - second DNA sequence to be aligned
        
        Returns: a list of lists representing the alignment score matrix        
    """
    rows = len(seqA) + 1
    cols = len(seqB) + 1
    A = [range(0, (len(seqB)+1)*GAP, GAP)]
    for i in range(1,rows):
        A.append( [i*GAP] + [0]*len(seqB) )
        for j in range(1,cols):
            A[i][j] = max( 
                          A[i-1][j-1] + match_score(seqA[i-1],seqB[j-1]), \
                          A[i-1][j] + GAP, \
                          A[i][j-1] + GAP)
    return A



def trace(A, seqA, seqB):
    ''' produces optimal alignment path through matrix.     
    Args:
    A - matrix
    seqA - string of dna base sequence, labelling rows of matrix A
    seqB - string of dna base sequence, labelling columns of matrix A
    
    Returns: tuple of two strings, the alignment sequence for seqA and 
              seqB respectively.
    
    '''
    #make input seqs lists, and add initial placeholder value
    seqA =  list(seqA) 
    seqB =  list(seqB)
    seqA.insert(0, '+')
    seqB.insert(0, '-')
    
    #initialize output alignment sequences
    a_aligned = ''
    b_aligned = ''
    
    #initialize iterating values i and j to begin at bottom right of matrix
    i = len(seqA) -1
    j = len(seqB) -1
              
    #build output alignment sequences until the trace line is at the top-left
    #    of the matrix
    while i != 0 or j != 0:
        
        #if the trace line is at the top edge of the matrix (i=0), go left
        if i == 0:
            #add gap to seqA alignment
            a_aligned = '-' + a_aligned
            #add base to seqB alignment
            b_aligned = seqB[j] + b_aligned
            j = j - 1
        
        #if the trace line is at the left edge of the matrix (j=0), go up
        elif j == 0:
            a_aligned = seqA[i] + a_aligned
            b_aligned = '-' + b_aligned
            i = i - 1        
       
        #if bases match, trace line goes diagonal, including both bases in 
        #   respective output alignnment sequences
        elif seqA[i] == seqB[j]:
            #add bases to both alignments
            a_aligned = seqA[i] + a_aligned
            b_aligned = seqB[j] + b_aligned
            j = j -1
            i = i -1
            
        #if top-left is greater than left and greater than top:
        elif A[i-1][j-1] > A[i][j-1] and A[i-1][j-1] > A[i-1][j]:
            #add respective base to both aligned seq
            a_aligned = seqA[i] + a_aligned
            b_aligned = seqB[j] + b_aligned
            i = i - 1
            j = j - 1
                            
        #if top greater than top-left and greater than left:
        elif A[i-1][j] > A[i-1][j-1] and A[i-1][j] > A[i][j-1]:
            #add base to A
            a_aligned = seqA[i] + a_aligned
            #add gap to B
            b_aligned = '-' + b_aligned
            i = i - 1
               
       
        #if left  is greater than top-left and greater than top:
        elif A[i][j-1] > A[i-1][j-1] and A[i][j-1] > A[i-1][j]:
            #add respective base to both aligned seq
            a_aligned = '-' + a_aligned
            b_aligned = seqB[j] + b_aligned
            j = j - 1
            
        #if there is some unforeseen combination, this is safety net
        #  it will proceed as if it were a match
        else:
            a_aligned = seqA[i] + a_aligned
            b_aligned = seqB[j] + b_aligned
            i = i -1
            j = j -1
                     
     
        
    #concatenates lists of aligned sequences into strings
    a_aligned = ''.join(a_aligned)
    b_aligned = ''.join(b_aligned)
    
    #returns tuple of two strings: seqA alignment and seqB alignment
    return a_aligned, b_aligned


if __name__ == "__main__":
    P = align2seq("ATCCTGGCG", "ACCCGTCCCG")
    print trace(P, "ATCCTGGCG", "ACCCGTCCCG")
    print hairpin("aatgcatcgxxxcgatgcat", 3)
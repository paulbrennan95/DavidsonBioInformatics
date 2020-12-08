
#Notes:
# NEED TO: 
#            fix alignment scoring system: 2 consecutive mismatches is extra penalty
#            find all grnas, only illustrate those with lowest off-target scores 
#             create comprehensive score for off-target specificity
#              organize the output grnas to line up with their specific region
#              
# 
# Potentially unsolvable problems:
#               change font to courier, or something with same letter distances
#               are windows scrollable?

from graphics import *
from button import *

def find_g_rna(dna):
    ''' Searches for eligble PAM sequences in desired Cas9 cut site dna region 
    Args-
    dna - string of sequence of interest for Cas9 cutting
    
    Returns: list of eligible g_rnas to determine off-target binding
    
    '''
    
    
    #start search for PAMs after base 20      
    pam_loc = 21
    #initialize list of all possible gRNAs in sequence
    grna_seqs = []
    for i in range(dna.count('gg')):
        x = dna.find('gg', pam_loc)
        #if any PAM was found, add its gRNA to list
        if x == -1:
            break
        else:
            pam_loc = x + 1
            #add 20 bases before PAM plus 3 bases of PAM to gRNA for list
            grna_seqs.append(dna[x-21:x+2])
   
    return grna_seqs

def segment(grna):
    seg_list = []
    seg_list.append(grna[:10])
    seg_list.append(grna[10:15])
    seg_list.append(grna[15:20])
    return seg_list

def align_score(site, g_rna):
    
    #inputs need to be 20  bp strings
    if len(site) != 23 or len(g_rna) != 23:
        return
    
    site = segment(site)
    g_rna = segment(g_rna)
    
    
    align_score = 0
    align_score_list = [5, 70, 50]
    
    for i in range(3):
        for j in range(len(site[i])):
            if site[i][j] == g_rna[i][j]:
                    align_score += align_score_list[i]
            else:
                    align_score -= align_score_list[i]
                    
    
  
    #FIXME: include double penalty for consecutive mismatches
    #maximum align score is 650
    
    return align_score




    
#print segment("aaacgactacgaatcgatcgatcactag")

def comp_ref(ref_file, grna):
    with open(ref_file, 'r') as rf:
        #remove fasta header
        rf.readline()
        #make list of 
        ref = rf.readlines()
        
        #remove newline characters
        new_ref = []
        
        for seq in ref:
            if seq[-1] == '\n':            
                seq = seq[0:len(seq)-1]
            new_ref.append(seq)    
                
        #join list of seqs into one string
        ref = ''.join(new_ref)
        
        #make all bases lowercase
        ref = ref.lower()
    
                    
        #need to iterate through ref genome scanning every 20bp sequence for align_scores
        #only record seq with align score greater than half max align score
        #subract off target align_scores from 
        
        align_score_list =[-1000, -1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000, -1000]
        
        for j in range(2):
            #on second time around, check the reverse complement strand
            if j == 1:
                ref = rev(comp(ref))
            for i in range(len(ref)-20):
                if align_score(ref[i:i+23], grna) > min(align_score_list):
                                        
                    align_score_list.append(align_score(ref[i:i+23], grna))
                    align_score_list = sorted(align_score_list)
                    align_score_list.remove(min(align_score_list))
        
        #remove the initial alignment, which is always perfect match
        align_score_list.remove(max(align_score_list))
        #return align_score_list           
        return int(round(100 - sum(align_score_list)/100.0))
      
    
def rev(dna):
    return dna[::-1]

def comp(seq):
    comp_dict = {'g':'c','c':'g','a':'t','t':'a'}
    comp_list = []
    for base in seq:
        comp_list.append(comp_dict[base])
    
    return ''.join(comp_list)
        
def make_menu(w):
    # Make the three buttons
    exit_bl = [92, 92]
    exit_tr = [98, 98]
    exit = Button(exit_bl, exit_tr, w, "Exit")
    
    find_bl = [65, 92]
    find_tr = [90, 98]
    find = Button(find_bl, find_tr, w, "Find Guide RNAs")    
    
    # Make the target sequence input box
    seq_inbox = Entry(Point(15, 95), 20)
    seq_inbox.anchor = Point(15, 95)
    seq_inbox.setText("atgccccgcgatcagctagctgggagctacgactactagctagctagc")
    seq_inbox.setFace("courier")
    seq_inbox.draw(w)
    si_label = Text(Point(12, 98), "Enter Target Sequence:")
    
    si_label.draw(w)    
    
    
    # Make the reference genome input box
    ref_inbox = Entry(Point(1595), 15)
    ref_inbox.anchor = Point(50, 95)
    ref_inbox.setText("short_ecoli.txt")
    ref_inbox.setFace("courier")
    ri_label = Text(Point(50, 98), "Enter Reference Genome File:")
    ri_label.draw(w)    
    ref_inbox.draw(w)      
    
  
    
    # Make the output file input box
    seq_inbox = Entry(Point(15, 95), 20)
    seq_inbox.anchor = Point(15, 95)
    seq_inbox.setText("output_file.txt")
    seq_inbox.setFace("courier")
    seq_inbox.draw(w)
    si_label = Text(Point(12, 98), "Enter Output File Name:")
    
    si_label.draw(w)    
    
    return exit, find, ref_inbox, seq_inbox

    
def main():
    # Make and set up the window
    w = GraphWin("CRISPR Guide RNA Finder", 800, 800)
    w.setCoords(0, 0, 100, 100)
    exit, find, ref_inbox, seq_inbox = make_menu(w)
    
    # Wait for the mouse and then get the input
    #    Need this so user has time to type and then click something BEFORE we grab the text
    p = w.getMouse()    
    
    # Repeat until user clicks the Exit button
    while not exit.clicked(p):
              
        
        if find.clicked(p):
            dna = str(seq_inbox.getText())    
            dna = dna.lower()
            ref_genome = str(ref_inbox.getText())
            #if any(dna) not in ['a','c','g','t']:
                #output = "Error, please enter DNA bases only"
        
            #else:
            if 1 == 1:    
                
            
                print "Input DNA Query:", dna
                output = "Results of CRISPR guide RNA search: Score is (100 - (avg(align_score)/100))\n"
                
                
                #find list of all possible gRNA sites in the users dna
                grnas = find_g_rna(dna) + find_g_rna(rev(comp(dna)))
                
                
                output = output + "\n|%-10s | %-25s | %-15s | %-5s|" %("Strand", "gRNA Sequence", "Location", "Score") + '\n'+ "_"*66
                output_list = []
                for grna in grnas:
                    #if the grna is found on the dna, it is on the direct strand, 
                    #   otherwise, it is on the reverse strand
                    
                   
                    
                    if dna.find(grna) > 0:
                        strand = "direct"
                        dna1 = dna
                    else:
                        strand = "reverse"
                        dna1 = rev(comp(dna))
                        #output includes: strand, gRNA sequence, gRNA location, off-target score
                    output_list.append([strand, grna, dna1.find(grna), comp_ref(ref_genome, grna)])
                   
                output_list.sort(key=lambda x: int(x[3]), reverse = True)
                for grna in output_list:
                    output = output + '\n' + "|%-10s | %-25s | %-15s | %-5s|" % (grna[0], grna[1], grna[2], grna[3])
                    
                with open(output_file, 'w') as new_file:
                    new_file.write(output)
                    
                
                p = w.getMouse() 
                
                # Clear the previous message, so new one will be readable
                info.undraw()                
      
    w.close()
        
#def output(grna_info, query_dna):
    #'''Accepts list of lists containing (query_dna, strand, grna, grna_loc, off_targ_sites)
    
    
    #'''
    ## size of one character is 15 units 
    ##   size of window should be length of query dna plus 15 or so on each side 
    
    ## size of window will be:
    ##                          width = 12 times the length of the query_dna 
    ##                          height = 20 times the number of grnas
    
    #win = GraphWin("Guide RNA Search Results", len(query_dna)*9, 50*len(grna_info) )
    
    #win.setCoords(0,0,100, 100)
    #query_seq = Text(Point(50,50), "Query:   " + query_dna)    
    #query_seq.draw(win)
    
    #d = 0
    #r = 0
    
    #for grna in grna_info:
        #print grna
        #if grna[1] == 'direct':
            #x = Text(Point(50, 53 + d*3 ), grna[2])    
            #x.draw(win)   
            #d = d + 1
        #if grna[1] == 'reverse':
            #print "XXXXX"
            #x = Text(Point(50, 47 - r*3), grna[2])    
            #x.draw(win)   
            #r = r + 1
    
    
    
    
    #if win.getMouse():
        #win.close()
     
       
                
#TEST PRINT Statements
#entry_win = GraphWin("Enter Sequence and Genome File", 200, 200)
#seq_entry = Entry(Point(50, 150), 100)
#file_entry = Entry(Point(50, 50), 100)
#main(seq_entry, file_entry)



main()
#main("GGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCgGGTattacatg", "short_ecoli.txt")
#print comp_ref('short_ecoli.txt', "AATGTCGATCGCCATTATGG")
#print align_score('actgcgcgtacgatcgatc', 'actgcgcgtacgatcgatc')



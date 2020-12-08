#looking for validity of DNA sequence
dna = raw_input('Enter DNA sequence: ')
dna = dna.lower()
x = 0
valid = True
while x < len(dna):
    if dna[x] != 'a' and dna[x] != 'c' and dna[x] != 'g' and dna[x] != 't':
        valid = False
    x = x + 1 
if not valid:
    print 'Warning: DNA sequence is not valid'

#  all 6 reading frames
if valid:   
    # z will be zero for forward strand, 1 for reverse strand
    z = 0
    while z < 2:
        if z == 0:
            orientation = 'direct'
        else: 
            orientation = 'reverse'
            dna = dna[::-1]
        #all 3 orfs on either forward or reverse
        y = 0
        while y <= 2:
            #looking for single orf
            orf1 = False
            i = y
            while i < len(dna) - 2:
                if dna[i:i+3] == 'atg':
                    start_orf1 = i
                    i = len(dna) - 2
                    #if start codon, look for stop codon
                    x = start_orf1 + 3
                    while x < len(dna) - 2:
                        if dna[x:x+3] == 'tag' or dna[x:x+3] == 'tga' or dna[x:x+3] == 'taa':
                            stop_orf1 = x
                            x = len(dna) - 2
                            orf1 = True
                        else:
                            x = x + 3        
                else:
                    i = i + 3
            if orf1:
                print '>ORF number 1 in reading frame %d on the %s strand extends from base' % ((y+1), orientation), \
                      start_orf1 + 1, 'to base', stop_orf1 + 4, '.'
                print dna[int(start_orf1):int(stop_orf1+3)]
            if not orf1:
                print "No ORFs were found in reading frame %d on the %s strand." % ((y+1), orientation) 
            y = y + 1
       
        z = z + 1


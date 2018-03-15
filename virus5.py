#!/usr/bin/env python



import bisect
from itertools import permutations,combinations

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
    

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match    
        
def kset_generator(seq,k):
    kmers = set()
    for idx in range(0,len(seq) - k + 1):
        #print seq[idx:idx + k]
        kmers.add(seq[idx:idx + k])
        
    return kmers

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest    

def pick_maximal_overlap(reads,k):
    reada,readb = None, None
    best_olen = 0
    for a,b in permutations(reads,2):
        olen = overlap(a,b, min_length = k)
        if olen > best_olen:
            reada, readb = a,b
            best_olen = olen
    return reada, readb, best_olen
    
def greedy_scs(reads,k):
    read_a, read_b, olen = pick_maximal_overlap(reads,k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads,k)
    return ''.join(reads)
    
def scs_list(ss):
    """ Returns a list of shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest    
    
def COLLAPSE(dna):
    long_one_length = 0   
    for item in dna:
        if len(item) > long_one_length:
            long_one_length = len(item)   
         
    #print long_one_length
    for item in dna:
        if len(item) == long_one_length:
            longest = item
            
    collapse  = [ x for x in dna if longest.find(x) == 0]
    #print collapse
    collapse2  = [ x for x in dna if len(x) == long_one_length]
    #print
    #print collapse2
    #print
    singletons = [ x for x in dna if x not in collapse]
    #print singletons
    collapse2 = collapse2 + singletons
    #print
    #print collapse2
    return collapse2
    
Sequences, quals = readFastq('ads1_week4_reads.fq')
print "The number of reads is "
print len(Sequences)
print
OVERLAPS = 0  
MATCHES = set()
STUFF ={}  
STUFF2 ={} 
SEQS = []
klen = 50 
SEQS2 = Sequences
for seq in Sequences:
    ks = kset_generator(seq,klen)
    
    for item in ks:
        STUFF.setdefault(item, []).append(seq)  
        
strong_suffix = [s for s in Sequences if len(STUFF.get(s[-klen:])) > 1 ] 
SINGLE = [s for s in Sequences if s not in strong_suffix ]
print
print "This one has no overlapping suffix with any other read "
print SINGLE
print   
print " The number of reads with a strong suffix is "
print len(strong_suffix)
print
Expendable = []

Coalesce  = []
my_list = []
for a in strong_suffix:
    #for b in STUFF.get(a[-klen:])
    #for b in Sequences:
    for b in STUFF.get(a[-klen:]):
        ol = overlap(a,b,klen)
        if ol > 50 and ol < 100:
            new_read = a + b[ol:]
            #print new_read
            my_list.append(new_read)
            #print
    new_ones = COLLAPSE(my_list)

print len(new_ones)
print
print new_ones[0]
print len(new_ones[0])
#new_ones = new_ones  + SINGLE
print
#print new_ones
print len(new_ones)
trunc_new_ones  = []
HUGS = {}
LONG_ONES = []
for item in new_ones:
    HUGS[item] = len(item)
for k,v in HUGS.items():
    print k, v
    if v > 110:
        LONG_ONES.append(k)
        
print
print len(LONG_ONES)
print
print LONG_ONES[0]
print len(LONG_ONES[0])
print

idx = 0

while idx < len(LONG_ONES):
    trunc_new_ones.append(LONG_ONES[idx])
    idx += 2
    
trunc_new_ones = trunc_new_ones  + SINGLE
print "In the array trunc_new_ones there are " + str(len(trunc_new_ones)) + " reads"

# At this point maybe run greedy_scs(reads,k) on the new_ones array
# preassemble some of the reads together

klen = 20
viral = greedy_scs(trunc_new_ones,klen)
print len(viral)
print
print viral.count('A')
print
print viral.count('T')
print
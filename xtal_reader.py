
# author Scott Lew

import re
import csv

class pcrystal():
    def __init__(self, id, seq):
        self.id  = id
        self.seq = seq
        self.label = 0

my_xtals = []
words = []
status = "arrow"
litmus = "diffraction-quality crystals"
litmus2 = "crystallization failed"
xtalfile = open("pxtal.txt","r")
for pxtal in xtalfile.read().split(">"):
    if len(pxtal) > 0:
        my_xtals.append(pxtal)
print()
for item in my_xtals:
    #print(item),
    id_seq = item.split("\n")
    pid = id_seq[0]
    seq = id_seq[1]
    # To deal with strange blank lines in crystals text file
    if len(pid) > 0:
        words = pid.split("#")
    #print(type(pid))
    #print(len(pid))
    if (len(words) > 0):
        #print(words[0])
        pid2 = words[0]
        #print(words[1])
        status = words[1]
        #print(status)
        #print(litmus)
        #print(type(status))
        #xtal = pcrystal(pid2,seq)
        if (re.search(litmus,status)):
            print("Yes!")
            xtal = pcrystal(pid2,seq)
            xtal.label = 1
        elif(re.search(litmus2,status)):
            print("No!")
            xtal = pcrystal(pid2,seq)
        else:
            print("Not sure about this protein")

    if (len(words) > 0):
        with open('xtal_2.csv', "a") as csvfile:
            csvwriter = csv.writer(csvfile,  delimiter=',')
            #my_data = [xtal.id,xtal.seq,xtal.label]
            csvwriter.writerow((xtal.id,xtal.seq,xtal.label))
            print ("Points written sucessfully to file")
        #print(xtal.id,xtal.seq,xtal.label)
    print()

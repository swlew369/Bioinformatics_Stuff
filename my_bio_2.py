
# My first Biopython Script
# Scott Lew
import csv
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#define pcrystal class
class pcrystal():
    def __init__(self, id, seq):
        self.id  = id
        self.seq = seq
        self.label = 0
        self.pI = 0
        self.MW = 0
        self.size = self._size()
        self.phobic = self.percent_phobic()
        self.phillic = self.percent_phillic()

    def _size(self):
        # Borrowed this code from Mr. Joao Henriques
        return len(self.seq)

    def percent_phobic(self):
        phobic = 0
        Amino_Acids = list(self.seq.upper())
        for aa in Amino_Acids:
            if aa in ['L','A','M','V','I','W','F']:
                phobic = phobic + 1
        return 100*float(phobic)/self.size

    def percent_phillic(self):
        phillic = 0
        Amino_Acids = list(self.seq.upper())
        for aa in Amino_Acids:
            if aa in ['D','K','R','S','E','H','N','Q']:
                phillic = phillic + 1
        return 100*float(phillic)/self.size
#print()
#print("Protein Analysis Started")
# open csv file and read in data row by rows
# process each row for pI, MW, etc
with open('xtal_2.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        pid = row[0]
        pseq = row[1].upper()
        pseq = pseq.replace('X','G')
        xtal = pcrystal(pid,pseq)
        xtal.label = int(row[2])
        analysed_seq = ProteinAnalysis(pseq)
        xtal.pI = analysed_seq.isoelectric_point()
        xtal.MW = analysed_seq.molecular_weight()
        #print()
        #print(pseq)
        #print(xtal.id)
        #print(xtal.seq.upper())
        #print(xtal.label)
        #print ('protein isoelectric point is %.2f' % xtal.pI)
        #print ('protein molecular weight is %.2f' % xtal.MW)
        #print ('protein length is ', xtal.size)
        #print ('percent hydrophobic is %.2f' % xtal.phobic)
        #print ('percent hydrophillic is %.2f' % xtal.phillic)
        #if (xtal.label == 0):
            #print("Crystallization Failed!")
        #else:
            #print("Protein Crystallized!")
        #print()
        with open('xtal_output2.csv', "a") as csvfile:
            csvwriter = csv.writer(csvfile,  delimiter=',')
            csvwriter.writerow((xtal.id,xtal.seq,xtal.label,xtal.pI,xtal.MW,xtal.size,xtal.phobic,xtal.phillic))
print()
print("Protein Analysis Finished")

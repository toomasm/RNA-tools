import pysam
from Bio import SeqIO

""" Code that divides the reads of an sam file between two new files according to nucleotide following
the mapped three prime end of the read. If the nucleotide is not A the read is true and is saved in the
correct reads file. If the nucleotide is A we have no means to know if the read has been correctly trimmed or not
and it is saved in "shady" reads file."""

#Opens the input files. SAM for mapped reads and fasta for template.
sam = pysam.Samfile(raw_input("Enter the sam file path: "), "r")
for fasta in SeqIO.parse(raw_input("Enter the path for your template sequence here: "), "fasta"):
    pass
#output = raw_input("Give a name for output sam file: ")

#output_ok = pysam.Samfile( "ok_" + output, "w", template = sam)
#output_shady = pysam.Samfile("shady_" + output, "w", template = sam)
#The output files between witch the data is divided.
output_ok = pysam.Samfile( "ok_.sam" , "w", template = sam)
output_shady = pysam.Samfile("shady_.sam" , "w", template = sam)


print fasta[-1]
def Sam_3prime_checker(sam_file):
#Code that divides the three prime seq reads according to nucleotide following the three prime end. 
#bacterial genome is circular and code takes this in account. If the read is negative and you map a read 
#to start position of the template file then nucleotide downstream is the last position in template. If read is positive
#and its three prime end is the last position of the template, the following nucletodie will be first nucleotide of the template.
    for read in sam: 
        if read.is_reverse == True: #Checks wether the read has positive or negative orientation
            if read.pos == 0: #Takes into account the circularity of bacterial DNA
                if fasta[-1] == "A":
                    output_shady.write(read)
                else:
                    output_ok.write(read)
            elif fasta[read.pos - 1] == "A":
                output_shady.write(read)
            else:
                output_ok.write(read)
        else:
            if (read.pos+ read.qlen) == len(fasta):
                if fasta[0] == "A":
                    output_shady.write(read)
                else:
                    output_ok.write(read)
            elif fasta[read.pos + read.qlen] == "A":
                output_shady.write(read)
            else:
                output_ok.write(read)
    print "Code has finished!"

Sam_3prime_checker(sam)            

sam.close()
output_ok.close()
output_shady.close()

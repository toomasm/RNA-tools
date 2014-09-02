from Bio import SeqIO

""" Program that trims all "A" nucleotides from the 3 prime end of the fastq reads. It removes "A" nucleotides till it reaches a 
none A nucleotide."""
#/media/toomas/DATA1/fastq test/3prime_fastq_test.fastq

#The input data path entry and creation of output file for trimmed reads.
path_in = raw_input("Enter the fastq file path ") 
#data_input = open(path_in, "r")
path_out = raw_input("Enter the output file name ")
output = open(path_out, "w")

#This code at first parses fastq file. Then it starts to iterate backwards over every characater in the read.
#It counts the number of A nucleotides at the end of the read, removes as many from the 3 prime end and writes the read 
#to a new file.

for rec in SeqIO.parse(path_in, "fastq"):#iterates over reads.
    count = 0
    for c in rec.seq[::-1]: #Iterates backwards over nucleotides in read.
        if c == "A": #counts the number of "A"-s at the 3 prime end.
            count += 1
        else:
            break
    if count == 0: #assures that the read wont be deleted if no A-s are at the end.
        pass
    else:
        rec = rec[:-count] #removes nucleotides from 3 prime end. if -count is 0 it deletes the read so it has to be atlest 1.
    SeqIO.write(rec, output, "fastq")

#data_input.close() 
output.close()

num_lines = sum(1 for line in open(path_in))
print "Input file has {0} lines".format(num_lines)
num_lines = sum(1 for line in open(path_out))
print "Output file has {0} lines".format(num_lines)

from Bio import SeqIO
import sys

""" Program that trims your specifified amount of nucleotides from the 5 prime end of the fastq reads. """

#The input data path entry and creation of output file for trimmed reads.
path_in = raw_input("Enter the fastq file path ") 
data_input = open(path_in, "r")
path_out = raw_input("Enter the output file name ")
output = open(path_out, "w")

#Asks you how much nucleotides would you like to trim from 5 prime end (enter an int) and checks if your entry is valid.
try:
    trimming_number = int(raw_input("How many nucleotides would you like to trim from 5 prime end? Please enter a single number: "))
except ValueError:
    print "This is not a single number!"
    sys.exit()

#Just says what the program is about to do. {0} is a placeholder for info behind .format.
print "This code will now trim {0} nucleotides from the five prime end of all your reads".format(trimming_number)

#This code parses fastq file and then iterates over every separate read removing x-nucleotides from the 5 prime end.
for rec in SeqIO.parse(data_input, "fastq"):
    rec = rec[trimming_number:] #removes needed amount of nucleotides from 5 prime end.
    SeqIO.write(rec, output, "fastq")                              
  
data_input.close()
output.close()

#shows how many reads there were in the input file and how many reads are in the output file. 

num_lines = sum(1 for line in open(path_in))
print "Input file has {0} lines".format(num_lines)
num_lines = sum(1 for line in open(path_out))
print "Output file has {0} lines".format(num_lines)
print "code has finished!"  
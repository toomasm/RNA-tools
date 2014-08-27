import collections
import csv
import pysam

def get_five_prime_positions(bam_file):
    '''Algorithm that makes a list that contains five prime end position of every read.
    
    Checks if the read is mappd to + or - strand. 
    
    If it is on + strand it adds the start position of the read to the list and adds 1 to it 
    (because start position is given in 0 based coordinates).
    
    When the read is on - strand it adds the end position to the list
    (not adding 1 this time because this coordinate is 1 based).

    When the strand information is neither "-" or "+" skips the read.'''

    five_prime_positive = []
    five_prime_negative = []
    for read in bam_file:
        if read.alen == read.rlen and read.is_unmapped == False:
            if read.is_reverse == False:
                five_prime_positive.append(read.pos+1)
            elif read.is_reverse == True:
                five_prime_negative.append(read.pos + read.rlen)
    #print('get_five_prime_positions found {} five prime ends.'.format(len(five_prime)))

    #The part where we count how many five prime ends we have in specific genome positions.
    five_prime_count_positive = collections.Counter(five_prime_positive)
    five_prime_count_negative = collections.Counter(five_prime_negative)
    #print('The counter found {} five prime end mapping positions.'.format(len(five_prime_count)))

    #write the results to csv file.
    with open(raw_input("Name of your csv output file "), "w") as output:
        writer = csv.writer(output)
        for key, value in five_prime_count_positive.items():
            writer.writerow([key, value, "p"])
        for key, value in five_prime_count_negative.items():
            writer.writerow([key, value, "-"])

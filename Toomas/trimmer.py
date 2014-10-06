from Bio import SeqIO
import sys
import pysam



def read_trimmer(input_filename, output_file, trim_from_5p_end, Flag_3prime):
    """
    Program that trims your specified amount of nucleotides from the 5 prime end of the fastq reads.

    This code also allows a step for polyA tail trimming.
    If Flag "Flag_3prime" is set, it trims all "A" nucleotides from the 3' end of the fastq reads.
    It removes "A" nucleotides till it reaches a none A nucleotide.
    If the seq quality of 3' end is of bad quality, it checks the nucleotide following the bad seq quality nucleotides,
    if it is "A", continue trimming. If not "A"- disregard the read.
    Exception is made if only one bad quality nucleotide is located at the three prime end.
    """

    #Prints out what operatios the code will perfom.
    if trim_from_5p_end >= 1 and Flag_3prime == False:
        print "Trimming {0} nucleotides from 5' end.".format(trim_from_5p_end)
    elif not trim_from_5p_end >= 1 and Flag_3prime == True:
        print "Trimming poly A tail from 3' end."
    elif trim_from_5p_end >= 1 and Flag_3prime == True:
        print "Trimming {0} nucleotides from 5' end and A nucleotides (possible polyA tail) from 3' end".format(trim_from_5p_end)

    #Checks the quality of 3' end.
    def three_prime_quality_assessor(fastq_read):
        
        quality =(fastq_read.letter_annotations["phred_quality"])
        count_quality = 0 #Bad quality nucleoties count in 3' end
        #Counts how many nucleotides in 3 prime end have low seq quality.
        for nucleotide_score in quality[::-1]:
            if nucleotide_score <= 2:
                count_quality += 1
            else:
                break
        return count_quality
        
    #Counts the "A" nucleotides in three prime end.    
    def three_prime_A_counter(fastq_read):

        count_nucleotide = 0  #"A"  nucleotide count in 3' end
        for c in fastq_read.seq[::-1]: #Iterates backwards over nucleotides in read.
            if c == "A": #counts the number of "A"-s at the 3 prime end.
                count_nucleotide += 1
            else:
                break
        return count_nucleotide
        
    #Removes "A" nucleotides from the three prime end.
    def polyA_remover(rec, count_A):

        if count_A == 0: #Assures that the read wont be deleted if no A-s are at the end.
            rec = rec #No nucleotide is removed and read stays as is.
        else:
            rec = rec[:-count_A] #Removes nucleotides from 3 prime end. If "A" count is 0 it deletes the read so it has to be atleast 1.

        if len(rec) < 20: #If read length after trimming is shorter then 20 makes a flag True to not write the records to output file.
            too_short_read = bool(True)
        else:
            too_short_read = bool(False)
            
        return rec, too_short_read
            
    #Trims the reads according to input.
    def trimming(input_filename, output_file, trim_from_5p_end):
        #Some variables for counting the reads added and removed from.
        count_reads_added = 0
        count_reads_not_added = 0
        count_PolyA_after_badq = 0
        for rec in SeqIO.parse(input_filename, "fastq"):
            if trim_from_5p_end >= 1:
                rec = rec[trim_from_5p_end:] #removes needed amount of nucleotides from 5 prime end.

            if Flag_3prime == True:

                #Counts how many nucleotides in 3 prime end have low seq quality.
                count_q = three_prime_quality_assessor(rec)
                #If seq quality is good commense "A" trimming.
                if count_q == 0:
                    #Counts how many "A" nucleotides there are in 3' end.
                    count_A = three_prime_A_counter(rec)
                    (rec, do_not_add_read) = polyA_remover(rec, count_A)

                #If more than one nucleotides in the 3' end have a bad seq quality it is impossible to be sure about the polyA tail.
                #These reads are left out.
                #If "A" follows the bad quality nucleotides, trim the bad quality nucleotides and all following "A"-s in row.
                #Save read as shady. 
                elif count_q >= 1 and (str(rec.seq[-(count_q+7):-(count_q)]) == "AAAAAAA"):
                    #Trims bad quality nucleotides from three prime end.
                    rec = rec[:-count_q]
                    count_PolyA_after_badq += 1
                    #Counts how many "A" nucleotides there are in 3' end.
                    count_A = three_prime_A_counter(rec)
                    (rec, do_not_add_read) = polyA_remover(rec, count_A)
                else:#Sets the flag for not adding a read if read does not meet any of the previously set conditions.
                    #For example "A" does not follow a faulty nucleotide.
                    do_not_add_read = bool(True)
                    #Checks whether to add read or not. This flag will not add read to output if it only composes of polyA tail.

                if do_not_add_read == True:
                    count_reads_not_added += 1
                    continue
            count_reads_added += 1        
            SeqIO.write(rec, output_file, "fastq")
    
        print 'Reads written to output: {0}'.format(count_reads_added)
        print 'Reads disgarded: {0}'.format(count_reads_not_added)
        print 'Possible reasons: too short length after trimming or bad seq quality in three prime end outside of polyA tail.'
        print 'Poly A after bad quality nucleotides: {0}'.format(count_PolyA_after_badq)
        print 'Trimming has finished!'

    with open(input_filename,'r') as input_filename:
        with open(output_file, 'w') as output_file:
            trimming(input_filename, output_file, trim_from_5p_end)
    
def three_prime_checker(bam_file, ref_genome_fasta, ok_reads_bam_filename, shady_reads_bam_filename):
    """
    Code that divides the reads of a SAM file between two new files according to nucleotide following
    the mapped three prime end of the read. If the nucleotide is not A the read is true and is saved in the
    correct reads file. If the nucleotide is A we have no means to know if the read has been correctly trimmed or not
    and it is saved in "shady" reads file.
    """
    
    """
    Code that divides the three prime seq reads according to nucleotide
    following the three prime end.
    Bacterial genome is circular and code takes this in account.
    If the read is negative and you map a read to start position of the template file then nucleotide downstream is the last
    position in template. If read is positive and
    its three prime end is the last position of the template, the following nucleotide will be first nucleotide of the template.
    """

    def is_shady(read, fasta):

        if read.is_reverse:
            test, indexA, indexB = read.pos == 0, -1, read.pos -1
        else:
            test, indexA, indexB = (read.pos+ read.qlen) == len(fasta), 0, read.pos + read.qlen

        if test: #Takes into account the circularity of bacterial DNA
            if str(fasta.seq[indexA]) == "A":
                return True
        elif str(fasta.seq[indexB]) == "A":
            return True
            
        return False
                  
    for fasta in SeqIO.parse(ref_genome_fasta, "fasta"):
        pass        

    #Opens the input files. SAM for mapped reads and fasta for template.
    with pysam.Samfile(bam_file, "rb") as bam:
        with pysam.Samfile(ok_reads_bam_filename , "wb", template = bam) as output_ok:
            with pysam.Samfile(shady_reads_bam_filename , "wb", template = bam) as output_shady:
                for read in bam:
                    #If read is of expected length write it as ok read.
                    #Else check for polyA tail on bad nucleotides.
                    if read.rlen < 99 and is_shady(read, fasta):
                        output_shady.write(read)
                    else:
                        output_ok.write(read)
    print "3 prime divider has finished!"
    return True


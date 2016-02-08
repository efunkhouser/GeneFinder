# -*- coding: utf-8 -*-
"""
Here it is! Week one (first six functions) complete.
I added at least one unit test to all except find_all_ORFs and find_all_ORFs_both_strands, as I found those unit tests sufficient.

@author: ELEANOR FUNKHOUSER

"""

import random
import doctest
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """
    if nucleotide == "A": #Pairs: A and T, C and G (2nd strand has complement of base)
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "G":
        return "C"
    else:
        print "That's not a nucleotide or is the wrong format; get_complement is returning None" #Just in case
    

doctest.run_docstring_examples(get_complement, globals())

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("ACCTTTGGGG")
    'CCCCAAAGGT'
    """
    revComplement = "" #initialize as empty string - will be adding onto this in the loop
    for index in range(len(dna)):
        index = index + 1
        compNucleotide = get_complement(dna[-1*index]) #Get the complement of the reverse (5') end
        revComplement = revComplement + compNucleotide #Add the complmenent nucleotide onto the reverse string
    return revComplement

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("TAACG")
    ''
    """
    i = 0 #initialize index to count bases
    dnaString = "" #initialize DNA; will add on 3 at a time in while loop
    while (i+2)<len(dna):
        codon = dna[i:(i+3)] #reading codon by codon
        if codon in ['TAG','TAA','TGA']:
            return dnaString #does not include the stop codon, just everything leading up to it
            break
        dnaString = dnaString + codon
        i = i + 3
    return dna #If no stop codon is found in frame

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("CATGCATTGATG")
    ['ATCCAT']
    """
    startCodon = 'ATG'
    i = 0 #For counting bases
    ORF_list = [] #Initializing the list that this function will return

    while (i+3)<=len(dna):
        codon = dna[i:(i+3)]
        if codon==startCodon:
            oneORF = rest_of_ORF(dna[i:]) #Get the ORF up to but not including the stop codon
            ORF_list.append(oneORF)
            i = i + len(oneORF) #On to the next ORF - no nested shit here
            ### NOTE: This resetting of i assumes that the dna given to find_all_ORFs_oneframe DOES HAVE A STOP CODON SOMEWHERE
        else:
            i = i+3 #Else, look at the next codon in the frame
    return ORF_list




def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ans = []
    ORF_1 = find_all_ORFs_oneframe(dna[:]) #returns list: ORFs in this frame
    ORF_2 = find_all_ORFs_oneframe(dna[1:]) #ORFs in 2nd frame
    ORF_3 = find_all_ORFs_oneframe(dna[2:]) #ORFs in 3rd frame

    ans = ans + ORF_1 + ORF_2 + ORF_3

    return ans


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    forwardDNA = dna
    reverseDNA = get_reverse_complement(dna) #Gives the 3'-5' other strand: bases reversed and complemented

    ans = [] #return list

    ans = ans + find_all_ORFs(forwardDNA) + find_all_ORFs(reverseDNA) #add to the list all ORFs from the original & opposite strands

    return ans


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    listofAll = find_all_ORFs_both_strands(dna)
    longestORF = ''
    for dnaString in listofAll:
        if len(dnaString) > len(longestORF):
            longestORF = dnaString #reassign
    return longestORF
###TODO: write another unit test here



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longestString = ''
    for trial in range(num_trials):
        longest_inTrial = longest_ORF(shuffle_string(dna))
        if len(longest_inTrial) > len(longestString):
            longestString = longest_inTrial
    return len(longestString)
### TODO: Figure out validation that isn't unit testing


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    i = 0
    aaString = ''
    while (i+3) <= len(dna):
        codon = dna[i:i+3]
        amino_acid = aa_table[codon]
        aaString += amino_acid
        i = i+3
    return aaString


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    minLength = longest_ORF_noncoding(dna,num_trials=1500)
    possibleORFs = find_all_ORFs_both_strands(dna) #outputs a list
    aa_List = []

    for each_ORF in possibleORFs:
        if len(each_ORF) >= minLength:
            aa_List.append(coding_strand_to_AA(each_ORF))

    if aa_List == []:
        print "Sorry, this DNA strand appears to be a piece of junk."
    else:
        return aa_List



if __name__ == "__main__":
    import doctest
    doctest.testmod()

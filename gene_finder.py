# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

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
    """
    if nucleotide == "A":
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "G":
        return "C"
    else:
        print "That's not a nucleotide or is the wrong format; get_complement is returning None"
    

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
    """
    revComplement = "" #initialize as empty string - will be adding onto this in the loop
    for index in range(len(dna)):
        index = index + 1
        compNucleotide = get_complement(dna[-1*index]) #Get the complement of the reverse (5') end
        revComplement = revComplement + compNucleotide
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
    """
    i = 0 #initialize index counter
    dnaString = "" #initialize DNA; will add on 3 at a time in while loop
    while (i+2)<len(dna):
        codon = dna[i:(i+3)] #reading codon by codon
        if codon=='TAG' or codon=='TAA' or codon=='TGA':
            return dnaString
            break
        dnaString = dnaString + codon #if not stop, add it to the return string
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
    """
    startCodon = 'ATG'
    i = 0
    ORF_list = []
    while (i+3)<=len(dna):
        codon = dna[i:(i+3)]
        if codon==startCodon:
            oneORF = rest_of_ORF(dna[i:])
            ORF_list.append(oneORF) #Add the ORF to the list
            i = i + len(oneORF) #On to the next one - no nested shit here
            ### NOTE: This resetting of i assumes that the dna given to find_all_ORFs_oneframe DOES HAVE A STOP CODON SOMEWHERE
        else:
            i = i+3
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
    ORF_1 = find_all_ORFs_oneframe(dna[:]) #returns list
    ORF_2 = find_all_ORFs_oneframe(dna[1:])
    ORF_3 = find_all_ORFs_oneframe(dna[2:])

    ans = ans+ORF_1 #repetitive for debuggin purposes
    ans = ans+ORF_2 #i'ma leave it for clarity's sake though
    ans = ans+ORF_3

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
    reverseDNA = get_reverse_complement(dna)

    ans = []

    ans = ans + find_all_ORFs(forwardDNA)
    ans = ans + find_all_ORFs(reverseDNA)

    return ans


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()

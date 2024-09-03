#!/usr/bin/env python

# Author: Amelia Dayton <daytonamelia@gmail.com>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.4x"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')
codon_dict = {
#       xTx           xCx           xAx           xGx
# Txx
    'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys',     # TxT
    'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys',     # TxC
    'TTA': 'Leu', 'TCA': 'Ser', 'TAA': '---', 'TGA': '---',     # TxA
    'TTG': 'Leu', 'TCG': 'Ser', 'TAG': '---', 'TGG': 'Trp',     # TxG
# Cxx
    'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg',     # CxT
    'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',     # CxC
    'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',     # CxA
    'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',     # CxG
# Axx
    'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser',     # AxT
    'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',     # AxC
    'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',     # AxA
    'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',     # AxG
# Gxx
    'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly',     # GxT
    'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',     # GxC
    'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',     # GxA
    'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'      # GxG
}

def convert_phred(letter: str) -> int:
    """Converts a single character into a +33 phred score"""
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Returns average quality score of input phred string as a float"""
    phred_sum = 0
    for i, score in enumerate(phred_score):
        phred_sum += convert_phred(score)
    return phred_sum/len(phred_score)

def validate_base_seq(seq:str, RNAflag:bool=False) -> bool:
    '''Returns True if case-insensitive string is composed of only As, Ts (or Us if RNAflag), Gs, Cs.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    lst_len = len(lst)
    if lst_len%2==1: # odd length list
        median_val = lst[int(lst_len/2 - 0.5)]
    else: # even length list
        median_val = (lst[int(lst_len/2 - 0.5)]+lst[int(lst_len/2 + 0.5)])/2
    return median_val

def revcomp(DNA:str) -> str:
    '''Returns the reverse complement of a DNA sequence.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used the proper sequence?"
    DNAtable = str.maketrans("ATCG", "TAGC")
    return DNA[::-1].translate(DNAtable)

def oneline_fasta(infile:str, outfile: str) -> None:
    '''Given a fasta file wth reads that span multiple lines, write a new fasta file with reads that only span oneline.'''
    first_line = True
    with open(infile, "r") as rf, open(outfile, "w") as wf:
        for line in rf:
            line = line.strip("\n")
            if line.startswith(">"):
                if first_line:
                    wf.write(line+"\n")
                    first_line = False
                else:
                    wf.write("\n"+line+"\n")
            else:
                wf.write(line)

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Passed convert_phred tests!")

    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    assert qual_score("II") == 40
    print("Passed qual_score tests!")

    assert validate_base_seq("ATCGATATAT") == True
    assert validate_base_seq('ABCDEF') == False
    assert validate_base_seq('AUCGUUCUGA') == False
    assert validate_base_seq('AUCGUUCUGA', True) == True
    print("Passed validate_base_seq tests!")

    assert gc_content("CGCG") == 1
    assert gc_content("ATGC") == 0.5
    assert gc_content("ATATATATA") == 0.0
    print("Passed gc_content tests!")

    assert calc_median([1,2,4,5]) == 3
    assert calc_median([1,3,7]) == 3
    print("Passed calc_median tests!")

    assert revcomp("ATCG") == "CGAT"
    assert revcomp("ATATAT") == "ATATAT"
    assert revcomp("GCGCAT") == "ATGCGC"
    print("Passed revcomp tests!")
    
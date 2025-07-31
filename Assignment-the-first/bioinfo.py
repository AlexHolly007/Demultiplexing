#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

#import test_validate_base_seq, test_calc_median, test_gc_content, test_qual_score, test_convert_phred


__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    '''Write your  own doc string'''
    return sum(convert_phred(x) for x in phred_score) / len(phred_score)


DNAbases = set('ATGCatgc')
RNAbases = set('AUGCaugc')
def validate_base_seq(DN_input:str, RNAflag:bool = False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(DN_input).issubset((RNAbases if RNAflag else DNAbases))


def gc_content(DNA):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"

    return (DNA.count("G")+DNA.count("C"))/len(DNA)


def calc_median(lst: list):
    '''Given a sorted list, returns the median value of the list'''
    length = len(lst)
    return (lst[(length//2)-1] + lst[(length//2)]) / 2 if ((length%2) == 0) else lst[(length//2)]


def oneline_fasta(inputfile: str, outfile: str):
    '''docstring'''
    with open(inputfile, "r") as fi, open(outfile, "w") as fo:
        first_line_check = True
        for line in fi:
            if line.startswith(">"):
                if first_line_check:
                    fo.write(line)
                else:
                    fo.write(f"\n{line}")
            else:
                fo.write(line.strip())

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    # test_calc_median.test_calc_median()
    # test_convert_phred.test_convert_phred()
    # test_gc_content.test_gc_content()
    # test_qual_score.test_qual_score()
    # test_validate_base_seq.test_validate_base_seq()
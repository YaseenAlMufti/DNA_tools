# DNA_tools
DNA Analysis tools written in Python 3

# Initialization
```
import DNAtoRNA
dtr = DNAtoRNA.DNA_Analysis()
```
# Methods & Usage

### DNA_to_RNA(seq)
This method converts a DNA Sequence (provided as a string) to an RNA sequence string. This method removes any whitespaces and/or empty lines.                                                                                                    

Example:
```
inp = "AAAAGCGCCGCGCGAAATTTTCCCGCGCGCTTTTTCGCTCGCGCTCGCGCTCGCCGT"
print(dtr.DNA_to_RNA(inp))
//output: AAAAGCGCCGCGCGAAAUUUUCCCGCGCGCUUUUUCGCUCGCGCUCGCGCUCGCCGU
```

### RNA_to_Amino(RNAseq, reading_Frame)
This method converts an RNA sequence(provided as a string) to three letter Amino Acid sequence returning the output and the number of Amino Acids in the sequence.

The reading frame can be customized using the second argument: reading_Frame = 1(default) or 2 or 3.

Example:
```
print(dtr.DNA_to_Amino(dtr.DNA_to_RNA(inp), 1))
//output: Lys-Ser-Ala-Ala-Arg-Asn-Phe-Pro-Ala-Arg-Phe-Phe-Ala-Arg-Ala-Arg-Ala-Arg-Arg 19
```

### three_to_One(AminoSeq)
This method translates a three letter Amino Acid sequence to one letter Amino Acid sequence.

Example:
```
print(dtr.three_to_One(dtr.RNA_to_Amino(dtr.DNA_to_RNA(inp), 1)))
//output: KSAARNFPARFFARARARR 19
```

### get_three_to_five(sequence)
Returns the reverse of a sequence 5' to 3'.

Example:
```
print(dtr.get_three_to_five(inp))
//output: TGCCGCTCGCGCTCGCGCTCGCTTTTTCGCGCGCCCTTTTAAAGCGCGCCGCGAAAA
```

### get_complementary(sequence)
Returns the complementary strands of the given strand (DNA or RNA sequence can both be used).

Example:
```
print(dtr.get_complementary(inp))
//output: TTTTCGCGGCGCGCTTTAAAAGGGCGCGCGAAAAAGCGAGCGCGAGCGCGAGCGGCA
```

### find_CpG(DNAsequence)
Returns the ranges where CpG islands might be found in a large DNA sequence. This method uses the sliding window algorithm to detect potential CpG islands.

Window length and cutoff configs:

window length: 200 bases

gc content cutoff: > 50%

observed/expected ratio cutoff: > 0.60

Example:
```
with open("sequence.txt", "r") as fin:
    r = str(fin.read())
    with open("islands.txt", "w") as fout:
        fout.writelines(item + "\n" for item in dtr.find_CpG(r))
//input: sequence.txt (contains a large DNA sequence(provided as an example in files))
//output: islands.txt (contains a large data of ranges of potential CpG islands)
```

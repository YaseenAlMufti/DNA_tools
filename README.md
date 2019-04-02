# DNA_tools
DNA Analysis tools written in Python 3

# Basic Usage
import DNAtoRNA
dtr = DNAtoRNA.DNA_Analysis()

# Method DNA_to_RNA(seq)
This method converts a DNA Sequence (provided as a string) to an RNA sequence string. This method removes any whitespaces and/or empty lines.                                                                                                    

Example:
inp = "AAAAGCGCCGCGCGAAATTTTCCCGCGCGCTTTTTCGCTCGCGCTCGCGCTCGCCGT"
print(dtr.DNA_to_RNA(inp))
//output: AAAAGCGCCGCGCGAAAUUUUCCCGCGCGCUUUUUCGCUCGCGCUCGCGCUCGCCGU

# Method RNA_to_Amino(RNAseq, reading_Frame)
This method converts an RNA sequence(provided as a string) to three letter Amino Acid sequence

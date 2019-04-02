
import re


class DNA_Analysis:
    amino_acids = {
        "uuu": "Phe-",
        "uuc": "Phe-",
        "uua": "Leu-",
        "uug": "Leu-",
        "cuu": "Leu-",
        "cuc": "Leu-",
        "cua": "Leu-",
        "cug": "Leu-",
        "auu": "Ile-",
        "auc": "Ile-",
        "aua": "Ile-",
        "aug": "Met-",
        "guu": "Val-",
        "guc": "Val-",
        "gua": "Val-",
        "gug": "Val-",
        "ucu": "Ser-",
        "ucc": "Ser-",
        "uca": "Ser-",
        "ucg": "Ser-",
        "ccu": "Pro-",
        "ccc": "Pro-",
        "cca": "Pro-",
        "ccg": "Pro-",
        "acu": "Thr-",
        "acc": "Thr-",
        "aca": "Thr-",
        "acg": "Thr-",
        "gcu": "Ala-",
        "gcc": "Ala-",
        "gca": "Ala-",
        "gcg": "Ala-",
        "uau": "Tyr-",
        "uac": "Tyr-",
        "uaa": "Stop-",
        "uag": "Stop-",
        "cau": "His-",
        "cac": "His-",
        "caa": "Gin-",
        "cag": "Gin-",
        "aau": "Asn-",
        "aac": "Asn-",
        "aaa": "Lys-",
        "aag": "Lys-",
        "gau": "Asp-",
        "gac": "Asp-",
        "gaa": "Glu-",
        "gag": "Glu-",
        "ugu": "Cys-",
        "ugc": "Cys-",
        "uga": "Stop-",
        "ugg": "Trp-",
        "cgu": "Arg-",
        "cgc": "Arg-",
        "cga": "Arg-",
        "cgg": "Arg-",
        "agu": "Ser-",
        "agc": "Ser-",
        "aga": "Arg-",
        "agg": "Arg-",
        "ggu": "Gly-",
        "ggc": "Gly-",
        "gga": "Gly-",
        "ggg": "Gly-"
    }

    amino_Acids_One_Letter = {
        'Trp-': "W",
        'Ala-': "A", 
        'Asn-': "N", 
        'Phe-': "F", 
        'Met-': "M", 
        'Arg-': "R", 
        'Gly-': "G", 
        'Glu-': "E", 
        'His-': "H", 
        'Ser-': "S", 
        'Lys-': "K", 
        'Leu-': "L", 
        'Asp-': "D", 
        'Thr-': "T", 
        'Tyr-': "Y", 
        'Ile-': "I", 
        'Stop-': "*", 
        'Pro-': "P", 
        'Gin-': "Q", 
        'Cys-': "C", 
        'Val-': "V"
    }

    def DNA_to_RNA(self, DNA_Sequence):
        r = "".join(DNA_Sequence.split())
        return re.sub("t", "u", r.lower()).upper()


    def RNA_to_Amino(self, RNA_Sequence, reading_Frame=1):
        if reading_Frame == 1:
            patternreg = re.compile("|".join(map(re.escape, self.amino_acids.keys())))
            r = patternreg.sub(
                lambda match: self.amino_acids[match.group(0)], RNA_Sequence.lower())
            if r[-1] == "-":
                r = r[:-1]
            c = r.split("-")
            # g = re.sub("Stop-", "Stop\n", r)
            return r + " " + str(len(c))
        elif reading_Frame == 2:
            RNA_Sequence = RNA_Sequence[1:]
            patternreg = re.compile("|".join(map(re.escape, self.amino_acids.keys())))
            r = patternreg.sub(
                lambda match: self.amino_acids[match.group(0)], RNA_Sequence.lower())
            if r[-1] == "-":
                r = r[:-1]
            c = r.split("-")
            # g = re.sub("Stop-", "Stop\n", r)
            return r + " " + str(len(c))
        elif reading_Frame == 3:
            RNA_Sequence = RNA_Sequence[2:]
            patternreg = re.compile("|".join(map(re.escape, self.amino_acids.keys())))
            r = patternreg.sub(
                lambda match: self.amino_acids[match.group(0)], RNA_Sequence.lower())
            if r[-1] == "-":
                r = r[:-1]
            c = r.split("-")
            # g = re.sub("Stop-", "Stop\n", r)
            return r + " " + str(len(c))


    def three_to_One(self, Amino_Sequence):
        Amino_Sequence = re.sub(r'[0-9]+', '', Amino_Sequence)
        Amino_Sequence = Amino_Sequence[:-1] + "-"
        patternreg = re.compile(
            "|".join(map(re.escape, self.amino_Acids_One_Letter.keys())))
        r = patternreg.sub(
            lambda match: self.amino_Acids_One_Letter[match.group(0)], Amino_Sequence)
        if r[-1] == "-":
            r = r[:-1]
        return r + " " + str(len(r))


    def find_CpG(self, Sequence):
        r = "".join(Sequence.split())
        if "t" in r:
            r = self.DNA_to_RNA(r)
        bases = len(r)
        return self.get_islands(r, bases)


    def get_islands(self, r, Nbases):
        istart = 0
        CpG_islands = []
        while istart + 201 <= Nbases:
            window_length = 200
            c_content = len(re.findall("C", r[istart:istart + 201]))
            g_content = len(re.findall("G", r[istart:istart + 201]))
            gc_content = round(
                (((c_content + g_content) / window_length) * 100), 2)
            observed_gc = len(re.findall("CG", r[istart:istart + 201]))
            expected_gc = (((c_content + g_content) / 2)**2) / window_length
            observed_expected_ratio = round((observed_gc / expected_gc), 2)
            if gc_content > 50 and observed_expected_ratio > 0.60:
                CpG_islands.append(
                    f"CpG island at region {istart + 1} to {istart + 200} (Obs/Exp = {observed_expected_ratio} and %GC = {gc_content})")
            print(c_content, " ", g_content, " ", gc_content, " ",
                observed_gc, " ", expected_gc, " ", observed_expected_ratio)
            istart += 1
        print(len(CpG_islands))
        return CpG_islands


    def get_complementary(self, Sequence):
        Sequence = "".join(Sequence.split()).upper()
        complementary_Sequence = ""
        if "U" in Sequence:
            for char in Sequence:
                if char == "A":
                    complementary_Sequence = complementary_Sequence + "U"
                elif char == "U":
                    complementary_Sequence = complementary_Sequence + "A"
                elif char == "C":
                    complementary_Sequence = complementary_Sequence + "G"
                elif char == "G":
                    complementary_Sequence = complementary_Sequence + "C"
                else:
                    complementary_Sequence = complementary_Sequence + char
        else:
            for char in Sequence:
                if char == "A":
                    complementary_Sequence = complementary_Sequence + "T"
                elif char == "T":
                    complementary_Sequence = complementary_Sequence + "A"
                elif char == "C":
                    complementary_Sequence = complementary_Sequence + "G"
                elif char == "G":
                    complementary_Sequence = complementary_Sequence + "C"
                else:
                    complementary_Sequence = complementary_Sequence + char
        return complementary_Sequence


    def get_three_to_five(self, Sequence):
        Sequence = "".join(Sequence.split()).upper()
        return Sequence[::-1]

    def test(self):
        inp = "AAAAGCGCCGCGCGAAATTTTCCCGCGCGCTTTTTCGCTCGCGCTCGCGCTCGCCGT"
        print(self.DNA_to_RNA(inp))

if __name__ == "__main__":
    a = DNA_Analysis()
    print("This is a test!")
    a.test()
    print("Test finished!")


# inp = "AAAAGCGCCGCGCGAAATTTTCCCGCGCGCTTTTTCGCTCGCGCTCGCGCTCGCCGT"
# a = DNA_Analysis()

# print(a.DNA_to_RNA(inp))
# print(a.RNA_to_Amino(a.DNA_to_RNA(inp2)))
# print(a.three_to_One(a.RNA_to_Amino(a.DNA_to_RNA(inp2), 3)))

# with open("sequence.txt", "r") as fin:
#     r = str(fin.read())
#     print(a.three_to_One(a.RNA_to_Amino(a.DNA_to_RNA(r), 1)))
#     print(a.find_CpG(r))
#     with open("islands.txt", "w") as fout:
#         fout.writelines(item + "\n" for item in a.find_CpG(r))

# print(a.get_complementary(inp))

# print(a.get_three_to_five(inp))

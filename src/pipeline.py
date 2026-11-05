from Bio import Align
from Bio.Seq import Seq

# Reference and sample DNA (starts with ATG, ends with a stop codon)
reference = Seq("ATGGCCATTGTAATGGGCCGCAAAGGGTGCCCGTAG")
sample    = Seq("ATGGCCATTGTAATGGGCCACAAAGGATGCCCGTAG")

# Step 1: align the two sequences
aligner = Align.PairwiseAligner()
aligner.mode = "global"
alignment = aligner.align(reference, sample)[0]
print("Alignment score:", alignment.score)


def classify(ref_aa, new_aa):
    if ref_aa == new_aa:
        return "synonymous"
    if new_aa == "*":
        return "nonsense"
    return "missense"


# Step 2: find mutations and annotate their effect on the protein
print("\nVariants found:")
for (r_start, r_end), (q_start, _) in zip(*alignment.aligned):
    for i in range(r_end - r_start):
        r_pos, q_pos = r_start + i, q_start + i
        ref_base, new_base = reference[r_pos], sample[q_pos]
        if ref_base != new_base:
            codon_start = (r_pos // 3) * 3
            ref_codon = reference[codon_start:codon_start + 3]
            new_codon = sample[codon_start:codon_start + 3]
            ref_aa = str(ref_codon.translate())
            new_aa = str(new_codon.translate())
            print(
                f"Position {r_pos + 1}: {ref_base}>{new_base} | "
                f"codon {codon_start // 3 + 1}: {ref_codon}>{new_codon} | "
                f"{ref_aa}>{new_aa} ({classify(ref_aa, new_aa)})"
            )

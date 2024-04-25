from guide import *
import numpy as np
import matplotlib.pyplot as plt

# s1 = Sequence('THISLINE-', Protein_Alphabet, gappy=True)
# s2 = Sequence('ISALIGNED', Protein_Alphabet)
# aln = Alignment([s1, s2])
# # print(aln)

b62 = readSubstMatrix("Files_WS2/blosum62.matrix", Protein_Alphabet)

seqs = readFastaFile("dao.fa", Protein_Alphabet)

my_student_number = 42066855
my_sequence_number = my_student_number % (len(seqs) - 1)
print('I will be using sequence number %d, which is %s' % (my_sequence_number, seqs[my_sequence_number].name))

def exercise_4(my_sequence_number, seqs, sub_matrix):
    """
    find two putative (assumed to be related) homologous (proteins sharing common ancestory)
    """

    # a bunch of temp vars to track the closest matches
    best1 = 0
    idx1 = 0
    best2 = 0
    idx2 = 0
    total = len(seqs)

    for idx in range(total):
        # Progress UI
        print(f"Running for {idx} / {total}")
        # do not match against itself - not a putative homolog
        if idx != my_sequence_number:
            seq = seqs[idx]
            # align my sequence and the current one and score
            aln = align(seqs[my_sequence_number], seq, sub_matrix, -7)
            percent = scoreAlignment(aln) / aln.alignlen

            # if it is the best, save it.
            # this means the current best is demoted the second best
            if best1 < percent:
                best2 = best1
                idx2 = idx1
                best1 = percent
                idx1 = idx

            # new second best, replace existing one
            elif best2 < percent:
                best2 = percent
                idx2 = idx

    print(f"Best match index is is {idx1}. Sequence is:\n\n{seqs[idx1]}")
    print(f"\n\nSecond best match index is is {idx2}. Sequence is:\n\n{seqs[idx2]}")
    with open("exercise_4.txt", "w") as f:
        f.write(f"{seqs[idx1].name} | percent: {best1}")
        f.write("\n")
        f.write(f"{seqs[idx2].name} | percent: {best2}")

def getConsensusForColumn(aln, colidx) -> str:
    symcnt = {}
    for seq in aln.seqs:
        mysym = seq[colidx]
        if mysym in symcnt:
            symcnt[mysym] += 1
        else:
            symcnt[mysym] = 1
    consensus = ""
    maxcnt = 0
    for mysym in symcnt:
        if symcnt[mysym] > maxcnt:
            maxcnt = symcnt[mysym]
            consensus = mysym
    return consensus


def exercise_5(aln: Alignment) -> Sequence:
    cols = len(aln[0])
    con_seq = ""
    # loop each col and grab the next char
    for i in range(cols):
        char = getConsensusForColumn(aln, i)
        # append it
        con_seq += char
    # return a Sequence instance
    return Sequence(con_seq, Protein_Alphabet, name='Consensus Sequence', gappy=True)

"""
Mine is index 83. sp|B1JGH2|MNMC_YERPY
Best match index is is 38.
sp|A4TM73|MNMC_YERPP | percent: 1.0

Second best match index is is 74.
sp|A9R7V8|MNMC_YERPG | percent: 1.0
"""
# exercise_4(my_sequence_number, seqs, b62)

aln = readClustalFile('dao.aln', Protein_Alphabet)
print('Loaded %d sequences into the alignment, which is %d columns wide' % (len(aln), aln.alignlen))
seq = exercise_5(aln)
print(f"Consensus Sequence for Exercise 5 is:\n\n{seq}")
# print(seq)

# print(len(aln[0]))

my_favourites = [
    "sp|Q4QKP6|MNMC_HAEI8",
    "sp|A5UCC9|MNMC_HAEIE",
    "sp|A5UEH1|MNMC_HAEIG",
    "sp|P44246|MNMC_HAEIN"
]

selected_seqs = []
for seq in aln.seqs:
    for name in my_favourites:
        if seq.name == name:
            selected_seqs.append(seq)

selected = Alignment(selected_seqs)


def exercise_6():
    fractional = selected.calcDistances('fractional')
    poisson = selected.calcDistances('poisson')
    gamma = selected.calcDistances('gamma')
    fig, ax = plt.subplots()

    print(fractional)
    print(poisson)
    print(gamma)
    ax.imshow(poisson, plt.cm.gray, interpolation='nearest')
    plt.yticks(np.arange(len(selected)), [s.name for s in selected])
    plt.title('Distances')
    plt.show()


# exercise_6()
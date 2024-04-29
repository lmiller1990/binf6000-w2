from guide import *

"""
Some useful utils
"""


def readfile(path):
    with open(path, "r") as file:
        return file.readlines()


def exercise_1(query_seq):
    """
    Usage:
    exercise_1("TMM25_HUMAN")
    """

    def process(txt):
        """make it a bit cleaner"""
        data = []
        for line in txt:
            data.append([s.strip() for s in line.split(" ") if s])
        return data

    def getid(lines):
        for line in lines:
            if line[0] == "ID":
                return line[1]

    def get_aa_seq(lines, signal_start, signal_end):
        aa_start = 0
        aa_end = 0
        for lineno, line in enumerate(lines):
            if line[0] == "SQ":
                # everything from SQ to end is protein sequence
                aa_start = lineno
            if line[0] == "//":
                aa_end = lineno

        aa = "".join(
            [item for sublist in lines[aa_start + 1 : aa_end] for item in sublist]
        )
        # minus 1, lists start at 0
        return aa[signal_start - 1 : signal_end]

    def get_signal_range(lines: list[list[str]]) -> tuple[int, int]:
        for line in lines:
            if line[0] == "FT" and line[1] == "SIGNAL":
                r = line[2]
                [s, e] = r.split("..")
                return int(s), int(e)

        raise RuntimeError("Did not find signal range")

    stream = getSequence(query_seq, format="txt")
    txt = stream.splitlines()
    lines = process(txt)
    seqid = getid(lines)
    sig_start, sig_end = get_signal_range(lines)
    aa_seq = get_aa_seq(lines, sig_start, sig_end)

    print(f"Currently processing the sequence with the name UNIPROT:{seqid}")
    print(f"There is a signal peptide ending at  position {sig_end}")
    print(f"The signal peptide looks like this: {aa_seq}")


def exercise_2():
    """Usage:
    exercise_2()
    """
    # list of `Sequence` class with useful information and methods about a given sequence
    seqs = []
    # list of ids representing an organism that contains DAO and has been reviewed
    names = searchSequences("family:DAO AND reviewed:true")
    #  count of occurrences of E. Coli DOA protein
    cnt = 0
    # get each sequence data from uniprot
    for name in names:
        seq = getSequence(name)
        # add to list for writing to file later
        seqs.append(seq)

        # ensure start of sequence is found in data
        start_index = seq.annot.find("OS=")

        # it's found!
        if start_index != -1:
            end_index = seq.annot.find("=", start_index + 3)

            # ensure end of sequence is found in data
            if end_index == -1:
                end_index = len(seq.annot)
            # python slicing to get annotation data
            species = seq.annot[start_index + 3 : end_index]
            # print name of seq
            print(seq.name, "\t", species)

            # if it's E. Coli, increment  count
            if species.startswith("Escherichia coli"):
                cnt = cnt + 1
        else:
            print(seq.name, "\tno species")

    # write to file
    writeFastaFile("dao.fa", seqs)

    if cnt > 0:
        print("Found %d exemplars in E. coli" % cnt)


def getConsensusForColumn(aln, colidx):
    # map of count of each symbol
    symcnt = {}
    for seq in aln.seqs:
        mysym = seq[colidx]
        if mysym in symcnt:
            # increment if already found
            symcnt[mysym] += 1
        else:
            # add to map if first time encountering
            symcnt[mysym] = 1
    consensus = ""
    maxcnt = 0
    # check each symbol and select the most frequently occurring
    for mysym in symcnt:
        if symcnt[mysym] > maxcnt:
            maxcnt = symcnt[mysym]
            consensus = mysym
    # return consensus char for this column
    return consensus


def exercise_5(aln):
    # get length of the sequence
    cols = len(aln[0])
    con_seq = ""
    # loop each col and grab the next char
    for i in range(cols):
        char = getConsensusForColumn(aln, i)
        # append it
        con_seq += char
    # return a Sequence instance
    return Sequence(con_seq, Protein_Alphabet, name="Consensus Sequence", gappy=True)


import pandas as pd


def get_yeasts():
    """Get a list of all the yeast species"""
    df = pd.read_csv("Files_WS2/sugars.csv")
    return list(df["Yeast"])


def exercise_7():
    """
    Curate the MalS.fa dataset
    """
    # Get all yeast from sugars.csv
    yeasts = get_yeasts()
    mals = readFastaFile("Files_WS2/MalS.fa", Protein_Alphabet)
    cnt = 0
    select = []
    # traverse mals dataset
    for seq in mals:
        # we are only interested in yeast also in the sugars.csv list
        if seq.name in yeasts:
            # white space -> _
            new_name = (seq.name + "_" + seq.annot).replace(" ", "_")
            seq.name = new_name
            cnt += 1
            select.append(seq)
        writeFastaFile("select.fa", select)
        print("Selected", cnt, "sequences")


voordeckers = [173, 231, 232, 233, 234, 294, 295, 324, 437]


def exercise_9():
    """
    This code print: YVGSLMQDE. This aligns with the figure 4 in the paper.
    """
    seq = "MVSMPN ... omitted for brevity ... NESA--"
    for a in voordeckers:
        print(seq[a - 1], end="")


def exercise_10():
    aln = readClustalFile("exercise7.aln", Protein_Alphabet)
    phy = runUPGMA(aln, "poisson")
    print(phy)  # print so we can copy into jstree
    return phy


def exercise_11():
    """Annotate tree with metabolized sugars"""
    df = pd.read_csv("Files_WS2/sugars.csv")
    yeasts = list(df["Yeast"])
    aln = readClustalFile("exercise7.aln", Protein_Alphabet)
    phy = runUPGMA(aln, "poisson")
    df = df.set_index(df["Yeast"])

    for col_name in df.columns[1:]:
        phy_copy = str(phy)
        t = 0
        f = 0
        for yeast in yeasts:
            can_metabolize = df[col_name][yeast]
            # add the labels using string replacement
            phy_copy = phy_copy.replace(yeast, f"{yeast} {can_metabolize}")
            # debugging - can see how many sugars can be metabolized
            if can_metabolize:
                t += 1
            else:
                f += 1

        # Prints: <Melizitose> => 8 / 14 = 57.14285714285714
        # Useful for debugging
        print(f"{col_name} => {t} / {f + t} = {(t / (t + f)) * 100}")


def exercise_12():
    residues = [173, 230, 231, 232, 233, 290, 291, 326, 437]

    aln = readClustalFile("exercise7.aln", Protein_wGAP)
    tree = runUPGMA(aln, "poisson")
    tree.putAlignment(aln)
    tree.parsimony()
    s = tree.strSites(residues)
    print(s)

from guide import *


def readfile(path):
    with open(path, "r") as file:
        return file.readlines()


def exercise_1(query_seq: str):
    def process(txt):
        """make it a bit cleaner"""
        # return [[s.strip() for s in line.split(" ") if s] for line in txt]
        data = []
        for line in txt:
            data.append([s.strip() for s in line.split(" ") if s])
        return data

    def getid(lines):
        for line in lines:
            if line[0] == "ID":
                return line[1]

    def get_aa_seq(lines: list[list[str]], signal_start: int, signal_end: int):
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


def exercise_3():
    s = """
    The fasta file format contains genome sequence data.
    Format is a '>' followed by a description of the sequence (on the same line).
    Following line is the sequence.
    Multiple sequences can be included in a single file. Empty lines are ignored.

    The example provided will fail because the second sequence has leading white space
    before the > character.
    """

    print(s)


# Execute all the exercises!
# exercise_1("TMM25_HUMAN")

# exercise_2()

exercise_3()

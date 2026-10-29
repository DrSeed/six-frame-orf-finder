# Six-frame open reading frame finder using the standard genetic code.
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def _orfs_one_strand(seq, min_aa, starts):
    found = []
    for frame in range(3):
        start = None
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if start is None and codon in starts:
                start = i
            elif start is not None and CODON_TABLE.get(codon) == "*":
                aa_len = (i - start) // 3
                if aa_len >= min_aa:
                    found.append((frame, start, i + 3, aa_len))
                start = None
    return found


def find_orfs(seq, min_aa=100, starts=("ATG",)):
    """Return ORFs (ATG..stop) on both strands, with coordinates on the + strand."""
    seq = seq.upper()
    n = len(seq)
    orfs = []
    for frame, s, e, aa in _orfs_one_strand(seq, min_aa, set(starts)):
        orfs.append({"strand": "+", "frame": frame + 1, "start": s, "end": e, "length_aa": aa,
                     "protein": translate(seq[s:e - 3])})
    rc = revcomp(seq)
    for frame, s, e, aa in _orfs_one_strand(rc, min_aa, set(starts)):
        orfs.append({"strand": "-", "frame": -(frame + 1), "start": n - e, "end": n - s, "length_aa": aa,
                     "protein": translate(rc[s:e - 3])})
    return sorted(orfs, key=lambda o: o["start"])

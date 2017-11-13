"""pdt.py.

Haystack
Primer finding (needle) in contigs (haystack)

Example:
    Haystack will discover primers with specific settings in a given fasta file
    that contains one or more contigs. Example usage is shown below:
        $ python haystack.py file.fasta # will display result in line
        $ python haystack.py file.fasta > test.out # will display in file

        Analyzing: ./BB3-1_Anello_contigs_labeled.fasta
        Contig: BB3-1_cap3-contigs.fa19
        [{'combo': (250, 269, 288),
        'ga_f': 'CCGAAGGTGAGTGAAACCACCGCGAAGAGCTAGAGGATCCCCG',
        'ga_r': 'CCGAATTGCCCCTTGACTACGGCGAAGAGCTTTTGAGTCGACCTG',
        'primer_f': 'CCGTAGTCAAGGGGCAATTCGG',
        'primer_r': 'CGGTGGTTTCACTCACCTTCGG',
        'primers_info': [[57.65, 0.0, -2.3], [57.8, -0.8, -2.77]],
        'virus_hed_dg': -2.16}]

        $

--------------------------------------------------------------------------------

Attributes:
    MAX_CONTIGS (int): Max number of contigs to print for each fasta file.
        Assumes that the fasta files are arranged from strongest read on the
        top of the document. Default is 3.

    MIN_LEN (int): Minimum length of primers in bp. Default is 20.

    MAX_LEN (int): Maximum length of primers in bp. Default is 27.

    CONTIG_BOUNDARY (int): Padding in evaluating possible primers.
        Padding applies to both beginning and end of the contig files.
        Default is 50.

    GC_PADDING (int): Number of S nucleotides (either C or G) that are at the
        beginning or end of the primer.

    MIN_TM (float): Mininum melting temperature of primer in deg C.
        Default is 56.

    MIN_TM (float): Maximum melting temperature of primer in deg C.
        Default is 58.

    MIN_HP_DG (float): Minimum delta G of hairpin formation in kcal/mol.
        Default is -9.

    MIN_HD_DG (float): Minimum delta G of homodimer formation in kcal/mol.
        Default is -9.

    MIN_HED_DG (float): Minimum delta G of heterodimer formation in kcal/mol.
        Default is -9.

Todo:
    * Hookup argument parsing

"""

import argparse
import itertools
import pprint as pp
import primer3
from Bio import SeqIO


MAX_CONTIGS = 3

MIN_LEN = 20
MAX_LEN = 27
CONTIG_BOUNDARY = 50
GC_PADDING = 3
MIN_TM = 56
MAX_TM = 58

MIN_HP_DG = -9  # hairpin dG
MIN_HD_DG = -9  # homodimer dG
MIN_HED_DG = -9 # heterodimer dG

VECTOR_F = 'CGAAGAGCTAGAGGATCCCCG'
VECTOR_R = 'CGAAGAGCTTTTGAGTCGACCTG'


def open_file(file_name):
    with open(file_name, "r") as fasta_file:
        seqs = list(SeqIO.parse(fasta_file, "fasta"))[0:MAX_CONTIGS]

        return seqs

def seq_walk(i, idxs, valid, valid_targets, rev=False):
    if i == len(idxs):
        return valid

    direction = -1 if rev else 1
    for mark in range(i, 0 if rev else len(idxs), direction):
        distance = abs(idxs[i] - idxs[mark]) + GC_PADDING - 1

        if distance > MAX_LEN:
            break
        elif distance >= MIN_LEN:
            try:
                valid_targets.add(idxs[mark])
                valid[idxs[i]].append(idxs[mark])
            except KeyError:
                valid[idxs[i]] = [idxs[mark]]

            distance = abs(idxs[i] - idxs[mark]) + GC_PADDING - 1

    return seq_walk(i + 1, idxs, valid, valid_targets)


def seq_walker(idxs):
    valid_fwd_idxs = {}
    valid_fwd_targets = set()

    seq_walk(0, idxs, valid_fwd_idxs, valid_fwd_targets)

    valid_middle_idxs = [x for x in valid_fwd_targets if x in valid_fwd_idxs]
    valid_back_idxs = [x for x in valid_middle_idxs if x in valid_fwd_idxs]

    all_combos = []
    for back in valid_back_idxs:
        for mid in valid_fwd_idxs[back]:
            try:
                combo_setup = [[back], [mid], valid_fwd_idxs[mid]]
                combo_nested = list(itertools.product(*combo_setup))

                all_combos += combo_nested
            except KeyError:
                continue

    return all_combos


def find_GC(seq):
    wrong_nt = ["A", "T"]
    idxs = []

    for i in range(CONTIG_BOUNDARY, len(seq) - CONTIG_BOUNDARY, 1):
        triplet = seq[i:i+GC_PADDING]
        if not any(nt in triplet for nt in wrong_nt):
            idxs.append(i)

    return idxs


def calculate_primer_info(combo, seq):
    virus_f = seq[combo[1]:combo[2] + GC_PADDING]
    virus_r = seq[combo[0]:combo[1] + GC_PADDING].reverse_complement()

    # from TL
    ga_f = str(virus_r.reverse_complement()) + VECTOR_F
    ga_r = str(virus_f.reverse_complement()) + VECTOR_R

    primers = [str(virus_f), str(virus_r)]

    primers_info = []
    for primer in primers:
        tm = primer3.calcTm(primer)
        hp = primer3.calcHairpin(primer).dg / 1000.0
        hod_dg = primer3.calcHomodimer(primer).dg / 1000.0

        if tm < MIN_TM or tm > MAX_TM:
            return None
        elif hp < MIN_HP_DG:
            return None
        elif hod_dg < MIN_HD_DG:
            return None

        primers_info.append([round(x, 2) for x in [tm, hp, hod_dg]])

    virus_hed_dg = primer3.calcHeterodimer(primers[0], primers[1]).dg / 1000.0

    if virus_hed_dg < MIN_HD_DG:
        return None

    final_info = {
        "combo": combo,
        "primer_f": primers[0],
        "primer_r": primers[1],
        "primers_info": primers_info,
        "virus_hed_dg": round(virus_hed_dg, 2),
        "ga_f": ga_f,
        "ga_r": ga_r
    }

    return final_info


def main():
    parser = argparse.ArgumentParser(description='Design primers from fasta contigs.')

    parser.add_argument('file', nargs='*', help='fasta file to read')
    parser.add_argument('--minTM', type=float,
        help='minimum Tm (deg C)', default=MIN_TM)
    parser.add_argument('--maxTM', type=float,
        help='maximum Tm (deg C)', default=MAX_TM)
    parser.add_argument('--hp', type=float,
        help='minimum hairpin delta G (kcal/mol)', default=MIN_HP_DG)
    parser.add_argument('--hoD', type=float,
        help='minimum homodimer delta G (kcal/mol)', default=MIN_HD_DG)
    parser.add_argument('--heD', type=float,
        help='minimum heterodimer delta G (kcal/mol)', default=MIN_HED_DG)

    args = parser.parse_args()

    for file_name in args.file:
        try:
            print("Analyzing: %s" % file_name)
            seqs = open_file(file_name)
        except FileNotFoundError:
            print("Error: %s not found" % (file_name))
            continue
        except IsADirectoryError:
            print("Error: cannot accept a directory. Please use a wildcard or input files only.")
            continue

        for i in range(0, MAX_CONTIGS):
            try:
                print("Contig: %s" % seqs[i].id)
                idxs = find_GC(seqs[i].seq.upper())

                all_combos = seq_walker(idxs)
                all_info = []
                for combo in all_combos:
                    info = calculate_primer_info(combo, seqs[i].seq)
                    if info:
                        all_info.append(info)

                all_info_sorted = sorted(all_info, key=lambda k: k['primers_info'][0], reverse=True)
                pp.pprint(all_info_sorted[:3])

                print("\n")

            except IndexError:
                break

        print("-" * 80)
main()

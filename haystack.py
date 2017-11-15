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

        $ [...]

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


def open_file(file_name, args):
    with open(file_name, "r") as fasta_file:
        seqs = list(SeqIO.parse(fasta_file, "fasta"))[0:args.maxContigs]

        return seqs

def seq_walk(i, idxs, valid, args, rev=False):
    if i == len(idxs):
        return valid

    direction = -1 if rev else 1
    upper_limit = 0 if rev else len(idxs)

    for mark in range(i, upper_limit, direction):
        distance = abs(idxs[i] - idxs[mark]) + args.gcPad - 1

        if distance > args.maxLen:
            break
        elif distance >= args.minLen:
            try:
                valid[idxs[i]].append(idxs[mark])
            except KeyError:
                valid[idxs[i]] = [idxs[mark]]

            distance = abs(idxs[i] - idxs[mark]) + args.gcPad - 1

    return seq_walk(i + 1, idxs, valid, args)


def seq_walker(idxs, args):
    valid_idxs = {}

    seq_walk(0, idxs, valid_idxs, args)

    all_combos = []
    for back in valid_idxs:
        for mid in valid_idxs[back]:
            try:
                combo_setup = [[back], [mid], valid_idxs[mid]]
                combo_nested = list(itertools.product(*combo_setup))

                all_combos += combo_nested
            except KeyError:
                continue

    return all_combos


def find_GC(seq, args):
    wrong_nt = ["A", "T"]
    idxs = []

    for i in range(args.contigBoundary, len(seq) - args.contigBoundary, 1):
        triplet = seq[i:i+GC_PADDING]
        if not any(nt in triplet for nt in wrong_nt):
            idxs.append(i)

    return idxs


def calculate_primer_info(combo, seq, args):
    virus_f = seq[combo[1]:combo[2] + args.gcPad]
    virus_r = seq[combo[0]:combo[1] + args.gcPad].reverse_complement()

    # from TL
    ga_f = str(virus_r.reverse_complement()) + VECTOR_F
    ga_r = str(virus_f.reverse_complement()) + VECTOR_R

    primers = [str(virus_f), str(virus_r)]

    primers_info = []
    for primer in primers:
        tm = primer3.calcTm(primer)
        hp = primer3.calcHairpin(primer).dg / 1000.0
        hod_dg = primer3.calcHomodimer(primer).dg / 1000.0

        if tm < args.minTM or tm > args.maxTM:
            return None
        elif hp < args.hp:
            return None
        elif hod_dg < args.hoD:
            return None

        primers_info.append([round(x, 2) for x in [tm, hp, hod_dg]])

    virus_hed_dg = primer3.calcHeterodimer(primers[0], primers[1]).dg / 1000.0

    if virus_hed_dg < args.heD:
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
    parser.add_argument('--minLen', type=int,
                        help='minimum primer length (bp)', default=MIN_LEN)
    parser.add_argument('--maxLen', type=int,
                        help='maximum primer length (bp)', default=MAX_LEN)

    parser.add_argument('--contigBoundary', type=int,
                        help='contig boundary (bp)', default=CONTIG_BOUNDARY)
    parser.add_argument('--maxContigs', type=int,
                        help='max number of contigs to analyze per fasta file',
                        default=MAX_CONTIGS)
    parser.add_argument('--gcPad', type=int,
                        help='number of G or C nucleotides at front and end of primer',
                        default=GC_PADDING)

    parser.add_argument('--minTM', type=float,
                        help='minimum Tm (deg C)', default=MIN_TM)
    parser.add_argument('--maxTM', type=float,
                        help='maximum Tm (deg C)', default=MAX_TM)

    parser.add_argument('--hp', type=float,
                        help='minimum hairpin delta G (kcal/mol)',
                        default=MIN_HP_DG)
    parser.add_argument('--hoD', type=float,
                        help='minimum homodimer delta G (kcal/mol)',
                        default=MIN_HD_DG)
    parser.add_argument('--heD', type=float,
                        help='minimum heterodimer delta G (kcal/mol)',
                        default=MIN_HED_DG)

    args = parser.parse_args()

    time_start = time.clock()

    for file_name in args.file:
        try:
            print("Analyzing: %s" % file_name)
            seqs = open_file(file_name, args)
        except FileNotFoundError:
            print("Error: %s not found" % (file_name))
            continue
        except IsADirectoryError:
            print("Error: cannot accept a directory. Please use a wildcard or input files only.")
            continue

        for i in range(0, args.maxContigs):
            try:
                print("Contig: %s" % seqs[i].id)
                idxs = find_GC(seqs[i].seq.upper(), args)

                all_combos = seq_walker(idxs, args)
                all_info = []
                for combo in all_combos:
                    info = calculate_primer_info(combo, seqs[i].seq, args)
                    if info:
                        all_info.append(info)

                all_info_sorted = sorted(all_info, key=lambda k: k['primers_info'][0], reverse=True)
                pp.pprint(all_info_sorted[:3])

                print()

            except IndexError:
                break
            except RecursionError:
                print("Error: Fasta file may not be formatted correctly.")
                break

        print("-" * 80)
main()

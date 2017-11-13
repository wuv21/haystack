# haystack
Primer finding (needle) in contigs (haystack)

# Example usage :
Haystack will discover primers with specific settings in a given fasta file that contains one or more contigs. Example usage is shown below:

```bash
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
```

Fasta files can be batch imported into the program with wildcard characters.

```bash
# will analyze all fasta files in current folder with surrounding 5 G or C nucleotides on either end of primer
$ python haystack.py --gcPad 5 *.fasta
```

## Arguments and Flags
\# *file*
> Required - Must be a fasta file or a collection of fasta files (by wildcard specification)

### Contig area to look
\# *--contigBoundary*
> Optional (int) - Contig boundary for haystack to start searching in bp. Default is 50bp.

\# *--maxContigs*
> Optional (int) - Maximum number of contigs to analyze per fasta file. Default is 3.

### Primer general settings
\# *--minLen*
> Optional (int) - Minimum primer length in bp. Default is 20bp.

\# *--maxLen*
> Optional (int) - Maximum primer length in bp. Default is 27bp.

\# *--gcPad*
> Optional (int) - Number of G or C nucleotides to have for beginning and end of. Default is 3.

\# *--minTm*
> Optional (float) - Minimum melting temperature in deg C. Default is 56 deg C.

\# *--maxTm*
> Optional (float) - Maximum melting temperature in deg C. Default is 58 deg C.

### Primer interaction settings
\# *--hp*
> Optional (float) - Minimum hairpin delta G in kcal/mol. Default is -9 kcal/mol.

\# *--hoD*
> Optional (float) - Minimum homodimer delta G in kcal/mol. Default is -9 kcal/mol.

\# *--heD*
> Optional (float) - Minimum heterodimer delta G in kcal/mol with forward and reverse primer. Default is -9 kcal/mol.
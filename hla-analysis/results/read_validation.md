# Targeted read validation checkpoint

DDBJ Slurm array 20604953 examined primary MHC-recruited reads after excluding secondary, supplementary, duplicate-marked and QC-failed reads. Six tasks completed; three failed because the local CRAM files for HG01252, HG01358 and HG01943 do not exist. No missing data are treated as negative evidence. Specificity job 20605048 completed.

Genomic 31-mer probes compare annotated assembly bases with the closest selected reference CDS. Alternative probes carry all alternative CDS SNPs within the window, but retain assembly-derived intronic context. Counts require all 31 bases at Phred ≥20 and count each read-name fragment once across mates. Probes are checked against both donor MHC assemblies. This is not a genome-wide mapping or allele-phasing assay. Reference-selection ties and alignments are retained in `read_marker_alignment_audit.tsv`.

## Coding candidates

- HLA-PROT-C-5948b39ec166038b: all tested distinguishing SNPs have at least five supporting fragments and a unique match in the donor MHC assemblies (NA18620). The full allele and novelty remain unconfirmed.
- HLA-PROT-DRB1-ad48d0a393823ed2: all tested distinguishing SNPs have at least five supporting fragments and a unique match in the donor MHC assemblies (HG02976). The full allele and novelty remain unconfirmed.
- HLA-PROT-DRB1-ba0192aed98e11e3: all tested distinguishing SNPs have at least five supporting fragments and a unique match in the donor MHC assemblies (NA19159). The full allele and novelty remain unconfirmed.

## Experimental disagreements

NA18943 HLA-A and DRB1 have identical coding sequences in the HPRC and JaSaPaGe assemblies; local distinguishing SNPs have 9–15 supporting fragments and unique donor-MHC matches. NA19007 HLA-A has 14 supporting fragments at its distinguishing SNP, with no alternative-probe fragments. This supports the assembly bases without proving why the historical experimental labels differ.

NA18608 DRB1 has supporting reads at the tested bases, but several probes occur in two or three places in the donor MHC assemblies. Its diploid allele discrepancy remains unresolved; these markers cannot assign all relevant bases to the intended DRB1 copy. Alternative counts can originate from the other haplotype or paralogs and are not automatically contradictory evidence.

Raw counts, coordinates, alleles, probe sequences and specificity locations are preserved. The ≥5-fragment rule is an exploratory ranking threshold; it has not been calibrated for sensitivity or specificity. Nearby markers can count the same fragments and are not independent confirmations.

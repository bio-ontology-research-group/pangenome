# Post-hackathon validation and association plan

UK Biobank individual WGS/WES and phenotypes are available to the project, but
access is not operational during the hackathon. No UKB analysis has run here.
Freeze candidates and this analysis specification before examining phenotypes.

## Candidate definitions and technical gate

Use `results/coding_candidates.tsv` protein sequence hashes and stable candidate
IDs. These are database-absent protein candidates, not established new alleles.
To make each genotypable, supplement the protein definition with phased genomic
sequence, coordinate/reference-build anchors, distinguishing nucleotide variants,
paralog-specific markers, and an explicit representation of reference and alternate
haplotypes. A protein can arise from different nucleotide haplotypes; retain that
mapping rather than assuming a single causal DNA change.

`results/within_type_contrasts.tsv` IDs identify comparisons conditional on matching
labels. They are not yet deployable variant IDs. Split each selected comparison
into explicit alternate sequence paths or a defined CN/arrangement event, and freeze
the resulting event definition. Record assembly provenance, independent donor and
family support, gaps/clipping, current IPD matches and literature classification.
The known HG03139 DRB5-absence result is a positive control. It is additional to
DRB1–DQ labels, but conventional typing including DRB5 already describes it.

Before UKB association testing, validate a short list using matched source reads,
alternative assemblies or informative pedigrees. Define locus-specific unique
markers and test them against all panel paralogs. For a deletion, require supported
junctions where accessible and calibrated depth; missing annotation alone fails
the gate. For apparent coding changes, inspect local read alignment, allele balance,
coverage, phasing and competing paralogous alignments. Keep unvalidated candidates
in a separate exploratory list.

## Genotyping feasibility

| Candidate | Preferred data and measurement | Main limitation |
|---|---|---|
| Coding allele | WGS with locus-aware competitive alignment/local assembly and discriminating variants; WES where targeted exons suffice | WES cannot prove complete-gene identity; paralogs and capture bias |
| Flanking sequence haplotype | WGS sequence/path assignment with unique anchors | WES usually lacks required sequence; no regulatory function established |
| DRB gene CN/arrangement | WGS depth calibrated to controls plus junction/path support | Whole-locus similarity or gene annotation is insufficient for CN |
| C4 forms | Dedicated WGS copy-number and sequence-marker assay with positive controls | C4 calls alone do not resolve the complete RCCX arrangement |
| Full RCCX structure | Explicit CYP21/TNX/STK19/C4 module definitions; targeted validation for complex cases | May remain unresolvable with ordinary short reads |

Fit/calibrate the genotyping procedure using phenotype-blind controls. Measure
call rate, uncertainty, batch effects, depth dependence and ancestry-specific
performance. Use matched WGS/WES and any orthogonal truth to estimate concordance;
do not use agreement between two methods on the same biased reads as independent
validation. Include carriers and noncarriers, and keep assay-development samples
out of its final technical evaluation. Apply any technical exclusion criteria
before phenotype testing.

## Association tests

1. Choose a small, biologically justified primary phenotype set for validated
   candidates, with frozen case definitions, exclusions and covariates. Treat any
   phenome-wide screen as secondary with its own multiplicity correction.
2. Test candidate dosage or copy number with a phenotype-appropriate model,
   accounting for relatedness, ancestry, age, sex and relevant sequencing/batch
   effects. Select the mixed-model or rare-variant implementation after the data
   structure is known. Inspect carrier counts and separation before using
   asymptotic tests; do not promise power without frequencies and case counts.
3. Compare a conventional-HLA baseline to baseline plus candidate. The baseline
   must include the available classical loci, including DRB3/4/5 where reliable,
   and suitable local MHC ancestry/haplotype covariates. Avoid a redundant full
   set of perfectly collinear allele indicators. Report conditional effect sizes
   and confidence intervals, and evaluate predictive gain on held-out individuals
   if making a prediction claim.
4. Repeat within comparable HLA backgrounds where enough carriers exist. Check
   LD and genotype uncertainty: an unidentifiable conditional effect is not
   evidence that the candidate has no biological effect. Use burden tests only
   for prespecified, defensible functional groupings, not arbitrary graph clusters.
5. Evaluate transport across sufficiently represented ancestry groups and seek
   replication outside discovery donors. Separate lack of genotypability, lack of
   power and lack of association. A UKB association alone does not establish
   molecular function or clinical utility.

Deliverables are a frozen genotyping specification, technical validation metrics,
carrier-frequency/power assessment, prespecified association specification, and
conditional association/replication results. Only the specification and candidate
prioritisation belong to the current hackathon analysis.

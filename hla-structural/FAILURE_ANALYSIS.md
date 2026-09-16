# Remaining structural-typing errors

The corrected model has 15 full-signature errors among 88 additional donors.

- non_C4_prototype_signature_only: 1
- C4_order_or_chromosome_distribution: 1
- C4_marginal_form_dosage: 10
- C4_joint_form_phase_with_matching_marginal_dosage: 3

2 errors involve a true structural signature absent from the family-excluded reference catalogue. The method still selected a unique optimum: a numerical optimum must not be called a confident genotype.

The sampled marker vocabulary supplies conserved total-C4 dosage markers but no qualifying markers for the joint C4AL/C4AS/C4BL/C4BS copy counts. Those joint forms combine A/B sequence identity with long/short status at separated positions; failure to find a universal 31-mer for a joint form does not mean the underlying component variants are unmeasurable. The independent joint-form marker baseline is therefore an identifiability diagnostic, not a fair substitute for a dedicated C4 assay.

Next assay development should explicitly include separate A/B diagnostic sites, long/short insertion junctions, and module-boundary evidence, then infer their phase. Preserve ambiguous arrangement sets where short reads do not distinguish them. Validate with independent long reads or targeted assays, including candidate structures absent from the reference. Freeze any revised assay and test a new untouched sample set; these 88 donors have now been inspected and are not a fresh validation set for further tuning.

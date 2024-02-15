Code for "Reducing mapping reference and lineage bias in Mycobacterium tuberculosis"

1. Simulations and variant calling:
> ART Illumina simulations: 1_art_sims.sh
> Mapping only pipeline: 2_mapping.sh
> Hybrid mapping and assembly pipeline: 3_assembly_mapping.sh

2. Analysis of bias
> VCF comparison: 4_compareVCfs.R
> Build null distribution: 5_null_dist.R
> Infer error model: 6_errorModel.R
> Calculate and plot bias: 7_mappingBias.R
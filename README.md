# Early seeding of future Richter transformation in chronic lymphocytic leukemia - Repository

This repository contains scripts, notebooks, and data accompanying the study "Early seeding of future Richter transformation in chronic lymphocytic leukemia" published in Nature Medicine (Nadeu, Royo, Massoni-Badosa, Playa-Albinyana, Garcia-Torre, et al. 2022).

## Abstract

Richter transformation (RT) is a paradigmatic evolution of chronic lymphocytic leukemia (CLL) into a very aggressive large B-cell lymphoma conferring a dismal prognosis. The mechanisms driving RT remain largely unknown. We have characterized the whole genome, epigenome and transcriptome, combined with single-cell DNA/RNA sequencing analyses and functional experiments, of 19 CLL developing RT. Studying 54 longitudinal samples covering up to 19 years of disease course, we uncovered minute subclones carrying genomic, immunogenetic, and transcriptomic features of RT-cells already at CLL diagnosis, which were dormant for up to 19 years before transformation. We also identified new driver alterations, discovered a novel mutational signature (SBS-RT), recognized an oxidative phosphorylation (OXPHOS)<sup>high</sup>-B-cell receptor signaling (BCR)<sup>low</sup> transcriptional axis in RT, and showed that OXPHOS inhibition reduces the proliferation of RT cells. These findings demonstrate the early seeding of subclones driving advanced stages of cancer evolution and uncover therapeutic targets for the, once expanded, lethal RT.


## Analyses

- **Mutational signatures:** Code to reproduce the extraction/assignment and fitting of mutational signatures as well as the characterization of the novel SBS-RT can be found [here](https://github.com/ferrannadeu/RichterTransformation/tree/main/MutationalSignatures) 
- **DNA methylation:** Code used to normalize DNA methylation data can be found [here](https://github.com/Duran-FerrerM/DNAmeth_arrays). Code to calculate the tumor cell content, CLL epitypes, and epiCMIT can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome).
- **H3K27ac and ATAC-seq:** Code to reproduce H3K27ac and ATAC-seq analyses can be found [here](https://github.com/ferrannadeu/RichterTransformation/tree/main/H3K27ac_ATAC-seq)
- **Bulk RNA-seq:** Code to reproduce bulk RNA-seq analyses can be found [here](https://github.com/ferrannadeu/RichterTransformation/tree/main/bulkRNA-seq)
- **Single-cell RNA-seq:** Code to reproduce scRNA-seq analyses can be found [here](https://github.com/massonix/richter_transformation).


## Contributors

Ferran Nadeu, Romina Royo, Ramon Massoni-Badosa, Beatriz Garcia-Torre, Martí Duran-Ferrer.

# Metagenomics-tools-atlas

A community-maintained inventory of tools and resources for **short-read shotgun metagenomics** — from raw reads to ecological inference.

---

## About

This list was assembled to give researchers a single reference point for the metagenomics software landscape. It covers the full analysis arc: QC, host removal, taxonomic profiling, functional annotation, assembly, binning, MAG quality assessment, strain-level analysis, statistical testing, and more.

**Scope:** The primary focus is Illumina short-read shotgun metagenomics. Long-read and hybrid tools are included where they complement short-read workflows. Amplicon (16S/ITS) tools are out of scope unless they serve a clear dual-use role.

**Disclaimer:** This inventory is a mix of human-curated additions and AI-assisted search, with human checking throughout. It is not exhaustive and reflects the state of the field at the time of compilation. If your tool is not here, that is an oversight — not a judgement. Please suggest it via an [Issue](../../issues/new). See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

Share freely. Use as needed.

---

## Table of Contents

- [Antimicrobial Resistance](#antimicrobial-resistance)
- [Assembly](#assembly)
- [Batch Correction](#batch-correction)
- [Benchmarking & Simulation](#benchmarking--simulation)
- [BGCs / Natural Products / AMPs](#bgcs--natural-products--amps)
- [Binning](#binning)
- [Causal Inference & Mediation](#causal-inference--mediation)
- [CAZyme & Carbohydrate Metabolism](#cazyme--carbohydrate-metabolism)
- [Community Typing / Mixture Models](#community-typing--mixture-models)
- [Compositional Modelling](#compositional-modelling)
- [Contamination Detection](#contamination-detection)
- [Coverage & Abundance](#coverage--abundance)
- [Data Access & Public Resources](#data-access--public-resources)
- [Data Containers & Infrastructure](#data-containers--infrastructure)
- [Data Retrieval / Public Metagenomes](#data-retrieval--public-metagenomes)
- [Data Wrangling](#data-wrangling)
- [Dietary Reconstruction](#dietary-reconstruction)
- [Differential Abundance](#differential-abundance)
- [Diversity Analyses](#diversity-analyses)
- [DNA / Protein Foundation Models](#dna--protein-foundation-models)
- [Eukaryotic / Fungal / Parasite Profiling](#eukaryotic--fungal--parasite-profiling)
- [Functional Profiling (Read-Based)](#functional-profiling-read-based)
- [Fungal Functional Ecology](#fungal-functional-ecology)
- [Gene Prediction & MAG Annotation](#gene-prediction--mag-annotation)
- [Gene Prediction & Viral Annotation](#gene-prediction--viral-annotation)
- [Generative Sequence / Protein Design](#generative-sequence--protein-design)
- [Genome Size Estimation](#genome-size-estimation)
- [Genome-Scale Metabolic Modelling](#genome-scale-metabolic-modelling)
- [Hierarchical Testing](#hierarchical-testing)
- [Host Read Removal](#host-read-removal)
- [In Situ Growth Rate / Replication Rate](#in-situ-growth-rate--replication-rate)
- [Interactive Visualisation](#interactive-visualisation)
- [K-mer / Reference-Free Comparison](#k-mer--reference-free-comparison)
- [Longitudinal & Time-Series Analysis](#longitudinal--time-series-analysis)
- [Machine Learning & Prediction](#machine-learning--prediction)
- [MAG Dereplication & Genome Catalogues](#mag-dereplication--genome-catalogues)
- [MAG Quality Assessment](#mag-quality-assessment)
- [Meta-Analysis & Batch Adjustment](#meta-analysis--batch-adjustment)
- [Metabolite / Metagenome Inference](#metabolite--metagenome-inference)
- [Metagenomics Workflows & Pipelines](#metagenomics-workflows--pipelines)
- [Microbial GWAS & Pangenome Association](#microbial-gwas--pangenome-association)
- [Microbial Load](#microbial-load)
- [Microbiome Biological-Age / Clock Models](#microbiome-biological-age--clock-models)
- [Mobilome / Plasmids / MGEs](#mobilome--plasmids--mges)
- [Multi-Omics Integration](#multi-omics-integration)
- [Network & Co-occurrence Inference](#network--co-occurrence-inference)
- [Phage / Virus Discovery](#phage--virus-discovery)
- [Phage Host Prediction](#phage-host-prediction)
- [Phage Therapy & Viral Workflows](#phage-therapy--viral-workflows)
- [Phylogenetics & Evolution](#phylogenetics--evolution)
- [Prophage Detection](#prophage-detection)
- [R / Bioconductor Analysis Frameworks](#r--bioconductor-analysis-frameworks)
- [Read QC & Preprocessing](#read-qc--preprocessing)
- [Reference Databases — AMR](#reference-databases--amr)
- [Reference Databases — BGCs](#reference-databases--bgcs)
- [Reference Databases — Curated Cohorts](#reference-databases--curated-cohorts)
- [Reference Databases — Function](#reference-databases--function)
- [Reference Databases — MAG Catalogues](#reference-databases--mag-catalogues)
- [Reference Databases — Phages / Viruses](#reference-databases--phages--viruses)
- [Reference Databases — Taxa](#reference-databases--taxa)
- [Reference Databases — Other](#reference-databases--other)
- [Sequencing Depth & Coverage](#sequencing-depth--coverage)
- [Simulation & Benchmarking](#simulation--benchmarking)
- [Source Tracking](#source-tracking)
- [Strain-Level Analysis](#strain-level-analysis)
- [Structural Metagenomics / Protein Structure at Scale](#structural-metagenomics--protein-structure-at-scale)
- [Taxon Set / Signature Enrichment](#taxon-set--signature-enrichment)
- [Taxonomic Placement & MAG Taxonomy](#taxonomic-placement--mag-taxonomy)
- [Taxonomic Profiling (Read-Based)](#taxonomic-profiling-read-based)
- [Temporal & Dynamic Modelling](#temporal--dynamic-modelling)
- [Viral & Virome Analysis](#viral--virome-analysis)
- [Visualisation](#visualisation)
- [Workflow & Pipeline Frameworks](#workflow--pipeline-frameworks)

---

## Antimicrobial Resistance

| Tool | Link | Description |
|---|---|---|
| ResFinder | [link](https://bitbucket.org/genomicepidemiology/resfinder) | DTU AMR gene detection with phenotype prediction (ResFinder4). |
| RGI / CARD | [link](https://github.com/arpcard/rgi) | Resistance Gene Identifier searches CARD's curated AMR ontology and homologs. |
| deepARG | [link](https://github.com/gaarangoa/deeparg) | Deep-learning ARG classifier for divergent resistance genes. |
| AMR++ | [link](https://github.com/Microbial-Ecology-Group/AMRplusplus) | Pipeline for profiling antimicrobial resistance genes in metagenomic sequencing data. |
| AMRFinderPlus | [link](https://github.com/ncbi/amr) | NCBI-maintained AMR + virulence gene detection from genomes and proteins. |
| hAMRonization | [link](https://github.com/pha4ge/hAMRonization) | Harmonises outputs from multiple AMR detection tools into a common schema. |
| Staramr | [link](https://github.com/phac-nml/staramr) | Scans bacterial genomes for AMR genes, point mutations and plasmid replicons. |
| ABRicate | [link](https://github.com/tseemann/abricate) | Mass screening of contigs against multiple AMR DBs (CARD, ResFinder, NCBI, etc.). |
| ARGs-OAP | [link](https://github.com/xinehc/args_oap) | Online and offline ARG profiling from metagenomic reads. |
| ARIBA | [link](https://github.com/sanger-pathogens/ariba) | Antimicrobial Resistance Identification By Assembly. |

## Assembly

| Tool | Link | Description |
|---|---|---|
| MetaQUAST | [link](https://github.com/ablab/quast) | Quality assessment of metagenomic assemblies. |
| metaSPAdes | [link](https://github.com/ablab/spades) | SPAdes mode for metagenomic assembly with repeat resolution and bubble handling. |
| OPERA-MS | [link](https://github.com/CSB5/OPERA-MS) | Hybrid metagenomic assembler for short and long reads. |
| metaFlye | [link](https://github.com/fenderglass/Flye) | Flye assembler in metagenomic mode for long-read (ONT/PacBio) data. |
| metaMDBG | [link](https://github.com/GaetanBenoitDev/metaMDBG) | Long-read metagenome assembler using minimizer de Bruijn graphs. |
| IDBA-UD | [link](https://github.com/loneknightpy/idba) | Iterative de Bruijn assembler for metagenomic data with uneven depth. |
| PenguiN | [link](https://github.com/soedinglab/plass) | Protein-guided nucleotide assembler emphasising recovery of coding regions and viral sequences. |
| MEGAHIT | [link](https://github.com/voutcn/megahit) | Memory-efficient de Bruijn graph assembler for large metagenomes. |
| hifiasm-meta | [link](https://github.com/xfengnefx/hifiasm-meta) | HiFi-specific haplotype-resolved metagenomic assembler. |

## Batch Correction

| Tool | Link | Description |
|---|---|---|
| PLSDAbatch | [link](https://bioconductor.org/packages/PLSDAbatch) | Multivariate PLS-DA based framework for batch-effect correction in microbiome data. |
| ConQuR | [link](https://github.com/wdl2459/ConQuR) | Batch correction for microbiome count data using quantile regression. |

## Benchmarking & Simulation

| Tool | Link | Description |
|---|---|---|
| CAMISIM | [link](https://github.com/CAMI-challenge/CAMISIM) | Simulation framework for metagenomic benchmarking with realistic communities. |
| InSilicoSeq | [link](https://github.com/HadrienG/InSilicoSeq) | Illumina metagenomic read simulator with empirical error models. |

## BGCs / Natural Products / AMPs

| Tool | Link | Description |
|---|---|---|
| AMPSphere | [link](https://ampsphere.big-data-biology.org/) | Atlas of antimicrobial peptides from global metagenomes. |
| ARTS | [link](https://arts.ziemertlab.com/) | Antibiotic Resistant Target Seeker finds resistance-linked BGCs. |
| BGC Atlas | [link](https://bgc-atlas.cs.uni-duesseldorf.de/) | Atlas of BGCs from MAGs across diverse environments. |
| antiSMASH | [link](https://github.com/antismash/antismash) | Gold-standard rule-based BGC identification and annotation. |
| BiG-SLiCE | [link](https://github.com/medema-group/bigslice) | Clustering of millions of BGCs into Gene Cluster Families (GCFs) at scale. |
| deepBGC | [link](https://github.com/Merck/deepbgc) | Deep-learning BGC prediction via Pfam embeddings and LSTM classifier. |
| GECCO | [link](https://github.com/zellerlab/GECCO) | Gene Cluster Prediction with Conditional Random Fields. |

## Binning

| Tool | Link | Description |
|---|---|---|
| MetaBAT2 | [link](https://bitbucket.org/berkeleylab/metabat) | Tetranucleotide-frequency + coverage clustering for MAG recovery. |
| magpurify2 | [link](https://github.com/apcamargo/magpurify2) | Decontamination / purification of MAGs using multi-feature anomaly detection. |
| SemiBin2 | [link](https://github.com/BigDataBiology/SemiBin) | Self-supervised deep contrastive learning binner with pretrained models. |
| CONCOCT | [link](https://github.com/BinPro/CONCOCT) | Composition + coverage Gaussian mixture model binner. |
| DAS_Tool | [link](https://github.com/cmks/DAS_Tool) | Consensus refinement of multiple binner outputs by SCG-based scoring. |
| MetaDecoder | [link](https://github.com/liu-congcong/MetaDecoder) | Reference-free metagenomic binning using sequence composition and coverage. |
| GraphMB | [link](https://github.com/MicrobialDarkMatter/GraphMB) | Graph neural network-assisted metagenomic binning using assembly graph and coverage. |
| VAMB | [link](https://github.com/RasmussenLab/vamb) | Variational autoencoder-based binner using composition and co-abundance. |
| COMEBin | [link](https://github.com/ziyewang/COMEBin) | Contrastive multi-view representation learning for metagenomic binning. |
| MaxBin2 | [link](https://sourceforge.net/projects/maxbin2/) | Expectation-maximisation binning using marker genes + tetranucleotide frequency. |

## Causal Inference & Mediation

| Tool | Link | Description |
|---|---|---|
| SparseMCMM | [link](https://github.com/chanw0/SparseMCMM) | Sparse compositional mediation modelling for microbiome mediators. |
| microbiome_mediation (medTest) | [link](https://github.com/KiRinHong/miMediation) | Compositional mediation testing via distance-based and SOBEL-style frameworks. |
| MicroBVS | [link](https://github.com/mkoslovsky/MicroBVS) | Bayesian variable selection for microbiome covariate models. |
| TwoSampleMR / MR-Base | [link](https://github.com/MRCIEU/TwoSampleMR) | Mendelian randomisation using GWAS summary stats; microbiome instruments via mbQTL. |

## CAZyme & Carbohydrate Metabolism

| Tool | Link | Description |
|---|---|---|
| dbCAN | [link](https://bcb.unl.edu/dbCAN2/) | Automated CAZyme annotation using the CAZy ecosystem. |

## Community Typing / Mixture Models

| Tool | Link | Description |
|---|---|---|
| DirichletMultinomial | [link](https://bioconductor.org/packages/DirichletMultinomial) | Dirichlet-multinomial mixture modelling for microbial community types. |

## Compositional Modelling

| Tool | Link | Description |
|---|---|---|
| fido | [link](https://jsilve24.github.io/fido/) | Bayesian multinomial logistic-normal modelling for microbiome count/compositional data. |

## Contamination Detection

| Tool | Link | Description |
|---|---|---|
| decontam | [link](https://github.com/benjjneb/decontam) | Statistical identification of contaminant features using prevalence and/or frequency patterns. |

## Coverage & Abundance

| Tool | Link | Description |
|---|---|---|
| Koverage | [link](https://github.com/beardymcjohnface/Koverage) | Snakemake/Snaketool workflow for read-coverage analysis across massive metagenomic or genomic datasets. |

## Data Access & Public Resources

| Tool | Link | Description |
|---|---|---|
| HoloFoodR | [link](https://bioconductor.org/packages/HoloFoodR) | R/Bioconductor interface to HoloFood resource data. |
| MGnifyR | [link](https://bioconductor.org/packages/MGnifyR) | R/Bioconductor interface for accessing and processing MGnify data. |
| microbiomeDataSets | [link](https://bioconductor.org/packages/microbiomeDataSets) | Bioconductor ExperimentHub package providing microbiome datasets in analysis-ready formats. |

## Data Containers & Infrastructure

| Tool | Link | Description |
|---|---|---|
| SummarizedExperiment | [link](https://bioconductor.org/packages/SummarizedExperiment) | Generic Bioconductor container for rectangular assay data with row and column metadata. |
| TreeSummarizedExperiment | [link](https://bioconductor.org/packages/TreeSummarizedExperiment) | Bioconductor S4 container for assays with hierarchical row/column structure such as taxonomy trees and phylogenies. |

## Data Retrieval / Public Metagenomes

| Tool | Link | Description |
|---|---|---|
| ENA Browser / enaBrowserTools | [link](https://github.com/enasequence/enaBrowserTools) | Searches/downloads reads and metadata from ENA. |
| SRA Toolkit | [link](https://github.com/ncbi/sra-tools) | Downloads and converts sequencing data from NCBI SRA. |

## Data Wrangling

| Tool | Link | Description |
|---|---|---|
| tidySummarizedExperiment | [link](https://bioconductor.org/packages/tidySummarizedExperiment) | Tidyverse-style manipulation interface for SummarizedExperiment objects. |

## Dietary Reconstruction

| Tool | Link | Description |
|---|---|---|
| MEDI | [link](https://github.com/Gibbons-Lab/medi) | Metagenomic Estimation of Dietary Intake; estimates food-derived DNA and nutrient composition from stool shotgun metagenomes. |

## Differential Abundance

| Tool | Link | Description |
|---|---|---|
| ANCOMBC / ANCOM-BC2 | [link](https://bioconductor.org/packages/ANCOMBC) | Bias-corrected compositional differential abundance testing, including ANCOM-BC and ANCOM-BC2 implementations. |
| dar | [link](https://bioconductor.org/packages/dar) | Bioconductor package for differential abundance testing in microbiome data. |
| metagenomeSeq | [link](https://bioconductor.org/packages/metagenomeSeq) | Zero-inflated Gaussian/cumulative-sum scaling method for sparse microbiome count data. |
| MaAsLin3 | [link](https://github.com/biobakery/maaslin3) | Multivariate association with linear models joint presence + abundance modelling. |
| MaAsLin2 | [link](https://huttenhower.sph.harvard.edu/maaslin/) | Multivariable association testing for microbiome features and metadata. |
| Songbird | [link](https://github.com/biocore/songbird) | Multinomial regression for log-ratio differentials. |
| ALDEx2 | [link](https://github.com/ggloor/ALDEx_bioc) | CLR transform + Dirichlet Monte Carlo for compositional differential abundance. |
| ZicoSeq | [link](https://github.com/jchen1981/GUniFrac) | Permutation-based DA with reference taxa and covariate adjustment. |
| corncob | [link](https://github.com/statdivlab/corncob) | Beta-binomial regression modelling abundance and dispersion separately. |
| radEmu | [link](https://github.com/statdivlab/radEmu) | Method for estimating fold changes from partially observed outcomes with applications in microbial metagenomics. |
| ZINQ | [link](https://github.com/wdl2459/ZINQ-L) | Zero-inflated quantile rank-sum test for DA. |
| LinDA | [link](https://github.com/zhouhj1994/LinDA) | Linear model on log-CLR with bias correction (intercept median). |
| fastANCOM | [link](https://github.com/ZRChao/fastANCOM) | Speed-optimised reimplementation of ANCOM. |
| lefser | [link](https://bioconductor.org/packages/lefser) | R/Bioconductor implementation of LEfSe-style metagenomic biomarker discovery. |
| MicrobiomeStat | [link](https://cran.r-project.org/package=MicrobiomeStat) | CRAN package containing microbiome statistical methods including LinDA workflows. |
| expam | [link](https://github.com/seansolari/expam) | Distance estimation using trees. |

## Diversity Analyses

| Tool | Link | Description |
|---|---|---|
| iNEXT | [link](https://chao.shinyapps.io/iNEXTOnline/) | Rarefaction/extrapolation of diversity with Hill numbers. |
| vegan | [link](https://cran.r-project.org/package=vegan) | R package for classical community ecology: ordinations, permutation tests, diversity indices. |
| breakaway | [link](https://github.com/adw96/breakaway) | Estimates species richness from frequency count tables. |
| phyloseq | [link](https://github.com/joey711/phyloseq) | R package for storing and analysing OTU/ASV tables with taxonomy, sample data, tree. |
| philr | [link](https://github.com/jsilve24/philr) | Phylogenetic isometric log-ratio transformation for microbiome compositional data. |
| scikit-bio | [link](https://github.com/scikit-bio/scikit-bio) | Python toolkit for diversity metrics, phylogenetics, sequence ops. |
| QIIME2 | [link](https://qiime2.org/) | End-to-end microbiome ecology platform: diversity, taxonomy, stats, plugins. |
| MicrobiomeAnalyst | [link](https://www.microbiomeanalyst.ca/) | Web platform for marker-gene + shotgun analyses with stats, visualisation, signatures. |

## DNA / Protein Foundation Models

| Tool | Link | Description |
|---|---|---|
| ProtTrans | [link](https://github.com/agemagician/ProtTrans) | Family of transformer protein LMs (BERT/T5/XLNet) trained on UniRef. |
| Evo 2 | [link](https://github.com/ArcInstitute/evo2) | 40B-parameter genome-scale LM trained across the tree of life (1M token context). |
| Evo | [link](https://github.com/evo-design/evo) | 7B-parameter genomic language model trained on bacterial/phage DNA. |
| ESM-2 | [link](https://github.com/facebookresearch/esm) | Protein language model family (8M–15B params) trained on UniRef. |
| HyenaDNA | [link](https://github.com/HazyResearch/hyena-dna) | Long-context DNA LM using Hyena operators (up to 1M tokens). |
| Nucleotide Transformer | [link](https://github.com/instadeepai/nucleotide-transformer) | InstaDeep transformer family pretrained on human and 850 species genomes. |
| DNABERT-2 | [link](https://github.com/MAGICS-LAB/DNABERT_2) | BERT-based DNA model with BPE tokenisation and multi-species training. |
| GenSLM | [link](https://github.com/ramanathanlab/genslm) | Genome-scale language model trained for pathogen surveillance (e.g. SARS-CoV-2 evolution). |

## Eukaryotic / Fungal / Parasite Profiling

| Tool | Link | Description |
|---|---|---|
| EukDetect | [link](https://github.com/allind/EukDetect) | Detects eukaryotic organisms in shotgun metagenomic datasets using marker genes. |
| EukRep | [link](https://github.com/patrickwest/EukRep) | Classifies assembled contigs as eukaryotic vs prokaryotic. |
| FUNGuild | [link](https://github.com/UMNFuN/FUNGuild) | Assigns ecological guilds to fungal taxa. |

## Functional Profiling (Read-Based)

| Tool | Link | Description |
|---|---|---|
| DIAMOND | [link](https://github.com/bbuchfink/diamond) | Ultra-fast translated sequence alignment against protein databases. |
| HUMAnN3 | [link](https://github.com/biobakery/humann) | Tiered profiling of gene families (UniRef) and MetaCyc pathways from shotgun reads. |
| eggNOG-mapper | [link](https://github.com/eggnogdb/eggnog-mapper) | Functional annotation via orthology to eggNOG; assigns KEGG, GO, COG, CAZy. |
| SUPER-FOCUS | [link](https://github.com/metageni/SUPER-FOCUS) | SEED subsystem-level functional profiling using DIAMOND/RAPSearch reduced DB. |
| MMseqs2 | [link](https://github.com/soedinglab/MMseqs2) | Fast protein/nucleotide search, clustering, and taxonomy workflows. |
| ShortBRED | [link](https://huttenhower.sph.harvard.edu/shortbred/) | Builds and profiles short, specific protein family markers in metagenomes. |
| MinPath | [link](https://omics.informatics.indiana.edu/MinPath/) | Infers parsimonious biological pathways from gene-family annotations. |
| MEGAN6 | [link](https://software-ab.cs.uni-tuebingen.de/download/megan6/) | Interactive LCA taxonomy + functional (SEED/KEGG/InterPro) analysis from BLAST/DIAMOND. |

## Gene Prediction & MAG Annotation

| Tool | Link | Description |
|---|---|---|
| Pyrodigal | [link](https://github.com/althonos/pyrodigal) | Python bindings/reimplementation for Prodigal-style gene prediction. |
| METABOLIC | [link](https://github.com/AnantharamanLab/METABOLIC) | Profiles metabolic and biogeochemical functional traits from genomes/MAGs. |
| Prodigal / metaProdigal | [link](https://github.com/hyattpd/Prodigal) | Predicts protein-coding genes in microbial genomes and metagenomic contigs. |
| Bakta | [link](https://github.com/oschwengers/bakta) | Rapid microbial genome and MAG annotation suite. |
| Prokka | [link](https://github.com/tseemann/prokka) | Rapid prokaryotic genome annotation. |
| DRAM | [link](https://github.com/WrightonLabCSU/DRAM) | Distilled and Refined Annotation of Metabolism for microbial genomes/MAGs. |

## Gene Prediction & Viral Annotation

| Tool | Link | Description |
|---|---|---|
| Prodigal-GV | [link](https://github.com/apcamargo/prodigal-gv) | Fork of Prodigal/pyrodigal with models for giant viruses and viruses using alternative genetic codes. |
| PHANOTATE | [link](https://github.com/deprekate/PHANOTATE) | Phage-specific gene prediction tool designed for phage genomes. |
| pyrodigal-gv | [link](https://pypi.org/project/pyrodigal-gv/) | Python extension/wrapper implementing Prodigal-GV models for viral and giant-virus gene prediction. |

## Generative Sequence / Protein Design

| Tool | Link | Description |
|---|---|---|
| ProteinMPNN | [link](https://github.com/dauparas/ProteinMPNN) | Message-passing neural network for sequence design from backbone. |
| ProGen2 | [link](https://github.com/enijkamp/progen2) | Scaled ProGen with larger model sizes and broader training. |
| ESM-IF | [link](https://github.com/facebookresearch/esm) | Inverse folding predicts sequence given backbone structure (ESM Inverse Folding). |
| RFdiffusion | [link](https://github.com/RosettaCommons/RFdiffusion) | Diffusion model for de novo protein backbone generation from RoseTTAFold. |
| ProGen | [link](https://github.com/salesforce/progen) | Conditional protein language model that generates functional sequences. |

## Genome Size Estimation

| Tool | Link | Description |
|---|---|---|
| MicrobeCensus | [link](https://github.com/snayfach/MicrobeCensus) | Estimates average genome size and genome-equivalent abundance from shotgun metagenomes. |

## Genome-Scale Metabolic Modelling

| Tool | Link | Description |
|---|---|---|
| AGORA2 | [link](https://www.vmh.life/) | Curated GEM resource of ~7,300 human-gut bacterial strains with drug metabolism. |
| APOLLO | [link](https://www.vmh.life/) | Curated GEM resource of ~7,300 human-gut bacterial strains with drug metabolism. |
| COMETS | [link](https://github.com/segrelab/comets) | Computational modelling of metabolic outputs over time. |
| CarveMe | [link](https://github.com/cdanielmachado/carveme) | Top-down GEM reconstruction using a universal template + genome evidence. |
| SMETANA | [link](https://github.com/cdanielmachado/smetana) | Community metabolic interaction analysis via constraint-based modelling. |
| BacArena | [link](https://github.com/euba/BacArena) | Agent-based simulation of microbial communities on a spatial grid using GEMs. |
| metaGEM | [link](https://github.com/franciscozorrilla/metaGEM) | Snakemake workflow from reads → MAGs → community GEMs. |
| gapseq | [link](https://github.com/jotech/gapseq) | Automated reconstruction of GEMs from genome with metabolite-pathway inference and gap-filling. |
| COBRApy | [link](https://github.com/opencobra/cobrapy) | Python library for constraint-based modelling (FBA, FVA, pFBA, etc.). |
| ModelSEED | [link](https://modelseed.org/) | Reaction/genome database + reconstruction service (also via KBase). |
| KBase | [link](https://www.kbase.us/) | DOE cloud platform for genome assembly, annotation, metabolic modelling. |

## Hierarchical Testing

| Tool | Link | Description |
|---|---|---|
| treeclimbR | [link](https://github.com/fionarhuang/treeclimbR) | Finds data-dependent resolution levels for hierarchical hypotheses, useful with taxonomic trees. |

## Host Read Removal

| Tool | Link | Description |
|---|---|---|
| hostile | [link](https://github.com/bede/hostile) | Curated human/host index plus CLI for fast, reproducible host-read depletion with read scrubbing. |
| Bowtie2 | [link](https://github.com/BenLangmead/bowtie2) | Short-read aligner commonly used to map and remove host (human) reads. |
| BWA-MEM2 | [link](https://github.com/bwa-mem2/bwa-mem2) | Faster reimplementation of BWA-MEM for short-read alignment to host reference. |
| minimap2 | [link](https://github.com/lh3/minimap2) | Versatile short/long-read aligner; popular for host removal with long-read or hybrid data. |

## In Situ Growth Rate / Replication Rate

| Tool | Link | Description |
|---|---|---|
| PTRC | [link](https://genie.weizmann.ac.il/software/bacterial_growth.html) | PTR Calculator — original Korem/Segal method estimating replication rate from coverage. |
| iRep | [link](https://github.com/christophertbrown/iRep) | Original index-of-replication tool estimating PTR from genome-wide coverage of draft genomes. |
| GRiD | [link](https://github.com/ohlab/GRiD) | Estimates peak-to-trough ratio (PTR) from short-read coverage vs reference genomes. |
| SMEG | [link](https://github.com/ohlab/SMEG) | Species-level growth estimation from MAGs using SNV-based replication peak detection. |
| CoPTR | [link](https://github.com/tyjo/coptr) | Coverage-based PTR estimation with statistical confidence intervals and reference DB. |
| DEMIC | [link](https://sourceforge.net/projects/demic/) | PTR estimation without complete reference using contig coverage patterns. |

## Interactive Visualisation

| Tool | Link | Description |
|---|---|---|
| iSEEtree | [link](https://bioconductor.org/packages/iSEEtree) | Interactive explorer for hierarchical data, including tree-structured microbiome features. |
| miaDash | [link](https://miadash-microbiome.2.rahtiapp.fi/) | Interactive Shiny dashboard exposing mia functionality for microbiome exploration without extensive R coding. |

## K-mer / Reference-Free Comparison

| Tool | Link | Description |
|---|---|---|
| Simka | [link](https://github.com/GATB/simka) | Reference-free comparative metagenomics using k-mer based ecological distances. |
| MASH screen / distance | [link](https://github.com/marbl/Mash) | MinHash sketching for genome/metagenome distance and containment-like comparisons. |

## Longitudinal & Time-Series Analysis

| Tool | Link | Description |
|---|---|---|
| microSTASIS | [link](https://bioconductor.org/packages/microSTASIS) | Microbiota stability assessment across longitudinal samples using iterative clustering. |
| miaTime | [link](https://microbiome.github.io/miaTime/) | R package for manipulation and analysis of microbial community time-series data using TreeSummarizedExperiment. |

## Machine Learning & Prediction

| Tool | Link | Description |
|---|---|---|
| IntegratedLearner | [link](https://github.com/ocbe-uio/IntegratedLearner) | R package for multi-omics classification and prediction workflows. |
| mikropml | [link](https://github.com/SchlossLab/mikropml) | R package for supervised machine learning pipelines on microbiome data. |
| SIAMCAT | [link](https://siamcat.embl.de/) | Statistical inference and machine learning framework for microbiome case-control datasets. |

## MAG Dereplication & Genome Catalogues

| Tool | Link | Description |
|---|---|---|
| dRep | [link](https://github.com/MrOlm/drep) | Dereplicates genomes/MAGs and compares genome quality/ANI. |
| FastANI | [link](https://github.com/ParBLiSS/FastANI) | Fast whole-genome average nucleotide identity estimation. |

## MAG Quality Assessment

| Tool | Link | Description |
|---|---|---|
| CheckV | [link](https://bitbucket.org/berkeleylab/checkv) | Quality assessment for viral genomes (completeness, contamination, provirus boundaries). |
| BUSCO | [link](https://busco.ezlab.org/) | Completeness assessment via lineage-specific single-copy orthologs. |
| CheckM2 | [link](https://github.com/chklovski/CheckM2) | Machine-learning predictor of MAG completeness and contamination. |
| GUNC | [link](https://github.com/grp-bork/gunc) | Chimerism detector flagging contigs from different lineages in a MAG. |

## Meta-Analysis & Batch Adjustment

| Tool | Link | Description |
|---|---|---|
| MMUPHin | [link](https://github.com/biobakery/MMUPHin) | Meta-analysis and batch correction framework for microbial community profiles. |

## Metabolite / Metagenome Inference

| Tool | Link | Description |
|---|---|---|
| MIMOSA2 | [link](https://github.com/borenstein-lab/MIMOSA2) | Links community metabolic potential to metabolomic variation using metabolic network information. |
| MelonnPan | [link](https://huttenhower.sph.harvard.edu/melonnpan/) | Predicts metabolite profiles from microbial community features using elastic-net models. |

## Metagenomics Workflows & Pipelines

| Tool | Link | Description |
|---|---|---|
| MAG_Snakemake_wf | [link](https://github.com/Finn-Lab/MAG_Snakemake_wf) | Snakemake workflow for recovery and quality assessment of prokaryotic MAGs from short-read host-associated metagenomic datasets. |
| SnakeWRAP | [link](https://github.com/jsede/SnakeWRAP) | Snakemake workflow automating multiple metaWRAP modules for high-throughput shotgun metagenome assembly and analysis. |
| OpusTaxa | [link](https://github.com/yenkaiC/OpusTaxa) | Open-source Snakemake workflow for short-read shotgun metagenomics combining QC, host removal, taxonomic profiling, assembly, functional analysis, SRA input, database setup, and harmonised cross-sample outputs. |

## Microbial GWAS & Pangenome Association

| Tool | Link | Description |
|---|---|---|
| Scoary | [link](https://github.com/AdmiralenOla/Scoary) | Pan-genome wide association testing on gene presence/absence. |
| Anpan | [link](https://github.com/biobakery/anpan) | Microbial gene-content association testing across strains within species. |
| treeWAS | [link](https://github.com/caitiecollins/treeWAS) | Phylogeny-aware association testing controlling for clonal structure. |
| Panaroo | [link](https://github.com/gtonkinhill/panaroo) | Quality-controlled pangenome inference with error correction. |
| PPanGGOLiN | [link](https://github.com/labgem/PPanGGOLiN) | Partitioned pangenome graphs combining gene families and synteny. |
| Pyseer | [link](https://github.com/mgalardini/pyseer) | Microbial GWAS framework supporting SNPs, k-mers, and gene presence/absence. |
| Roary | [link](https://github.com/sanger-pathogens/Roary) | Fast core/accessory pangenome construction. |
| PIRATE | [link](https://github.com/SionBayliss/PIRATE) | Pangenome construction across multiple identity thresholds. |

## Microbial Load

| Tool | Link | Description |
|---|---|---|
| Microbial Load Predictor (MLP) | [link](https://github.com/grp-bork/microbial_load_predictor) | Predicts fecal microbial load/cell density from relative-abundance taxonomic profiles. |
| DspikeIn | [link](https://bioconductor.org/packages/release/bioc/html/DspikeIn.html) | Estimating absolute abundance from microbial spike-in controls. |

## Microbiome Biological-Age / Clock Models

| Tool | Link | Description |
|---|---|---|
| Galkin gut microbiome aging clock | [link](https://github.com/InSilicoMedicine/microbiome_clock) | Deep neural network predicting host chronological age from gut 16S/shotgun profiles. |

## Mobilome / Plasmids / MGEs

| Tool | Link | Description |
|---|---|---|
| PlasmidFinder | [link](https://bitbucket.org/genomicepidemiology/plasmidfinder) | Replicon-based identification of plasmid incompatibility groups. |
| ICEberg | [link](https://db-mml.sjtu.edu.cn/ICEberg/) | Database of integrative and conjugative elements (ICEs). |
| PlasX | [link](https://github.com/michaelkyu/PlasX) | Plasmid sequence classifier via gene-family features and logistic regression. |
| MOB-suite | [link](https://github.com/phac-nml/mob-suite) | Plasmid reconstruction, typing, and mobility prediction from short reads. |
| PlasClass | [link](https://github.com/Shamir-Lab/PlasClass) | Classifies plasmid-derived sequences using machine learning. |

## Multi-Omics Integration

| Tool | Link | Description |
|---|---|---|
| DIABLO (mixOmics) | [link](http://mixomics.org/mixdiablo/) | Sparse multi-block PLS-DA for supervised multi-omics integration. |
| MultiAssayExperiment | [link](https://bioconductor.org/packages/MultiAssayExperiment/) | Bioconductor data structure for coordinated multi-omics objects. |
| MMvec | [link](https://github.com/biocore/mmvec) | Neural network estimating microbe–metabolite co-occurrence probabilities. |
| MOFA2 | [link](https://github.com/bioFAM/MOFA2) | Multi-Omics Factor Analysis (Bayesian group factor model). |

## Network & Co-occurrence Inference

| Tool | Link | Description |
|---|---|---|
| CCREPE | [link](https://bioconductor.org/packages/ccrepe/) | Compositional re-sampling for empirical p-values on similarity measures. |
| SparCC | [link](https://github.com/JCSzamosi/SparCC3) | Sparse Correlations for Compositional data via log-ratio approximation. |
| FlashWeave | [link](https://github.com/meringlab/FlashWeave.jl) | Julia implementation of fast conditional dependency inference for microbiome data. |
| NetCoMi | [link](https://github.com/stefpeschel/NetCoMi) | R framework for network construction and comparison across groups. |
| mLDM | [link](https://github.com/tinglab/mLDM) | Metagenomic Latent Dirichlet Multinomial model for taxa interactions and environmental factors. |
| SPIEC-EASI | [link](https://github.com/zdk123/SpiecEasi) | Sparse inverse covariance estimation for compositional microbiome networks. |
| SPRING | [link](https://github.com/GraceYoon/SPRING) | Semi-parametric rank-based approach for inferring microbial association networks. |

## Phage / Virus Discovery

| Tool | Link | Description |
|---|---|---|
| VIBRANT | [link](https://github.com/AnantharamanLab/VIBRANT) | Identification + annotation of viral genomes with HMM + neural network. |
| geNomad | [link](https://github.com/apcamargo/genomad) | Hybrid (sequence + alignment) classifier for plasmids and viruses. |
| DeepVirFinder | [link](https://github.com/jessieren/DeepVirFinder) | CNN-based viral sequence detection from contigs. |
| VirFinder | [link](https://github.com/jessieren/virfinder) | k-mer based prediction of viral contigs from metagenomic assemblies. |
| VirSorter2 | [link](https://github.com/jiarong/VirSorter2) | Multi-classifier viral sequence identification from metagenomic contigs. |
| MARVEL | [link](https://github.com/LaboratorioBioinformatica/MARVEL) | Phage MAG identification via gene-based features and random forest. |
| Cenote-Taker2 | [link](https://github.com/mtisza1/Cenote-Taker2) | Automated virus discovery and annotation from contigs. |
| DRAM-v | [link](https://github.com/WrightonLabCSU/DRAM) | Metabolism-aware annotation of viral contigs, including auxiliary metabolic genes. |

## Phage Host Prediction

| Tool | Link | Description |
|---|---|---|
| CRISPRCasTyper | [link](https://github.com/Russel88/CRISPRCasTyper) | Detects and types CRISPR-Cas systems. |
| SpacePHARER | [link](https://github.com/soedinglab/spacepharer) | CRISPR spacer phage matching tool with sensitivity for fragmented hits. |
| RaFAH | [link](https://sourceforge.net/projects/rafah/) | Random forest host prediction from viral proteomes. |
| iPHoP | [link](https://bitbucket.org/srouxjgi/iphop) | Integrated host prediction combining alignment, CRISPR, and ML signals. |
| CRISPRopenDB | [link](https://github.com/edzuf/CrisprOpenDB) | Curated CRISPR spacer database with web search for phage host links. |
| vHULK | [link](https://github.com/LaboratorioBioinformatica/vHULK) | Deep-learning host prediction from phage protein content. |
| PHIST | [link](https://github.com/refresh-bio/PHIST) | Fast exact k-mer matching between phage and host genomes. |

## Phage Therapy & Viral Workflows

| Tool | Link | Description |
|---|---|---|
| Sphae | [link](https://github.com/linsalrob/sphae) | Automated Snakemake/Snaketool workflow to assemble, annotate, and evaluate phage candidates for therapeutic suitability. |

## Phylogenetics & Evolution

| Tool | Link | Description |
|---|---|---|
| TreeTime | [link](https://github.com/neherlab/treetime) | Maximum-likelihood phylodynamic and time-scaled phylogeny inference. |
| Gubbins | [link](https://github.com/nickjcroucher/gubbins) | Detects recombination in bacterial whole-genome alignments. |

## Prophage Detection

| Tool | Link | Description |
|---|---|---|
| PhiSpy | [link](https://github.com/linsalrob/PhiSpy) | Predicts prophage regions in bacterial/archaeal genomes using similarity- and composition-based features. |

## R / Bioconductor Analysis Frameworks

| Tool | Link | Description |
|---|---|---|
| MicrobiotaProcess | [link](https://bioconductor.org/packages/MicrobiotaProcess) | Tidy-style framework for deep mining and visualisation of microbiome and ecological data. |
| microeco | [link](https://cran.r-project.org/package=microeco) | R framework for data mining and ecological analysis of microbiome communities. |
| mia | [link](https://microbiome.github.io/mia/) | Core miaverse package for microbiome data wrangling, transformations, diversity, ordination and downstream analysis using TreeSummarizedExperiment. |
| miaViz | [link](https://microbiome.github.io/miaViz/) | Visualisation package for microbiome data, especially TreeSummarizedExperiment-based objects. |
| bugsigdbr | [link](https://bioconductor.org/packages/bugsigdbr/) | R/Bioconductor client to query and analyse BugSigDB signatures. |
| TaxSEA | [link](https://bioconductor.org/packages/TaxSEA) | Taxon set enrichment analysis for ranked microbiome differential abundance or correlation results. |

## Read QC & Preprocessing

| Tool | Link | Description |
|---|---|---|
| SeqKit | [link](https://bioinf.shenwei.me/seqkit/) | Fast command-line toolkit for FASTA/Q manipulation and summary statistics. |
| Cutadapt | [link](https://cutadapt.readthedocs.io/) | Adapter trimming and demultiplexing for high-throughput sequencing reads. |
| KneadData | [link](https://github.com/biobakery/kneaddata) | bioBakery wrapper combining Trimmomatic QC with Bowtie2 host-read depletion. |
| MultiQC | [link](https://github.com/MultiQC/MultiQC) | Aggregates QC outputs (FastQC, fastp, Bowtie2, Kraken2, etc.) into a single HTML report. |
| fastp | [link](https://github.com/OpenGene/fastp) | All-in-one FASTQ QC, adapter trimming, polyG/polyX clipping, dedup, length/quality filtering. |
| FastQC | [link](https://github.com/s-andrews/FastQC) | Per-sample QC report on base quality, GC, adapter content, duplication. |
| Trimmomatic | [link](https://github.com/usadellab/Trimmomatic) | Adapter and quality trimming for Illumina reads with sliding-window and ILLUMINACLIP modes. |
| BBDuk (BBTools) | [link](https://sourceforge.net/projects/bbmap/) | K-mer based adapter trimming, contaminant filtering, quality trimming, kmer matching. |
| FastQ Screen | [link](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) | Screens reads against multiple reference genomes to detect contamination or sample composition issues. |

## Reference Databases — AMR

| Resource | Link | Description |
|---|---|---|
| ResFinder DB | [link](https://bitbucket.org/genomicepidemiology/resfinder_db) | DTU AMR reference DB; acquired-resistance focus. |
| CARD | [link](https://card.mcmaster.ca/) | Comprehensive Antibiotic Resistance Database with curated ARO ontology. |
| ARG-ANNOT | [link](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/) | Antibiotic Resistance Gene-ANNOTation database (Bordeaux). |

## Reference Databases — BGCs

| Resource | Link | Description |
|---|---|---|
| MIBiG | [link](https://mibig.secondarymetabolites.org/) | Minimum Information about a Biosynthetic Gene cluster curated BGC reference. |

## Reference Databases — Curated Cohorts

| Resource | Link | Description |
|---|---|---|
| GMrepo v2 | [link](https://gmrepo.humangut.info/) | Curated database of human gut metagenomes and phenotypes. |
| curatedMetagenomicData | [link](https://waldronlab.io/curatedMetagenomicData/) | Curated, uniformly processed human metagenome datasets for R/Bioconductor. |
| microbiome-metabolome-curated-data | [link](https://github.com/borenstein-lab/microbiome-metabolome-curated-data) | Paired datasets of metagenomes and metabolomes. |

## Reference Databases — Function

| Resource | Link | Description |
|---|---|---|
| eggNOG | [link](http://eggnog5.embl.de/) | Orthologous group database with functional annotations across taxa. |
| CAZy | [link](http://www.cazy.org/) | Carbohydrate-Active enZYmes database (GHs, GTs, PLs, CEs, AAs). |
| MetaCyc | [link](https://metacyc.org/) | Manually curated metabolic pathway database. |
| Pfam | [link](https://www.ebi.ac.uk/interpro/) | Protein family HMMs of conserved domains. |
| InterPro | [link](https://www.ebi.ac.uk/interpro/) | Integrated protein signatures across Pfam, SMART, PROSITE, PANTHER, etc. |
| KEGG | [link](https://www.kegg.jp/) | Curated pathway, ortholog (KO), and reaction database. |

## Reference Databases — MAG Catalogues

| Resource | Link | Description |
|---|---|---|
| HRGM2 | [link](https://www.hrgm.org/) | Human Reference Gut Microbiome non-redundant genome catalogue. |
| HROM | [link](https://www.hrom.org/) | Human oral MAG catalogue. |
| gcMeta | [link](https://gcmeta.wdcm.org/) | Multi-biome catalogue/database portal with public catalogue APIs and archive bundles. |
| GEM catalog | [link](https://portal.nersc.gov/GEM/) | Global Earth Microbiome catalogue of bacterial and archaeal MAGs across environments. |
| OMDB | [link](https://omdb.genzentrum.lmu.de/) | Ocean Microbiome Database with reconstructed genomes, genes and scaffold catalogues. |
| cFMD | [link](https://github.com/SegataLab/cFMD) | Curated Food Metagenomic Data resource with food shotgun metagenomes, MAGs and HUMAnN outputs. |
| LakePulse MAG Catalogue | [link](https://datadryad.org/) | Canadian freshwater and oligosaline lake MAG catalogue. |
| MRGM | [link](https://www.mrgm.org/) | Mouse gut MAG catalogue and companion profiling databases. |
| RUG2 | [link](https://datashare.ed.ac.uk/) | Cattle rumen MAG catalogue and companion annotations. |
| SMAG | [link](https://microbma.github.io/project/SMAG.html) | Global Soil MAG catalogue from thousands of soil metagenomes. |
| UHGG | [link](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0) | Unified Human Gastrointestinal Genome catalog (>200k genomes from human gut). |

## Reference Databases — Phages / Viruses

| Resource | Link | Description |
|---|---|---|
| IMG/VR | [link](https://img.jgi.doe.gov/vr/) | JGI's integrated viral genome database from cultivated and uncultivated sources. |
| MGV | [link](https://portal.nersc.gov/MGV/) | Metagenomic Gut Virus catalog of >180k high-quality gut viral genomes. |
| RVDB | [link](https://rvdb.dbi.udel.edu/) | Reference Viral DataBase curated viral protein/nucleotide reference. |
| Meta-virus Resource (MetaVR) | [link](https://metavirus-resource.org/) | Multi-biome database of uncultivated viral genomes, vOTUs and protein clusters. |

## Reference Databases — Taxa

| Resource | Link | Description |
|---|---|---|
| MetaPhlAn4 + SGB database | [link](https://github.com/biobakery/MetaPhlAn) | Species-level genome bin marker database underlying MetaPhlAn4. |
| GTDB | [link](https://gtdb.ecogenomic.org/) | Genome Taxonomy Database standardised prokaryotic taxonomy from genome phylogeny. |
| MGnify | [link](https://www.ebi.ac.uk/metagenomics/) | EMBL-EBI resource for metagenomic datasets, assemblies, genomes and analysis outputs. |
| NCBI RefSeq | [link](https://www.ncbi.nlm.nih.gov/refseq/) | Curated reference genomes and taxonomy maintained by NCBI. |

## Reference Databases — Other

| Resource | Link | Description |
|---|---|---|
| BacDive | [link](https://bacdive.dsmz.de/) | Bacterial diversity metadatabase with phenotypic and physiological traits. |
| MiMeDB | [link](https://mimedb.org/) | Microbe-metabolite interaction database. |
| mBodyMap | [link](https://mbodymap.microbiome.cloud/) | Body-site resolved human microbiome abundance and prevalence reference. |
| VFDB | [link](http://www.mgc.ac.cn/VFs/) | Virulence Factor Database curated virulence factor gene catalog. |

## Sequencing Depth & Coverage

| Tool | Link | Description |
|---|---|---|
| Mosdepth | [link](https://github.com/brentp/mosdepth) | Fast depth/coverage calculation from BAM/CRAM files. |
| CoverM | [link](https://github.com/wwood/CoverM) | Calculates read coverage of contigs, genomes and MAGs from metagenomic mapping. |
| Nonpareil | [link](https://nonpareil.readthedocs.io/) | Estimates metagenomic coverage and projects sequencing effort using read redundancy. |

## Simulation & Benchmarking

| Tool | Link | Description |
|---|---|---|
| miaSim | [link](https://microbiome.github.io/miaSim/) | Simulates longitudinal and time-series microbial community dynamics using ecological models such as gLV, Ricker, Hubbell neutral and consumer-resource models. |

## Source Tracking

| Tool | Link | Description |
|---|---|---|
| SourceTracker2 | [link](https://github.com/biota/sourcetracker2) | Bayesian microbial source tracking method. |
| FEAST | [link](https://github.com/cozygene/FEAST) | Fast expectation-maximisation source tracking for microbial communities. |

## Strain-Level Analysis

| Tool | Link | Description |
|---|---|---|
| StrainPhlAn4 | [link](https://github.com/biobakery/MetaPhlAn) | SNV-based strain phylogenies built from MetaPhlAn4 species marker reads. |
| StrainGE | [link](https://github.com/broadinstitute/StrainGE) | Reference-based strain detection and abundance at low coverage. |
| StrainFinder | [link](https://github.com/cssmillie/StrainFinder) | EM-based mixture model of strain frequencies at marker SNV sites. |
| MIDAS2 | [link](https://github.com/czbiohub-sf/MIDAS2) | Strain-level SNP and CNV profiling against UHGG species pangenomes. |
| SameStr | [link](https://github.com/danielpodlesny/SameStr) | Shared strain detection via SNV co-occurrence across samples. |
| inStrain | [link](https://github.com/MrOlm/inStrain) | Population genetics on MAGs: ANI, SNVs, popANI, conANI, gene content. |
| PStrain | [link](https://github.com/wshuai294/PStrain) | Iterative reconstruction of strain composition from SNV profiles. |
| TRACS | [link](https://github.com/gtonkinhill/tracs/) | Strain level tracking across samples. |

## Structural Metagenomics / Protein Structure at Scale

| Tool | Link | Description |
|---|---|---|
| ESM Metagenomic Atlas | [link](https://esmatlas.com/) | Database of >600M predicted structures from metagenomic ORFs (ESMFold). |
| ESMFold | [link](https://github.com/facebookresearch/esm) | Single-sequence protein structure prediction from ESM-2 protein LM. |

## Taxon Set / Signature Enrichment

| Tool | Link | Description |
|---|---|---|
| BugSigDB | [link](https://bugsigdb.org/) | Manually curated database of published microbial signatures across conditions. |
| EnrichM | [link](https://github.com/geronimp/enrichM) | Function enrichment in MAGs — KEGG/Pfam/CAZy/COG enrichment vs background. |
| MicrobiomeAnalyst TSEA | [link](https://www.microbiomeanalyst.ca/) | Taxon Set Enrichment Analysis module testing user lists against curated signatures. |

## Taxonomic Placement & MAG Taxonomy

| Tool | Link | Description |
|---|---|---|
| PhyloPhlAn | [link](https://github.com/biobakery/phylophlan) | Phylogenetic placement on universal marker tree (>400 markers); strain-level support. |
| skani | [link](https://github.com/bluenote-1577/skani) | Fast accurate ANI estimation for MAGs and large genome collections. |
| GTDB-Tk | [link](https://github.com/Ecogenomics/GTDBTk) | Classifies prokaryotic MAGs against the GTDB taxonomy using marker placements + ANI. |
| Mash | [link](https://github.com/marbl/Mash) | Original MinHash sketch tool for genome distance and screening. |
| CAT/BAT | [link](https://github.com/MGXlab/CAT_pack) | Contig Annotation Tool and Bin Annotation Tool for taxonomic classification of contigs and MAGs. |
| sourmash | [link](https://github.com/sourmash-bio/sourmash) | MinHash/FracMinHash sketches for fast ANI, classification, and search. |

## Taxonomic Profiling (Read-Based)

| Tool | Link | Description |
|---|---|---|
| CLARK | [link](http://clark.cs.ucr.edu/) | K-mer based taxonomic classification of metagenomic reads. |
| MetaPhlAn4 | [link](https://github.com/biobakery/MetaPhlAn) | Marker-gene based species- and SGB-level profiler using clade-specific markers. |
| Kaiju | [link](https://github.com/bioinformatics-centre/kaiju) | Protein-level read classification via translated BWT search; better for divergent taxa. |
| sylph | [link](https://github.com/bluenote-1577/sylph) | ANI-based metagenomic profiler using k-mer containment and zero-inflated statistics. |
| Centrifuge | [link](https://github.com/DaehwanKimLab/centrifuge) | BWT/FM-index based classifier compressing many references for taxonomic assignment. |
| Kraken2 | [link](https://github.com/DerrickWood/kraken2) | K-mer based exact-match LCA taxonomic classifier of reads against a hash database. |
| KrakenUniq | [link](https://github.com/fbreitwieser/krakenuniq) | Kraken with HyperLogLog unique-k-mer counts per taxon to flag false positives. |
| Bracken | [link](https://github.com/jenniferlu717/Bracken) | Bayesian re-estimation of species/genus abundance from Kraken2 read assignments. |
| Recentrifuge | [link](https://github.com/khyox/recentrifuge) | Interactive and statistical comparative analysis of Centrifuge/Kraken-style outputs. |
| mOTUs3 | [link](https://github.com/motu-tool/mOTUs) | Universal single-copy marker-gene-based operational taxonomic unit profiler. |
| KMCP | [link](https://github.com/shenwei356/kmcp) | Accurate metagenomic profiler with COBS-style k-mer index and confidence scoring. |
| SingleM | [link](https://wwood.github.io/singlem/) | Marker-gene based profiling of shotgun short- and long-read metagenomes, with strength for novel lineages. |

## Temporal & Dynamic Modelling

| Tool | Link | Description |
|---|---|---|
| BEEM | [link](https://github.com/csb5/BEEM) | Estimate gLV parameters from cross-sectional or short time-series data. |
| MDSINE2 | [link](https://github.com/gerberlab/MDSINE2) | Bayesian generalised Lotka–Volterra inference with regularisation and perturbations. |
| TIME | [link](https://web.rniapps.net/time/) | Web tool for time-resolved microbiome data exploration and trend testing. |

## Viral & Virome Analysis

| Tool | Link | Description |
|---|---|---|
| VirSorter (legacy) | [link](https://bitbucket.org/MAVERICLab/virsorter) | Original viral sequence detection tool based on viral hallmark genes and depletion of cellular genes. |
| ViralVerify | [link](https://github.com/ablab/viralVerify) | Classifies metagenomic contigs as viral, bacterial, or uncertain using gene prediction, HMMs and a Naive Bayes classifier. |
| PhaMers | [link](https://github.com/jondeaton/PhaMers) | Supervised k-mer based classifier for identifying bacteriophage sequences from metagenomic contigs. |
| PhaBOX | [link](https://github.com/KennthShang/PhaBOX) | Web server/workflow for phage contig identification, lifestyle prediction, taxonomic classification, and host prediction. |
| Phage Galaxy | [link](https://phage.usegalaxy.eu/) | Galaxy-based platform containing workflows and tools for phage assembly, annotation and comparative genomics. |
| Empathi | [link](https://github.com/Alexandre-Boulay/Empathi) | Embedding-based hierarchical phage protein annotation tool. |
| Pharokka | [link](https://github.com/gbouras13/pharokka) | Rapid, standardised annotation of bacteriophage genomes and metagenomes. |
| Phold | [link](https://github.com/gbouras13/phold) | Protein structure-informed bacteriophage genome annotation using ProstT5 and Foldseek against curated phage protein structures. |
| PhaGO | [link](https://github.com/JiaojiaoGuan/PhaGO) | Phage protein function annotation integrating protein embeddings with genomic context. |
| PhaVIP | [link](https://github.com/KennthShang/PhaVIP) | Phage virion protein classification using chaos game representation and vision transformer models. |
| Phynteny | [link](https://github.com/susiegriggo/Phynteny) | Synteny-aware phage protein function prediction using conserved gene order across phages. |
| multiPhATE2 | [link](https://github.com/carolzhou/multiPhATE2) | Automated functional annotation and comparison pipeline for phage genomes. |
| Phables | [link](https://github.com/Vini2/phables) | Resolves complete phage genomes from fragmented viral metagenome assembly graphs using graph decomposition. |
| VirHostMatcher-Net | [link](https://github.com/WeiliWw/VirHostMatcher-Net) | Network/deep-learning style host prediction for viruses/phages based on sequence features. |
| ViralQC | [link](https://github.com/KennthShang/ViralQC) | Assesses completeness and contamination of predicted viral contigs/bins using foundation-model-style sequence information. |
| PhaGCN | [link](https://github.com/KennthShang/PhaGCN) | Graph convolutional network classifier for taxonomic classification of phage contigs. |
| VITAP | [link](https://github.com/KennthShang/VITAP) | Viral Taxonomic Assignment Pipeline for DNA and RNA virus classification. |
| vConTACT2 | [link](https://github.com/Hocnonsense/vcontact2) | Gene-sharing network approach for viral taxonomy and clustering of viral genomes/contigs. |
| PRFect | [link](https://github.com/deprekate/prfect) | Machine-learning method to predict programmed ribosomal frameshifts in coding genes. |

## Visualisation

| Tool | Link | Description |
|---|---|---|
| Microviz | [link](https://david-barnett.github.io/microViz/) | R package for ordination, taxonomic bar charts, heatmaps on phyloseq objects. |
| GraPhlAn | [link](https://github.com/biobakery/graphlan) | Compact circular phylogenetic and taxonomic visualisation. |
| q2-qurro / Qurro | [link](https://github.com/biocore/qurro) | Interactive rank-plot visualisation for log-ratio analyses. |
| Pavian | [link](https://github.com/fbreitwieser/pavian) | Shiny app for interactive exploration of Kraken/KrakenUniq/Bracken reports. |
| Krona | [link](https://github.com/marbl/Krona) | Interactive hierarchical (multi-level pie) visualisation of taxonomic abundance. |
| iTOL | [link](https://itol.embl.de/) | Interactive Tree Of Life web-based phylogenetic tree annotation. |

## Workflow & Pipeline Frameworks

| Tool | Link | Description |
|---|---|---|
| ANVI'O | [link](https://anvio.org/) | Interactive platform for microbial pan-genomics, metagenomics, MAG curation and visualisation. |
| MetaWRAP | [link](https://github.com/bxlab/metaWRAP) | Modular MAG-recovery pipeline wrapping multiple binners and refiners. |
| ATLAS | [link](https://github.com/metagenome-atlas/atlas) | Snakemake metagenomics pipeline producing MAGs with binning, refinement, annotation. |
| Aviary | [link](https://github.com/rhysnewell/aviary) | Workflow for metagenome assembly, binning, quality control and genome recovery. |
| Sunbeam | [link](https://github.com/sunbeam-labs/sunbeam) | Snakemake metagenomics pipeline (QC, host removal, classification, assembly). |
| bioBakery3/4 | [link](https://huttenhower.sph.harvard.edu/biobakery_workflows/) | Huttenhower lab suite: KneadData, MetaPhlAn, HUMAnN, StrainPhlAn, MaAsLin. |
| nf-core/mag | [link](https://nf-co.re/mag) | Nextflow pipeline for shotgun metagenomics: QC, assembly, binning, taxonomy, AMR. |
| nf-core/taxprofiler | [link](https://nf-co.re/taxprofiler) | Nextflow pipeline for taxonomic profiling of shotgun/metabarcoding data with multiple classifiers. |
| Snakemake | [link](https://snakemake.readthedocs.io/) | General reproducible workflow-management system widely used in bioinformatics. |
| Nextflow | [link](https://www.nextflow.io/) | Workflow-management system underlying nf-core and many cloud/HPC metagenomic workflows. |

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) — suggestions via Issues are very welcome.

## License

[MIT](LICENSE)

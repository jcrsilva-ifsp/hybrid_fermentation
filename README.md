# hybrid-fermentation-metagenomics

Genome-resolved shotgun metagenomics pipeline for spontaneous/mixed fermentation systems — recovers **both bacterial and fungal** MAGs, functionally annotates them, and identifies fungal species by genome-level ANI rather than protein homology.

Developed at CCBL (FCFRP-USP), under supervision of Prof. Ricardo Roberto da Silva.

[Leia em Português](README.pt-BR.md)

## Motivation

Standard metagenomic binning pipelines (e.g. MetaWRAP with CONCOCT/MaxBin2/MetaBAT2) were built almost exclusively for prokaryotes. Fungal genomes — larger, with introns and different codon usage — tend to be fragmented or discarded during binning. In beverages where the fungal fraction is biologically central (yeasts in Pulque, Kombucha, Lambic), this produces a systematic blind spot: the pipeline never recovers a single fungal MAG, even when the DNA is present in abundance.

This repository implements and validates a fix: domain bifurcation (Tiara) before binning, so bacteria and fungi are processed by domain-appropriate tools on separate branches. Validated against three independent public datasets (Pulque, Kombucha, Lambic — see below).

A second, independent issue was found and fixed during validation: taxonomic identification of fungal MAGs by protein homology against Swiss-Prot systematically misidentifies fermentation yeasts that are under-represented in that database (e.g. *Brettanomyces bruxellensis* has essentially no nuclear-genome entries in Swiss-Prot). Genome-level ANI comparison (fastANI) against a small reference-genome set resolves this and does not depend on reference-database annotation coverage.

## Pipeline overview

| Phase | Tool(s) | Output |
|---|---|---|
| 1 | Trimmomatic (+ optional bmtagger, paired-end only) | Clean reads |
| 1c | Kraken2 + Bracken | Direct taxonomic profile |
| 2 | MEGAHIT | Assembled contigs |
| 3 | Tiara | Bifurcated contigs (prokaryotic / eukaryotic) |
| 4A | MetaWRAP (CONCOCT+MaxBin2+MetaBAT2) + CheckM | Bacterial MAGs |
| 4B | MetaBAT2 + BUSCO (fungi_odb12.2) | Fungal MAGs |
| 5 | Prodigal + PROKKA (bacteria); MetaEuk (fungi); Kofamscan/KEGG (both) | Functional annotation |
| 5b | fastANI vs. reference genome set | Fungal species ID |
| 6 | antiSMASH + dbCAN | BGC/CAZyme prospection (not configured in this version) |

Two script variants: `pipeline_hibrido_corrigido.sh` (paired-end, Illumina) and `pipeline_hibrido_singleend.sh` (single-end, e.g. Ion Torrent).

## Installation

Six isolated conda/micromamba environments are required (see `environments/*.yml`).

```bash
git clone <REPO_URL>
cd hybrid-fermentation-metagenomics

for ENV in tiara metaeuk busco kofamscan bracken fastani; do
    micromamba env create -p ./env_${ENV} -f environments/env_${ENV}.yml
done
```

MetaWRAP is **not** reliably installable via conda (the `metawrap` conda package name collides with an unrelated Python package). Install from source:

```bash
git clone https://github.com/bxlab/metaWRAP.git
# follow the repository's own installation instructions (config-metawrap)
```

### Reference data

| Data | Used by | Source |
|---|---|---|
| Swiss-Prot (MMseqs2/MetaEuk-indexed) | Fungal gene prediction by homology | UniProt |
| fungi_odb12.2 | Fungal MAG completeness (BUSCO) | BUSCO official downloads |
| KOfam profiles + ko_list | KEGG functional mapping | ftp.genome.jp/pub/db/kofam/ |
| Kraken2 PlusPF | Direct taxonomic profiling | Standard Kraken2 collection |
| Fungal reference genomes for ANI ID (Fase 5b) | Fungal species identification | NCBI (fetched automatically on first run — see below) |

Fase 5b downloads its reference genomes automatically via `datasets` (NCBI Datasets CLI, included in `env_fastani`) the first time it runs, using the accessions configured at the top of the script:

```bash
ANI_REFS["Saccharomyces_cerevisiae"]="GCF_000146045.2"
ANI_REFS["Brettanomyces_bruxellensis"]="GCF_011074885.1"
ANI_REFS["Zygosaccharomyces_rouxii"]="GCA_000026365.1"
```

Add or remove entries to match the yeasts relevant to your system (e.g. *Pichia*, *Hanseniaspora*, *Dekkera*, lager-adapted strains).

### Reproducing the exact environment used in this work

```bash
micromamba activate -p ./env_<name>
conda env export | grep -v "^prefix: " > environments/env_<name>.yml
```

## Usage

Edit `BASE_DIR` and every `/path/to/...` placeholder at the top of the script for your environment. Then submit one job per sample:

```bash
sbatch --export=ALL,SAMPLE_ID=<accession> --job-name=metag_<accession> pipeline_hibrido_corrigido.sh
```

**Strongly recommended**: run one pilot sample first, following its log (`$BASE_DIR/logs/<accession>.log`), before submitting in batch. SLURM resource directives (partition, CPUs, memory) are cluster-specific and must be re-checked on any new system.

## Validation

Tested against three independent public datasets of spontaneous/mixed fermentation:

- **Pulque** (Chacón-Vargas et al., 2020, *Scientific Reports*, BioProject PRJNA603591): bacterial MAG taxonomy consistent with the reference article; fungal MAGs recovered and correctly identified as *Saccharomyces cerevisiae* (ANI ~87-88% against the S288C reference — below the conventional ≥95% species threshold, likely reflecting real lineage divergence, not identification uncertainty); *Zymomonas mobilis*, absent from binned MAGs, confirmed present and dynamically consistent with the article via Kraken2+Bracken.
- **Kombucha** (Suhre et al., 2025, *Food Bioscience*): fungal MAGs, originally misidentified as *S. cerevisiae* by protein homology, correctly re-identified as *Brettanomyces bruxellensis* by ANI (many hits ≥94-98%) — consistent with the reference article, which reports *Brettanomyces* dominating the fermented community (>88% in some samples).
- **Lambic** (De Roos et al., 2020, *Frontiers in Microbiology*, ENA BioProject PRJEB28363; Ion Torrent single-end): bacterial succession (Acetobacter-type → Lactobacillales-type across fermentation time) reproduces the acidification-then-maturation dynamic described in the article; fungal ID by ANI: 3/6 samples confirmed *S. cerevisiae* (ANI ~99%), 1/6 corrected to *Brettanomyces bruxellensis* (ANI 98.45%), 2/6 unresolved against the current reference set (candidate: *Pichia membranifaciens*, reported in the article but not yet included as a reference).

## Known limitations

- `metawrap bin_refinement` fails frequently (upstream MetaWRAP bug); a CheckM-based fallback handles this automatically.
- Not every sample yields a fungal MAG — absence is a valid result, not a failure.
- Fase 5b's reference-genome set is necessarily incomplete; a MAG with no confident ANI match may simply need a reference genome not yet included.
- Fase 6 (antiSMASH + dbCAN) is not configured in this version.
- Kraken2 PlusPF does not include plant genomes; in samples with plant matrix (e.g. fruit-based Kombucha), a large fraction of reads may be "unclassified" for this reason, not due to method failure.

## Repository structure

```
.
├── README.md
├── README.pt-BR.md
├── pipelines/
│   ├── pipeline_hibrido_corrigido.sh       (paired-end)
│   ├── pipeline_hibrido_singleend.sh       (single-end)
├── environments/
│   ├── env_tiara.yml
│   ├── env_metaeuk.yml
│   ├── env_busco.yml
│   ├── env_kofamscan.yml
│   ├── env_bracken.yml
│   └── env_fastani.yml
└── LICENSE
```

## Authorship

Jean Carlos Rodrigues da Silva, postdoctoral researcher — CCBL, FCFRP-USP.
Supervision: Prof. Ricardo Roberto da Silva.

## License

TBD (see `LICENSE`).

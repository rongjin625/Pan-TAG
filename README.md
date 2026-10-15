<a id="installation"></a>

## Getting Started

Pan-TAG supports Linux x86-64 and requires Git and [Pixi](https://pixi.sh/). Install on a node with internet access:

```bash
cd /path/to/software
git clone https://github.com/rongjin625/Pan-TAG.git
cd Pan-TAG
pixi install
```

For Liftoff annotation mode, also install:

```bash
pixi install --environment liftoff
```

Check the installation:

```bash
pixi run pantag --version
```

## Table of Contents

- [Getting Started](#installation)
- [Introduction](#introduction)
- [Pan-TAG Construction](#construction)
  - [Core Workflow](#core)
    - [Preparing Inputs](#inputs)
    - [Generating Configuration](#configuration)
    - [Running the Workflow](#analysis)
  - [Downstream Analysis](#downstream)
    - [Sequence Extraction](#extraction)
    - [Sequence Clustering and Visualization](#clustering)
- [Key Parameters](#parameters)
- [Outputs](#outputs)
- [Resuming Runs and Troubleshooting](#troubleshooting)
- [Citation and License](#citation-license)

<a id="introduction"></a>

## Introduction

Pan-TAG is a Snakemake-based tool for analyzing tandemly arrayed gene (TAG) clusters across a pangenome. Given genome assemblies and annotations from the same species, it identifies TAGs in individual samples, reconciles them across samples using flanking-gene synteny, and rescans each sample to validate clusters and quantify gene copy numbers.

The main outputs are a cross-sample TAG table, `pangenome_TAGs.tsv`, and a copy-number matrix, `refined_TAG_count_matrix.tsv`. Pan-TAG also supports Liftoff annotation transfer, pseudogene recovery, within-cluster sequence extraction and clustering, PCoA, and gene arrangement plots.

<a id="construction"></a>

## Pan-TAG Construction

<a id="core"></a>

### Core Workflow

Generate per-sample TAGs, a cross-sample TAG table, and a copy-number matrix from genome assemblies and annotations.

<a id="inputs"></a>

#### Preparing Inputs

##### Create a run directory

Multiple datasets can share the same code and environment, but each dataset should have its own run directory:

```bash
CODE_ROOT=/path/to/software/Pan-TAG
RUN_ROOT=/path/to/pantag_runs/my_dataset

mkdir -p "$RUN_ROOT/config"
cd "$RUN_ROOT"
```

Run all subsequent commands from `$RUN_ROOT`. Absolute paths are recommended for input FASTA and GFF files.

##### Prepare the sample table

Create `$RUN_ROOT/config/samples.tsv` as a tab-separated file with exactly these columns, in this order:

```text
sample_id	fasta	gff
```

Requirements:

- The first data row must be the reference sample.
- The reference sample must have both FASTA and GFF files.
- Each `sample_id` must be unique.
- FASTA sequence names must match the corresponding GFF `seqid` values in column 1.
- All samples in a run must belong to the same species.

**Use existing annotations (provided)**

Provide a GFF file for every sample:

```tsv
sample_id	fasta	gff
reference	/data/reference.fa	/data/reference.gff3
sample01	/data/sample01.fa	/data/sample01.gff3
sample02	/data/sample02.fa	/data/sample02.gff3
```

Use `annotation.mode: provided`.

**Transfer a common reference annotation (liftoff)**

Provide a GFF only for the reference sample in the first row. Leave the third column empty for all other samples:

| sample_id | fasta | gff |
|---|---|---|
| reference | `/data/reference.fa` | `/data/reference.gff3` |
| sample01 | `/data/sample01.fa` | Leave empty |
| sample02 | `/data/sample02.fa` | Leave empty |

Use `annotation.mode: liftoff`. Pan-TAG builds a database from the reference GFF and runs Liftoff for each non-reference sample. The two annotation modes cannot be mixed within a run.

The feature list for Liftoff's `-f` option is read from `$CODE_ROOT/config/liftoff_features.txt`, which currently contains `gene` and `pseudogene`. Modify this project file only if additional top-level feature types need to be transferred.

<a id="configuration"></a>

#### Generating Configuration

Generate a run-specific configuration from `samples.tsv`; do not use the example data paths in the repository template.

For existing GFF files:

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag configure \
  --samples "$RUN_ROOT/config/samples.tsv" \
  --annotation-mode provided \
  --pseudogene-enabled \
  --template "$CODE_ROOT/config/config.yaml" \
  --output "$RUN_ROOT/config/config.yaml" \
  --results-dir "$RUN_ROOT/results" \
  --logs-dir "$RUN_ROOT/logs" \
  --work-dir "$RUN_ROOT/work"
```

To transfer the reference annotation, replace `--annotation-mode provided` with `--annotation-mode liftoff` and leave the non-reference GFF fields empty as described above.

To disable pseudogene recovery, replace `--pseudogene-enabled` with `--no-pseudogene-enabled`. List all configuration options with:

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag configure -h
```

`configure` checks the sample table format, sample IDs, and the existence of all specified FASTA/GFF files before writing the configuration. It does not run the analysis or replace validation of file contents and annotation relationships.

<details>
<summary>New species or annotation sources: GFF attribute mapping and filtering</summary>

The default template accommodates common GENCODE, Ensembl, and standard GFF3 annotations:

```yaml
annotation:
  gene_feature_types: [gene]
  pseudogene_feature_types: [gene, pseudogene]
  transcript_feature_types: [mRNA, transcript]
  gene_id_keys: [gene_id, ID]
  gene_name_keys: [gene_name, Name, gene_symbol, gene, geneName]
  transcript_id_keys: [transcript_id, ID]
  transcript_parent_gene_keys: [gene_id, Parent]
  biotype_keys: [gene_type, gene_biotype, biotype]
  coding_filter_mode: auto
  coding_values: [protein_coding]
  excluded_attribute_terms: [readthrough_gene, readthrough]
  pseudogene_values: [pseudogene, processed_pseudogene, unprocessed_pseudogene]
  strip_id_prefixes: ["gene:", "transcript:"]
  fallback_gene_name_to_id: true
```

For a new annotation source, check feature types in GFF column 3 and the ID, Name, Parent, and biotype attributes in column 9, then adjust `config.yaml` accordingly. A record is treated as a pseudogene only if its feature type is in `pseudogene_feature_types` and its biotype matches `pseudogene_values`. `pseudogenic_transcript` is excluded by default.

Available `coding_filter_mode` values:

- `auto`: apply biotype filtering if a recognized coding biotype is present; otherwise select genes with translated proteins.
- `biotype`: require a biotype matching `coding_values`.
- `translated_protein`: select genes based on the availability of translated proteins, without coding-biotype filtering.

All three modes require an associated nonempty translated protein and apply `excluded_attribute_terms`. The longest protein is retained for each gene.

</details>

<a id="analysis"></a>

#### Running the Workflow

##### Inspect the execution plan

Check the core workflow DAG first:

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" snakemake all \
  --directory "$RUN_ROOT" \
  --snakefile "$CODE_ROOT/workflow/Snakefile" \
  --configfile "$RUN_ROOT/config/config.yaml" \
  --cores 12 \
  --dry-run \
  --printshellcmds
```

##### Run the core workflow only

Generate the final Pan-cluster table and refined count matrix, without sequence extraction or clustering:

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" snakemake all \
  --directory "$RUN_ROOT" \
  --snakefile "$CODE_ROOT/workflow/Snakefile" \
  --configfile "$RUN_ROOT/config/config.yaml" \
  --cores 12 \
  --rerun-incomplete \
  --printshellcmds
```

##### Run from inputs to all downstream plots

Replace `snakemake all` above with `snakemake downstream_all`, keeping the other arguments unchanged. This runs from the original inputs through sequence extraction, clustering, PCoA, and structure plotting for all Pan-clusters.

When using a job scheduler, `--cores` must not exceed the CPUs allocated to the job.

<a id="downstream"></a>

### Downstream Analysis

After the core workflow finishes, select Pan-clusters of interest without editing `config.yaml`.

<a id="extraction"></a>

#### Sequence Extraction

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag extract \
  --config "$RUN_ROOT/config/config.yaml" \
  --tag Pan_cluster_0132 \
  --tag Pan_cluster_0433 \
  --cores 4
```

<a id="clustering"></a>

#### Sequence Clustering and Visualization

Group genes within each Pan-cluster by sequence similarity and generate PCoA, per-sample structure plots, and collapsed structure plots. No separate plotting command is needed.

Pan-TAG uses the silhouette score (ranging from −1 to 1) to evaluate within-group cohesion and between-group separation of sequence groups within each Pan-cluster. Higher scores generally indicate clearer separation. We use **silhouette score ≥ 0.90** as an empirical reference for well-separated groups, not as a threshold for determining correctness. Lower scores do not necessarily indicate incorrect grouping and should be examined alongside PCoA and gene clustering results.

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag cluster \
  --config "$RUN_ROOT/config/config.yaml" \
  --tag Pan_cluster_0132 \
  --tag Pan_cluster_0433 \
  --cores 4
```

Replace the repeated `--tag` options with `--all` to process all Pan-clusters. Both commands accept `--dry-run --printshellcmds` to preview the execution plan.

For example, cluster and visualize all Pan-clusters from an existing core run:

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag cluster \
  --config "$RUN_ROOT/config/config.yaml" --all --cores 12
```

`pantag cluster` requires existing core outputs, including the internal manifest and FASTA indexes in `work_dir`. It can generate missing sequences but does not regenerate missing core outputs. For a given selection, all selected clusters finish sequence clustering before structure plotting begins. Both structure plots are generated only for clusters with `status=ok` and valid PCoA coordinates. Clusters with too few sequences or no usable k-mers retain a status record but receive no structure plots.

<a id="parameters"></a>

## Key Parameters

Parameters are stored in the generated `config.yaml`:

| Configuration key | Description | Template value |
|---|---|---:|
| `liftoff.threads` | Threads per Liftoff job | 8 |
| `tag.diamond_threads` | Threads per DIAMOND blastp job | 4 |
| `tag.diamond_max_target_seqs` | Maximum reported non-self targets per query | 5 |
| `tag.blastp_*` | blastp identity, query coverage, and subject coverage thresholds (%) | 20 |
| `pseudogene.search_flank_bp` | Flanking search distance for pseudogene recovery (bp) | 50000 |
| `pseudogene.blastx_*` | blastx identity, query coverage, and subject coverage thresholds (%) | 20 |
| `flanks.gene_count` | Genes on each TAG flank used for QC/synteny | 10 |
| `merge.member_overlap_threshold` | Cross-sample member overlap threshold | 0.3 |
| `merge.require_same_seqid` | Restrict merging candidates to matching seqid names | false |
| `refine.flank_overlap_threshold` | Flanking-gene overlap threshold for refinement | 0.3 |
| `downstream.cluster.k` | Canonical k-mer length | 31 |
| `downstream.cluster.ubiquitous_threshold` | Remove k-mers occurring in a fraction of sequences greater than this value | 0.95 |
| `downstream.cluster.max_clusters` | Maximum number of groups tested; 0 uses the number of distinct gene names + 3, capped by sequence count | 0 |
| `downstream.cluster.silhouette_tolerance` | Allowed silhouette-score loss when preferring fewer groups | 0.05 |
| `downstream.cluster.min_sequences` | Minimum sequences required for clustering | 10 |
| `downstream.cluster.pcoa` | Generate PCoA plots | true |

DIAMOND blastp applies both query and subject coverage thresholds and excludes self-hits with `--no-self-hits`. For each sample, MCScanX inputs exclude nested genes fully contained within the preceding retained gene interval. Their count is recorded as `nested_genes_skipped` in `<sample>_prepare_qc.tsv`.

<a id="outputs"></a>

## Outputs

### Results directory

The example below assumes `results_dir: results`. Downstream files are generated only when the corresponding analysis is requested.

```text
results/
├── local_tags/
│   └── <sample_id>_TAGs.tsv
├── global/
│   ├── pangenome_TAGs.tsv
│   └── refined_TAG_count_matrix.tsv
└── downstream/
    ├── global_cluster_summary.tsv
    └── Pan_cluster_XXXX/
        ├── sequences/
        │   └── Pan_cluster_XXXX_sequences.fa
        └── clusters/
            ├── Pan_cluster_XXXX_pcoa.pdf
            ├── Pan_cluster_XXXX_pcoa_coordinates.tsv
            ├── Pan_cluster_XXXX_TAG_structure.pdf
            └── Pan_cluster_XXXX_TAG_structure_collapsed.pdf
```

- `<sample_id>_TAGs.tsv`: cleaned per-sample TAGs, including members, coordinates, density, and flanking genes; pseudogene columns are included when pseudogenes are recovered.
- `pangenome_TAGs.tsv`: final Pan-cluster table with member names/IDs, source sample, seqid, sample occupancy, and merged left/right flanking genes. `pan_cluster_size` counts distinct member names, not total gene copies across samples. `gene_names` and `gene_ids` are independently deduplicated and sorted sets, not positionally paired lists.
- `refined_TAG_count_matrix.tsv`: Pan-cluster × sample count matrix, with detection counts, gap-free sample counts, mean, standard deviation, and coefficient of variation.
- `Pan_cluster_XXXX_sequences.fa`: sequences from gap-free, non-span records in the Pan-cluster; headers include sample, gene name, and gene ID.
- `Pan_cluster_XXXX_pcoa.pdf`: two panels sharing the same coordinates, colored by original gene name on the left and sequence group (`Pan_cluster_XXXX-SGXX`) on the right.
- `Pan_cluster_XXXX_pcoa_coordinates.tsv`: per-sequence PCoA coordinates, gene information, sequence group, and assigned color.
- `Pan_cluster_XXXX_TAG_structure.pdf`: TAG structural haplotypes shown for individual samples. Target genes are colored by **Sequence Group inferred by sequence-based reclustering**, rather than by gene names in the original annotation; flanking genes are shown in gray.
- `Pan_cluster_XXXX_TAG_structure_collapsed.pdf`: structural haplotypes collapsed by copy number and gene arrangement, annotated with the corresponding sample counts. Arrows carry `SGxx` labels, and flanking genes are omitted. Structural haplotype details for each sample are recorded in `work/downstream/clusters/<global_tag_id>/TAG_compositions.tsv`.
- `global_cluster_summary.tsv`: cumulative statuses and key metrics for completed clustering analyses.

### Count matrix

Sample cells use the notation below. In `N*` and `N#`, `N` denotes the count (e.g., `3*` or `3#`):

| Value | Meaning |
|---|---|
| Positive integer | One valid block, with no ambiguous N bases in the target interval |
| `N*` | One valid block, with ambiguous N bases in the target interval |
| `N#` | Multiple blocks supported by flanking genes (span) |
| `NA` | No target genes found, or no block passed flanking-gene QC |

Refinement removes a candidate cluster only when its fraction of span samples is strictly greater than 0.9.

### Intermediate files and logs

With `outputs.debug: true`, the following file is also published:

```text
results/global/local_to_global.tsv
```

`TAG_sequence_manifest.tsv` is the internal manifest used for downstream extraction and clustering. It records member coordinates, QC status, and extraction eligibility; only eligible records are extracted. Retain `work_dir` if further extraction, clustering, or plotting is planned.

<a id="troubleshooting"></a>

## Resuming Runs and Troubleshooting

### Resume an interrupted run

Rerun the original command after a failure. With `--rerun-incomplete`, Snakemake regenerates incomplete outputs and reuses outputs that remain up to date.

If an interrupted run leaves a lock, first confirm that no other Snakemake process is using the run directory, then execute:

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" snakemake \
  --directory "$RUN_ROOT" \
  --snakefile "$CODE_ROOT/workflow/Snakefile" \
  --configfile "$RUN_ROOT/config/config.yaml" \
  --unlock
```

### Run multiple datasets concurrently

Datasets can share `CODE_ROOT` and the Pixi environment, but must use separate `RUN_ROOT`, `results_dir`, `logs_dir`, and `work_dir` paths. Never run two jobs against the same `work_dir`.

### Common messages

- `Nothing to be done`: requested outputs already exist and are up to date; this is not a failure.
- `IncompleteFilesException`: rerun the original command with `--rerun-incomplete`; do not mark unverified files as complete.
- `Directory cannot be locked`: another process is running, or a previous interrupted run left a lock.
- More jobs after editing `config.yaml`: inspect a `--dry-run`; changed inputs or parameters may trigger affected steps to rerun.

<a id="citation-license"></a>

## Citation and License

### Citation

Pan-TAG v0.1 is under development. Publication and citation details will be added when available.

### License

The repository does not yet include a project-level `LICENSE`. The upstream license for the customized MCScanX is provided in `src/pantag/mcscanx/LICENSE`.

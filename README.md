# ChIPseq_Processing

A SLURM ChIP-seq pipeline written in plain bash. Every step is its own `sbatch`
job, and the jobs are chained with `--dependency=afterok`, so a single command
submits the whole analysis for every sample in a samplesheet. It follows the
design of [nf-core/chipseq](https://nf-co.re/chipseq): same samplesheet, same
steps and a similar output layout. Filtering, QC metrics, signal tracks and
replicate reproducibility follow the
[ENCODE ChIP-seq pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline2)
and the [ENCODE standards](https://www.encodeproject.org/chip-seq/).
Defaults are set for the UVM VACC (`bluemoon` partition, hg38).

## Workflow

```
 per FASTQ (run)     trim ──► align
                     FastQC, Trim Galore     bowtie2 → sorted BAM
                                 │
 per library         filter ──► bamqc
 (sample_REPn)       merge lanes, ENCODE     phantompeakqualtools (NSC, RSC,
                     filters, MarkDuplicates, fragment length), bamCoverage
                     NRF/PBC1/PBC2, dedup     bigWig
                                 │   (+ the control library's filter job)
 per IP library      peaks
                     MACS main call (narrow/broad) + FRiP, relaxed call,
                     fold-enrichment & p-value bigWigs, plotFingerprint
                          │                                  │
 per IP sample       repro                     per antibody  consensus
                     IDR (narrow) / overlap    merged peak set +
                     (broad), pseudoreplicates, featureCounts matrix
                     conservative & optimal sets
                                 │
 once                multiqc (runs after all other jobs, even failed ones)
```

| Step | Script | What it does |
|------|--------|--------------|
| trim | `bin/01_trim.sh` | FastQC on the raw reads; Trim Galore (adapter + quality trimming) with FastQC on the trimmed reads |
| align | `bin/02_align.sh` | bowtie2 (`-X 2000` for paired-end), coordinate sort |
| filter | `bin/03_filter.sh` | Merges the runs of a library. ENCODE filtering: `-F 1804` (plus `-f 2` and `fixmate` for PE), MAPQ ≥ 30, Picard MarkDuplicates, library complexity (NRF, PBC1, PBC2), duplicate removal. Writes flagstat, idxstats and samtools stats. |
| bamqc | `bin/04_bam_qc.sh` | phantompeakqualtools on 15M subsampled reads (NSC, RSC, fragment length); CPM-normalised bigWig |
| peaks | `bin/05_call_peaks.sh` | Main MACS call: narrow uses q < 0.05, broad uses `--broad-cutoff 0.1`, as in nf-core. Blacklist filter, FRiP. Relaxed ENCODE call (p < 0.01, top 500k peaks). Fold-enrichment and −log10(p) signal bigWigs. deepTools fingerprint. |
| repro | `bin/06_reproducibility.sh` | ENCODE reproducibility: true replicates, pooled pseudoreplicates and self-pseudoreplicates. Comparisons use IDR (narrow) or naive overlap (broad). Gives conservative and optimal peak sets, the rescue and self-consistency ratios, and a pass/borderline/fail flag. |
| consensus | `bin/07_consensus.sh` | Merged peak set across all libraries of an antibody, a presence/absence matrix, and featureCounts read counts for differential binding |
| multiqc | `bin/08_multiqc.sh` | Summary tables of every metric, plus the MultiQC report |

Paired-end libraries are peak-called as `BAMPE`. For single-end libraries the
reads are shifted by the SPP fragment length estimate
(`--nomodel --extsize <fraglen>`, as ENCODE does), falling back to
`DEFAULT_FRAGLEN`.

## Setup (once)

1. **Software.** Create the conda environment (mamba or micromamba also work):
   ```bash
   conda env create -f environment.yml
   ```
   Then edit `setup_env()` in `conf/chipseq.config` so that each job
   activates it:
   ```bash
   setup_env() {
       source "$HOME/miniconda3/etc/profile.d/conda.sh"
       conda activate chipseq
   }
   ```
   You can use `module load ...` lines here instead if the tools come from
   modules.

   The environment installs **MACS3**, not MACS2. The current bioconda builds of
   MACS2 2.2.9.1 crash on newer Linux (glibc ≥ 2.31) with
   `undefined symbol: __log_finite`. MACS3 takes the same arguments and runs
   the same algorithm. To use a working MACS2 (for example a VACC module), set
   `MACS_CMD="macs2"`.

   IDR 2.0.4.2 calls `numpy.int`, which NumPy 1.24 removed, and MACS3 and
   deepTools need a newer NumPy. The pipeline therefore runs IDR through a
   small shim that restores the old aliases. You don't need to do anything for
   this. If you have IDR in its own environment, set `IDR_CMD` to that `idr`
   binary instead.

2. **References.** Check the paths in `conf/chipseq.config`:
   - `BOWTIE2_INDEX`: the bowtie2 index prefix. The default is your
     GRCh38_noalt_as index.
   - `CHROM_SIZES`: a `chrom<TAB>length` file for the same assembly.
   - `BLACKLIST`: the ENCODE blacklist, for example
     [hg38-blacklist.v2.bed](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz).
     Download and unzip it to the configured path, or set it to `""`.
   - `MACS_GSIZE` and `EFFECTIVE_GENOME_SIZE` for other species.

3. **SLURM.** Set the partition, account and per-step resources in the config.
   `SLURM_EXTRA_ARGS=(--mail-user=you@uvm.edu --mail-type=FAIL)` sends an email
   when a job fails.

## Samplesheet

Use the nf-core/chipseq (v2) format: a CSV file with a header row, where each
row is one FASTQ file or FASTQ pair. See
[`assets/samplesheet_example.csv`](assets/samplesheet_example.csv).

```csv
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate,peak_type
CTCF_WT,/path/CTCF_WT_R1_L001_1.fastq.gz,/path/CTCF_WT_R1_L001_2.fastq.gz,1,CTCF,INPUT_WT,1,narrow
CTCF_WT,/path/CTCF_WT_R1_L002_1.fastq.gz,/path/CTCF_WT_R1_L002_2.fastq.gz,1,CTCF,INPUT_WT,1,narrow
CTCF_WT,/path/CTCF_WT_R2_1.fastq.gz,/path/CTCF_WT_R2_2.fastq.gz,2,CTCF,INPUT_WT,2,narrow
H3K27me3_WT,/path/H3K27me3_WT_R1_1.fastq.gz,/path/H3K27me3_WT_R1_2.fastq.gz,1,H3K27me3,INPUT_WT,1,broad
INPUT_WT,/path/INPUT_WT_R1_1.fastq.gz,/path/INPUT_WT_R1_2.fastq.gz,1,,,,
INPUT_WT,/path/INPUT_WT_R2_1.fastq.gz,/path/INPUT_WT_R2_2.fastq.gz,2,,,,
```

| Column | Description |
|--------|-------------|
| `sample` | Sample name (letters, numbers, `.`, `_`, `-`). All replicates of a condition share it. |
| `fastq_1`, `fastq_2` | Gzipped FASTQ files. Leave `fastq_2` empty for single-end data. Relative paths are resolved from the directory you submit from. |
| `replicate` | Biological replicate number (1, 2, ...). Rows with the same sample and replicate are merged as technical replicates or lanes. |
| `antibody` | Target. Leave empty for input/control samples. Libraries that share an antibody form one consensus peak set. |
| `control` | The `sample` name of the input for this IP. Leave empty to call peaks without a control. |
| `control_replicate` | Which replicate of `control` to use. If empty, the same replicate number is used when it exists, otherwise replicate 1. |
| `peak_type` | Optional, `narrow` or `broad`. Default: `PEAK_TYPE_DEFAULT`. Use narrow for TFs, H3K4me3, H3K27ac and H3K9ac. Use broad for H3K27me3, H3K36me3, H3K9me3, H3K79me2 and H3K4me1. |

Libraries are named `<sample>_REP<replicate>`, and runs are named
`<library>_T<n>`. The pipeline validates the whole sheet before it submits
anything. It checks for missing files, duplicated FASTQs, mixed SE/PE runs
within a library, unknown controls and inconsistent antibodies.

## Running

```bash
# 1. validate and preview the job graph (submits nothing)
./submit_chipseq.sh -i samplesheet.csv -o /path/to/results --dry-run

# 2. submit everything
./submit_chipseq.sh -i samplesheet.csv -o /path/to/results

# use a project-specific config
./submit_chipseq.sh -i samplesheet.csv -o results -c my_project.config

# 3. monitor
squeue -u $USER
./pipeline_status.sh /path/to/results
```

Run the submit script from a login node; it only calls `sbatch`. Jobs that
depend on a failed job are cancelled automatically
(`--kill-on-invalid-dep=yes`). The MultiQC job still runs afterwards, so you
get a report of whatever finished.

**Resume after a failure.** Fix the problem (for example raise `ALIGN_TIME`),
then resubmit with `-r`:
```bash
./submit_chipseq.sh -i samplesheet.csv -o /path/to/results -r
```
Any step with a completion marker in `results/pipeline_info/done/` is
skipped, and everything else is resubmitted with the correct dependencies. To
force a step to rerun, delete its `.done` file. With `-r` you can also add new
samples to the samplesheet and only the new work runs. Reproducibility,
consensus and MultiQC must include the new samples, so delete their `.done`
files (`repro.*`, `consensus.*`) before resubmitting.

Each submission saves a copy of the samplesheet and config to
`pipeline_info/`, and the jobs read that copy, so editing the config never
affects a run already in the queue. `pipeline_info/jobs_<timestamp>.tsv`
lists every job ID, step, dependencies and log file.

## Output

```
results/
├── fastqc/raw/                        FastQC of the raw reads
├── trimgalore/                        trimming reports (+ fastqc/ of trimmed reads)
├── bowtie2/
│   ├── logs/                          bowtie2 alignment summaries per run
│   └── merged_library/
│       ├── <lib>.filt.dedup.bam(.bai) final filtered, de-duplicated BAMs
│       └── qc/                        flagstat, idxstats, stats, Picard metrics, <lib>.pbc.qc
├── phantompeakqualtools/              <lib>.spp.out (NSC/RSC), cross-correlation plots, fraglen
├── bigwig/                            <lib>.CPM.bigWig coverage tracks
├── deeptools/plotFingerprint/
├── macs2/
│   ├── narrow/ | broad/
│   │   ├── <lib>_peaks.narrowPeak / .broadPeak         raw MACS output
│   │   ├── <lib>_peaks.bfilt.narrowPeak / .broadPeak   blacklist filtered
│   │   ├── qc/<lib>.FRiP.txt
│   │   └── consensus/<antibody>/      consensus peaks (.bed, .saf, .boolean.tsv)
│   │                                  and featureCounts matrix
│   ├── relaxed/<lib>/                 relaxed (p < 0.01) peaks used for IDR
│   └── signal/                        <lib>.fc_signal.bigWig, <lib>.pval_signal.bigWig
├── reproducibility/<sample>/
│   ├── <sample>.optimal.<idr|overlap>.narrowPeak        recommended final peak set
│   ├── <sample>.conservative.<idr|overlap>.narrowPeak
│   ├── <sample>.pooled.relaxed.narrowPeak
│   ├── <sample>.reproducibility.qc.tsv
│   └── comparisons/                   every IDR/overlap comparison (+ IDR plots)
├── qc/
│   ├── library_qc_summary.tsv         one row per library, all metrics
│   └── reproducibility_summary.tsv    one row per sample
├── multiqc/multiqc_report.html
├── logs/<step>/                       SLURM logs
├── pipeline_info/                     config/samplesheet copies, job lists, done markers
└── work/                              intermediates (safe to delete when finished)
```

**Which peaks to use.** For TFs and narrow marks with replicates, use
`reproducibility/<sample>/<sample>.optimal.idr.narrowPeak`. It is ENCODE's
"IDR thresholded optimal" set, with IDR < 0.05. For broad marks, use
`<sample>.optimal.overlap.narrowPeak` (the ENCODE "replicated peaks") or the
per-replicate `macs2/broad/*.bfilt.broadPeak`. For differential binding
between conditions, use the featureCounts matrix in `macs2/<type>/consensus/`
(for example with DESeq2 or DiffBind).

## QC guide (ENCODE)

All of these values are in `qc/library_qc_summary.tsv` and in the MultiQC
report.

| Metric | Meaning | Target |
|--------|---------|--------|
| `final_reads` | Usable reads after filtering and de-duplication | TF ≥ 20M fragments per replicate; narrow histone ≥ 20M; broad histone ≥ 45M |
| `NRF` | Distinct / total reads | > 0.9 ideal, 0.8–0.9 acceptable, < 0.5 concerning |
| `PBC1` | One-read positions / distinct positions | > 0.9 none, 0.8–0.9 mild, 0.5–0.8 moderate, < 0.5 severe bottlenecking |
| `PBC2` | One-read / two-read positions | > 10 none, 3–10 mild, 1–3 moderate, < 1 severe |
| `NSC` | Normalised strand cross-correlation | > 1.05 |
| `RSC` | Relative strand cross-correlation | > 0.8 (≥ 1 good) |
| `FRiP` | Fraction of reads in (main) peaks | TF > 0.01 (ENCODE), usually much higher for good histone data |
| rescue ratio, self-consistency ratio | Replicate agreement | ≤ 2 (both ≤ 2 = pass, one > 2 = borderline, both > 2 = fail) |

## Configuration reference

Every option is documented in `conf/chipseq.config`. The most useful ones:

| Option | Default | Notes |
|--------|---------|-------|
| `PEAK_TYPE_DEFAULT` | `narrow` | Used when the samplesheet has no `peak_type` |
| `MACS_CMD` | `macs3` | or `macs2` |
| `MACS_QVAL` / `MACS_BROAD_CUTOFF` | 0.05 / 0.1 | Main peak calls (nf-core defaults) |
| `RELAXED_PVAL` / `RELAXED_NPEAKS` | 0.01 / 500000 | Relaxed calls for IDR (ENCODE defaults) |
| `IDR_THRESHOLD` | 0.05 | |
| `REPRO_METHOD` | `auto` | `idr`, `overlap`, or `auto` (IDR for narrow, overlap for broad) |
| `MAPQ_THRESHOLD` | 30 | |
| `KEEP_DUPLICATES` | `false` | |
| `FILTER_CHRS` | `""` | e.g. `"chrM"` |
| `RUN_SPP`, `RUN_FINGERPRINT`, `RUN_REPRODUCIBILITY`, `RUN_CONSENSUS` | `true` | Turn optional steps off |
| `CLEAN_INTERMEDIATES` | `true` | Delete trimmed FASTQs and per-run BAMs once merged |

## Differences from nf-core/chipseq

- nf-core blacklist-filters the BAM files. This pipeline follows ENCODE and
  blacklist-filters the peaks instead.
- nf-core has no IDR or pseudoreplicate analysis; this pipeline adds ENCODE's
  replicate reproducibility step.
- Not included: HOMER peak annotation, deepTools gene-body profiles, preseq and
  the DESeq2 report. The consensus featureCounts matrix is ready for DESeq2 or
  DiffBind.

#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -J hybrid_metagenomics
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --output=hibrido_process_%j.log

# ==========================================================================
# Hybrid metagenomics pipeline for spontaneous/mixed fermentation systems
# Genome-resolved recovery of BOTH bacterial and fungal MAGs, functional
# annotation, and ANI-based fungal species identification.
#
# Paired-end (Illumina) version. See pipeline_hibrido_singleend.sh for
# single-end (Ion Torrent) data.
# ==========================================================================

set -uo pipefail

# --------------------------------------------------------------------------
# 0. CONFIGURATION — edit all paths below for your environment
# --------------------------------------------------------------------------
MICROMAMBA_BIN_PEDRO="/path/to/micromamba"

# Tool environments (each isolated; activated in a subshell per call)
ENV_TIARA="/path/to/env_tiara"
ENV_METAEUK="/path/to/env_metaeuk"
ENV_BUSCO_PATH="/path/to/env_busco"
ENV_KOFAMSCAN="/path/to/env_kofamscan"
ENV_BRACKEN="/path/to/env_bracken"
ENV_FASTANI="/path/to/env_fastani"
PROKARYOTE_ENV="/path/to/metawrap-env"   # MetaWRAP, Prodigal, PROKKA, bwa, samtools, CheckM

# Binaries / reference databases
METAEUK_BIN="$ENV_METAEUK/bin/metaeuk"
METAEUK_REF_DB="/path/to/banco_proteinas/metaeuk_swissprotDB"
BUSCO_BIN="$ENV_BUSCO_PATH/bin/busco"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-fungi_odb12.2}"
BUSCO_DOWNLOAD_PATH="/path/to/busco_downloads/busco_downloads"
KOFAM_PROFILE="/path/to/kofam_db/profiles"
KOFAM_KO_LIST="/path/to/kofam_db/ko_list"
KRAKEN2_DB="${KRAKEN2_DB:-/path/to/banco_kraken}"
BRACKEN_READ_LEN="${BRACKEN_READ_LEN:-150}"          # match your read length; DB has 50/75/100/150/200/250/300
TRIMMOMATIC_JAR="/path/to/trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS_DIR="/path/to/trimmomatic/adapters"
HUMAN_REF_PREFIX=""                                   # bmtagger index prefix; leave empty to skip human decontamination

# ANI reference genomes for fungal species ID (Fase 5b). Add more entries
# as needed for your system; accessions below cover common fermentation
# yeasts (Saccharomyces cerevisiae, Brettanomyces bruxellensis,
# Zygosaccharomyces rouxii).
ANI_REF_DIR="/path/to/genomas_referencia_leveduras/candidatos_ani"
declare -A ANI_REFS
ANI_REFS["Saccharomyces_cerevisiae"]="GCF_000146045.2"
ANI_REFS["Brettanomyces_bruxellensis"]="GCF_011074885.1"
ANI_REFS["Zygosaccharomyces_rouxii"]="GCA_000026365.1"
ANI_MIN_IDENTITY=90   # below this, report as "no confident match" rather than best-of-bad-options

BASE_DIR="/path/to/project"
THREADS=20
MIN_CONTIG_LEN=1500

SRA_IDS=("${SAMPLE_ID:?Set SAMPLE_ID via --export=ALL,SAMPLE_ID=<accession>}")

RAW_DIR="$BASE_DIR/raw_reads"
QC_DIR="$BASE_DIR/qc_reads"
DECONT_DIR="$BASE_DIR/decontaminated_reads"
KRAKEN_DIR="$BASE_DIR/kraken_bracken"
ASSEMBLY_DIR="$BASE_DIR/assembly"
BIFURCATION_DIR="$BASE_DIR/bifurcation"
BIN_PRO_DIR="$BASE_DIR/binning_prokariotos"
BIN_EUK_DIR="$BASE_DIR/binning_eukariotos"
FUNC_DIR="$BASE_DIR/anotacao_funcional"
PROSPEC_DIR="$BASE_DIR/prospeccao_especifica"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_DIR" "$QC_DIR" "$DECONT_DIR" "$KRAKEN_DIR" "$ASSEMBLY_DIR" "$BIFURCATION_DIR" \
         "$BIN_PRO_DIR" "$BIN_EUK_DIR" "$FUNC_DIR" "$PROSPEC_DIR" "$LOG_DIR" "$ANI_REF_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

# check_step: verify a step produced non-empty output before continuing.
# Prints OK/ERROR and returns non-zero on failure so callers can `|| continue`.
check_step() {
    local step_name="$1" output_path="$2"
    if [ -s "$output_path" ]; then
        log "OK: step '$step_name' completed."
        return 0
    else
        log "ERROR: step '$step_name' produced no output ($output_path). Skipping sample."
        return 1
    fi
}

# --------------------------------------------------------------------------
# Fetch ANI reference genomes once (skips if already present)
# --------------------------------------------------------------------------
fetch_ani_references() {
    for name in "${!ANI_REFS[@]}"; do
        local acc="${ANI_REFS[$name]}"
        local fna
        fna=$(find "$ANI_REF_DIR" -path "*${acc}*" -name "*.fna" 2>/dev/null | head -1)
        if [ -z "$fna" ]; then
            log "Downloading ANI reference $name ($acc)..."
            (
                cd "$ANI_REF_DIR"
                datasets download genome accession "$acc" --include genome --filename "${acc}.zip" \
                    && unzip -o "${acc}.zip" -d "$acc" >/dev/null
            )
        fi
    done
}

# ==========================================================================
# MAIN LOOP
# ==========================================================================
for SRA_ID in "${SRA_IDS[@]}"; do
    SAMPLE_LOG="$LOG_DIR/${SRA_ID}.log"
    log "=== Processing $SRA_ID ===" | tee -a "$SAMPLE_LOG"

    # ---- Fase 1: QC (Trimmomatic) + optional human decontamination ----
    READS_R1="$QC_DIR/${SRA_ID}_R1_paired.fq"
    READS_R2="$QC_DIR/${SRA_ID}_R2_paired.fq"

    if [ -s "$READS_R1" ] && [ -s "$READS_R2" ]; then
        log "[Fase 1] Existing clean reads found, skipping Trimmomatic." | tee -a "$SAMPLE_LOG"
    else
        log "[Fase 1] Trimmomatic..." | tee -a "$SAMPLE_LOG"
        java -jar "$TRIMMOMATIC_JAR" PE -threads $THREADS \
            "$RAW_DIR/${SRA_ID}_1.fastq" "$RAW_DIR/${SRA_ID}_2.fastq" \
            "$QC_DIR/${SRA_ID}_R1_paired.fq" "$QC_DIR/${SRA_ID}_R1_unpaired.fq" \
            "$QC_DIR/${SRA_ID}_R2_paired.fq" "$QC_DIR/${SRA_ID}_R2_unpaired.fq" \
            ILLUMINACLIP:"$TRIMMOMATIC_ADAPTERS_DIR"/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
            2>&1 | tee -a "$SAMPLE_LOG"
        check_step "Trimmomatic" "$READS_R1" || continue

        if [ -n "$HUMAN_REF_PREFIX" ]; then
            log "[Fase 1b] bmtagger human read removal..." | tee -a "$SAMPLE_LOG"
            bmtagger.sh -b "${HUMAN_REF_PREFIX}.bitmask" -x "${HUMAN_REF_PREFIX}.srprism" \
                -T /tmp -q1 -1 "$READS_R1" -2 "$READS_R2" -o "$DECONT_DIR/${SRA_ID}" 2>&1 | tee -a "$SAMPLE_LOG"
            # bmtagger.sh output naming varies by version; match any *1*/*2* file.
            shopt -s nullglob
            CANDIDATES_R1=("$DECONT_DIR/${SRA_ID}"*1*.fastq)
            CANDIDATES_R2=("$DECONT_DIR/${SRA_ID}"*2*.fastq)
            shopt -u nullglob
            if [ ${#CANDIDATES_R1[@]} -ge 1 ] && [ ${#CANDIDATES_R2[@]} -ge 1 ]; then
                READS_R1="${CANDIDATES_R1[0]}"; READS_R2="${CANDIDATES_R2[0]}"
            else
                log "WARNING: bmtagger output not recognized, continuing without decontamination." | tee -a "$SAMPLE_LOG"
            fi
        fi
    fi

    # ---- Fase 1c: direct taxonomic profile (Kraken2 + Bracken) ----
    log "[Fase 1c] Kraken2+Bracken..." | tee -a "$SAMPLE_LOG"
    mkdir -p "$KRAKEN_DIR/$SRA_ID"
    if [ -d "$KRAKEN2_DB" ]; then
        (
            set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$ENV_BRACKEN"; set -u
            kraken2 --db "$KRAKEN2_DB" --threads $THREADS --paired --memory-mapping \
                --report "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_kraken2_report.txt" \
                --output "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_kraken2_output.txt" \
                "$READS_R1" "$READS_R2"
            bracken -d "$KRAKEN2_DB" -i "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_kraken2_report.txt" \
                -o "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_bracken_species.txt" \
                -w "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_bracken_report.txt" \
                -r "$BRACKEN_READ_LEN" -l S -t 10
        ) 2>&1 | tee -a "$SAMPLE_LOG"
    else
        log "WARNING: KRAKEN2_DB not found, skipping Fase 1c." | tee -a "$SAMPLE_LOG"
    fi

    # ---- Fase 2: assembly (MEGAHIT) ----
    log "[Fase 2] MEGAHIT..." | tee -a "$SAMPLE_LOG"
    rm -rf "$ASSEMBLY_DIR/$SRA_ID"
    megahit -1 "$READS_R1" -2 "$READS_R2" -t $THREADS --min-contig-len $MIN_CONTIG_LEN \
        -o "$ASSEMBLY_DIR/$SRA_ID" 2>&1 | tee -a "$SAMPLE_LOG"
    check_step "MEGAHIT" "$ASSEMBLY_DIR/$SRA_ID/final.contigs.fa" || continue

    # ---- Fase 3: domain bifurcation (Tiara) ----
    log "[Fase 3] Tiara..." | tee -a "$SAMPLE_LOG"
    mkdir -p "$BIFURCATION_DIR/$SRA_ID"
    (
        cd "$BIFURCATION_DIR/$SRA_ID" || exit 1
        set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$ENV_TIARA"; set -u
        tiara -i "$ASSEMBLY_DIR/$SRA_ID/final.contigs.fa" -o "tiara_out.txt" -t $THREADS --tf all
    ) 2>&1 | tee -a "$SAMPLE_LOG"

    PROK_FASTA=$(find "$BIFURCATION_DIR/$SRA_ID" -iname "bacteria*" | head -1)
    EUK_FASTA=$(find "$BIFURCATION_DIR/$SRA_ID" -iname "eukarya*" | head -1)
    check_step "Tiara (prokaryotic contigs)" "$PROK_FASTA" || continue

    # ---- Fase 4A: bacterial binning (MetaWRAP) ----
    log "[Fase 4A] Bacterial MAG recovery..." | tee -a "$SAMPLE_LOG"
    set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$PROKARYOTE_ENV"; set -u

    METAWRAP_READS_DIR="$BIN_PRO_DIR/$SRA_ID/_reads_link"
    mkdir -p "$METAWRAP_READS_DIR"
    ln -sf "$(readlink -f "$READS_R1")" "$METAWRAP_READS_DIR/${SRA_ID}_1.fastq"
    ln -sf "$(readlink -f "$READS_R2")" "$METAWRAP_READS_DIR/${SRA_ID}_2.fastq"

    metawrap binning -o "$BIN_PRO_DIR/$SRA_ID" -t $THREADS -a "$PROK_FASTA" \
        --metabat2 --maxbin2 --concoct \
        "$METAWRAP_READS_DIR/${SRA_ID}_1.fastq" "$METAWRAP_READS_DIR/${SRA_ID}_2.fastq" 2>&1 | tee -a "$SAMPLE_LOG"
    check_step "MetaWRAP binning" "$BIN_PRO_DIR/$SRA_ID/metabat2_bins" || continue

    BAC_BINS_DIR="$BIN_PRO_DIR/$SRA_ID/${SRA_ID}_refinado"
    metawrap bin_refinement -o "$BAC_BINS_DIR" -t $THREADS -c 90 -x 5 \
        -A "$BIN_PRO_DIR/$SRA_ID/metabat2_bins" \
        -B "$BIN_PRO_DIR/$SRA_ID/maxbin2_bins" \
        -C "$BIN_PRO_DIR/$SRA_ID/concoct_bins" 2>&1 | tee -a "$SAMPLE_LOG"

    if [ -s "$BAC_BINS_DIR/metawrap_90_5_bins.stats" ]; then
        log "OK: bin_refinement succeeded." | tee -a "$SAMPLE_LOG"
        BAC_BINS_FINAL_DIR="$BAC_BINS_DIR/metawrap_90_5_bins"
    else
        # metawrap bin_refinement has a known bug (consolidate_two_sets_of_bins.py);
        # fall back to CheckM directly on the metabat2 bins.
        log "WARNING: bin_refinement failed (known MetaWRAP bug). Falling back to CheckM directly." | tee -a "$SAMPLE_LOG"
        CHECKM_DIR="$BIN_PRO_DIR/$SRA_ID/${SRA_ID}_checkm_direto"
        mkdir -p "$CHECKM_DIR/bins_filtrados_corrigido"
        checkm lineage_wf -t $THREADS -x fa "$BIN_PRO_DIR/$SRA_ID/metabat2_bins" "$CHECKM_DIR" \
            --reduced_tree --tab_table -f "$CHECKM_DIR/checkm_results.txt" 2>&1 | tee -a "$SAMPLE_LOG"

        # Columns 12 (Completeness) / 13 (Contamination) in CheckM's tab_table output.
        awk -F'\t' -v comp=90 -v cont=5 \
            'NR>1 && $12+0 >= comp && $13+0 <= cont {print $1".fa"}' \
            "$CHECKM_DIR/checkm_results.txt" > "$CHECKM_DIR/bins_aprovados.txt"

        while read -r binfile; do
            [ "$binfile" = "bin.unbinned.fa" ] && continue   # not a real MAG
            cp "$BIN_PRO_DIR/$SRA_ID/metabat2_bins/$binfile" "$CHECKM_DIR/bins_filtrados_corrigido/" 2>/dev/null
        done < "$CHECKM_DIR/bins_aprovados.txt"

        BAC_BINS_FINAL_DIR="$CHECKM_DIR/bins_filtrados_corrigido"
        if [ -z "$(ls -A "$BAC_BINS_FINAL_DIR" 2>/dev/null)" ]; then
            log "WARNING: no bacterial bin met the quality threshold for $SRA_ID." | tee -a "$SAMPLE_LOG"
        fi
    fi

    # ---- Fase 4B: fungal binning (MetaBAT2) ----
    EUK_BINS_DIR=""
    if [ -s "$EUK_FASTA" ]; then
        log "[Fase 4B] Fungal MAG recovery..." | tee -a "$SAMPLE_LOG"
        mkdir -p "$BIN_EUK_DIR/$SRA_ID/bins"
        bwa index "$EUK_FASTA" 2>&1 | tee -a "$SAMPLE_LOG"
        bwa mem -t $THREADS "$EUK_FASTA" "$READS_R1" "$READS_R2" 2>>"$SAMPLE_LOG" | \
            samtools sort -@ $THREADS -o "$BIN_EUK_DIR/$SRA_ID/euk_alignment.bam"
        samtools index "$BIN_EUK_DIR/$SRA_ID/euk_alignment.bam"
        jgi_summarize_bam_contig_depths --outputDepth "$BIN_EUK_DIR/$SRA_ID/depth.txt" \
            "$BIN_EUK_DIR/$SRA_ID/euk_alignment.bam" 2>&1 | tee -a "$SAMPLE_LOG"
        metabat2 -i "$EUK_FASTA" -a "$BIN_EUK_DIR/$SRA_ID/depth.txt" \
            -o "$BIN_EUK_DIR/$SRA_ID/bins/euk_bin" -m 1500 2>&1 | tee -a "$SAMPLE_LOG"

        if [ -n "$(ls -A "$BIN_EUK_DIR/$SRA_ID/bins" 2>/dev/null)" ]; then
            EUK_BINS_DIR="$BIN_EUK_DIR/$SRA_ID/bins"
            log "[Fase 4B] BUSCO validation..." | tee -a "$SAMPLE_LOG"
            (
                set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$ENV_BUSCO_PATH"; set -u
                for bin in "$EUK_BINS_DIR"/*.fa; do
                    "$BUSCO_BIN" -i "$bin" -o "$(basename "$bin" .fa)" \
                        --out_path "$BIN_EUK_DIR/${SRA_ID}_busco" -l "$BUSCO_LINEAGE" \
                        -m genome --offline --download_path "$BUSCO_DOWNLOAD_PATH" -f 2>&1
                done
            ) | tee -a "$SAMPLE_LOG"
        else
            log "WARNING: MetaBAT2 produced no fungal bins for $SRA_ID." | tee -a "$SAMPLE_LOG"
        fi
    fi

    # ---- Fase 5: functional annotation ----
    log "[Fase 5] Functional annotation..." | tee -a "$SAMPLE_LOG"
    for bin in "$BAC_BINS_FINAL_DIR"/*.fa; do
        [ -f "$bin" ] || continue
        bin_name=$(basename "$bin" .fa)
        prodigal -i "$bin" -a "$FUNC_DIR/${SRA_ID}/${bin_name}_bac.faa" -p single
        prokka --outdir "$FUNC_DIR/${SRA_ID}/${bin_name}_prokka" --prefix "$bin_name" --cpus $THREADS --force "$bin"
    done

    for bin in "$EUK_BINS_DIR"/*.fa; do
        [ -f "$bin" ] || continue
        bin_name=$(basename "$bin" .fa)
        "$METAEUK_BIN" easy-predict "$bin" "$METAEUK_REF_DB" \
            "$FUNC_DIR/${SRA_ID}/${bin_name}_metaeuk" "$FUNC_DIR/${SRA_ID}/${bin_name}_tmp" --threads $THREADS
        rm -rf "$FUNC_DIR/${SRA_ID}/${bin_name}_tmp"
    done

    (
        set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$ENV_KOFAMSCAN"; set -u
        for faa in "$FUNC_DIR/${SRA_ID}"/*_bac.faa "$FUNC_DIR/${SRA_ID}"/*_metaeuk.fas; do
            [ -f "$faa" ] || continue
            out_name=$(basename "$faa" | sed 's/\.faa$//; s/_metaeuk\.fas$/_euk/')
            exec_annotation -o "$FUNC_DIR/${SRA_ID}/${out_name}_kegg.txt" "$faa" \
                -p "$KOFAM_PROFILE" -k "$KOFAM_KO_LIST" --cpu $THREADS -f mapper
        done
    ) 2>&1 | tee -a "$SAMPLE_LOG"

    # ---- Fase 5b: fungal species ID via genome-level ANI ----
    # Protein-homology ID (previously used here) is unreliable for
    # fermentation yeasts under-represented in Swiss-Prot (e.g. Brettanomyces
    # has essentially no nuclear-genome entries). ANI compares genome
    # sequence directly and does not depend on reference-DB annotation
    # coverage.
    if [ -n "$EUK_BINS_DIR" ] && [ -n "$(ls -A "$EUK_BINS_DIR" 2>/dev/null)" ]; then
        log "[Fase 5b] Fungal species ID (fastANI)..." | tee -a "$SAMPLE_LOG"
        (
            set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$ENV_FASTANI"; set -u
            fetch_ani_references

            REF_LIST="$ANI_REF_DIR/ref_list.txt"
            find "$ANI_REF_DIR" -name "*.fna" > "$REF_LIST"

            for bin in "$EUK_BINS_DIR"/*.fa; do
                [ -f "$bin" ] || continue
                bin_name=$(basename "$bin" .fa)
                out="$FUNC_DIR/${SRA_ID}/${bin_name}_ani.txt"
                fastANI -q "$bin" --rl "$REF_LIST" -o "$out" 2>>"$SAMPLE_LOG"
                if [ -s "$out" ]; then
                    best=$(sort -t$'\t' -k3,3 -rn "$out" | head -1)
                    ani=$(echo "$best" | cut -f3)
                    if awk -v a="$ani" -v m="$ANI_MIN_IDENTITY" 'BEGIN{exit !(a>=m)}'; then
                        log "  $bin_name -> $(echo "$best" | cut -f2 | xargs basename) (ANI ${ani}%)" | tee -a "$SAMPLE_LOG"
                    else
                        log "  $bin_name -> best hit ${ani}% ANI, below ${ANI_MIN_IDENTITY}% threshold; no confident match" | tee -a "$SAMPLE_LOG"
                    fi
                else
                    log "  $bin_name -> no match against any reference" | tee -a "$SAMPLE_LOG"
                fi
            done
        ) 2>&1 | tee -a "$SAMPLE_LOG"
    fi

    # ---- Fase 6: specific prospection (antiSMASH + dbCAN) ----
    log "[Fase 6] Specific prospection (antiSMASH + dbCAN)..." | tee -a "$SAMPLE_LOG"
    # Not configured in this version. Left as an explicit extension point.

    log "Sample $SRA_ID finished." | tee -a "$SAMPLE_LOG"
done

log "=========================================================="
log " PIPELINE FINISHED"
log " Check $LOG_DIR/*.log per sample for warnings/errors."
log "=========================================================="

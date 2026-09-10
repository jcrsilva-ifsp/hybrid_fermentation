#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -J hybrid_metagenomics_SE
#SBATCH --time=20:00:00
#SBATCH --mem=256G
#SBATCH --output=hibrido_process_%j.log

# ==========================================================================
# Hybrid metagenomics pipeline — single-end variant (e.g. Ion Torrent).
# See pipeline_hibrido_corrigido.sh for the paired-end (Illumina) version
# and for full inline documentation; only single-end-specific differences
# are commented here.
# ==========================================================================

set -uo pipefail

# --------------------------------------------------------------------------
# 0. CONFIGURATION — edit all paths below for your environment
# --------------------------------------------------------------------------
MICROMAMBA_BIN_PEDRO="/path/to/micromamba"

ENV_TIARA="/path/to/env_tiara"
ENV_METAEUK="/path/to/env_metaeuk"
ENV_BUSCO_PATH="/path/to/env_busco"
ENV_KOFAMSCAN="/path/to/env_kofamscan"
ENV_BRACKEN="/path/to/env_bracken"
ENV_FASTANI="/path/to/env_fastani"
PROKARYOTE_ENV="/path/to/metawrap-env"

METAEUK_BIN="$ENV_METAEUK/bin/metaeuk"
METAEUK_REF_DB="/path/to/banco_proteinas/metaeuk_swissprotDB"
BUSCO_BIN="$ENV_BUSCO_PATH/bin/busco"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-fungi_odb12.2}"
BUSCO_DOWNLOAD_PATH="/path/to/busco_downloads/busco_downloads"
KOFAM_PROFILE="/path/to/kofam_db/profiles"
KOFAM_KO_LIST="/path/to/kofam_db/ko_list"
KRAKEN2_DB="${KRAKEN2_DB:-/path/to/banco_kraken}"
BRACKEN_READ_LEN="${BRACKEN_READ_LEN:-200}"   # Ion Torrent reads run longer than Illumina; check actual length
TRIMMOMATIC_JAR="/path/to/trimmomatic-0.39.jar"

ANI_REF_DIR="/path/to/genomas_referencia_leveduras/candidatos_ani"
declare -A ANI_REFS
ANI_REFS["Saccharomyces_cerevisiae"]="GCF_000146045.2"
ANI_REFS["Brettanomyces_bruxellensis"]="GCF_011074885.1"
ANI_REFS["Zygosaccharomyces_rouxii"]="GCA_000026365.1"
ANI_MIN_IDENTITY=90

BASE_DIR="/path/to/project"
THREADS=20
MIN_CONTIG_LEN=1500

SRA_IDS=("${SAMPLE_ID:?Set SAMPLE_ID via --export=ALL,SAMPLE_ID=<accession>}")

RAW_DIR="$BASE_DIR/raw_reads"
QC_DIR="$BASE_DIR/qc_reads"
KRAKEN_DIR="$BASE_DIR/kraken_bracken"
ASSEMBLY_DIR="$BASE_DIR/assembly"
BIFURCATION_DIR="$BASE_DIR/bifurcation"
BIN_PRO_DIR="$BASE_DIR/binning_prokariotos"
BIN_EUK_DIR="$BASE_DIR/binning_eukariotos"
FUNC_DIR="$BASE_DIR/anotacao_funcional"
PROSPEC_DIR="$BASE_DIR/prospeccao_especifica"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_DIR" "$QC_DIR" "$KRAKEN_DIR" "$ASSEMBLY_DIR" "$BIFURCATION_DIR" \
         "$BIN_PRO_DIR" "$BIN_EUK_DIR" "$FUNC_DIR" "$PROSPEC_DIR" "$LOG_DIR" "$ANI_REF_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

check_step() {
    local step_name="$1" output_path="$2"
    if [ -s "$output_path" ]; then
        log "OK: step '$step_name' completed."; return 0
    else
        log "ERROR: step '$step_name' produced no output ($output_path). Skipping sample."; return 1
    fi
}

fetch_ani_references() {
    for name in "${!ANI_REFS[@]}"; do
        local acc="${ANI_REFS[$name]}"
        local fna
        fna=$(find "$ANI_REF_DIR" -path "*${acc}*" -name "*.fna" 2>/dev/null | head -1)
        if [ -z "$fna" ]; then
            log "Downloading ANI reference $name ($acc)..."
            ( cd "$ANI_REF_DIR" && datasets download genome accession "$acc" --include genome --filename "${acc}.zip" \
                && unzip -o "${acc}.zip" -d "$acc" >/dev/null )
        fi
    done
}

# ==========================================================================
# MAIN LOOP
# ==========================================================================
for SRA_ID in "${SRA_IDS[@]}"; do
    SAMPLE_LOG="$LOG_DIR/${SRA_ID}.log"
    log "=== Processing $SRA_ID (single-end) ===" | tee -a "$SAMPLE_LOG"

    # ---- Fase 1: QC (Trimmomatic, SE mode). No bmtagger: single-end
    # decontamination is uncommon and not needed for beverage samples. ----
    READS_R1="$QC_DIR/${SRA_ID}_SE.fq"
    if [ -s "$READS_R1" ]; then
        log "[Fase 1] Existing clean reads found, skipping Trimmomatic." | tee -a "$SAMPLE_LOG"
    else
        log "[Fase 1] Trimmomatic (SE)..." | tee -a "$SAMPLE_LOG"
        java -jar "$TRIMMOMATIC_JAR" SE -threads $THREADS \
            "$RAW_DIR/${SRA_ID}.fastq.gz" "$READS_R1" \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2>&1 | tee -a "$SAMPLE_LOG"
        check_step "Trimmomatic" "$READS_R1" || continue
    fi

    # ---- Fase 1c: Kraken2 + Bracken (no --paired) ----
    log "[Fase 1c] Kraken2+Bracken..." | tee -a "$SAMPLE_LOG"
    mkdir -p "$KRAKEN_DIR/$SRA_ID"
    if [ -d "$KRAKEN2_DB" ]; then
        (
            set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$ENV_BRACKEN"; set -u
            kraken2 --db "$KRAKEN2_DB" --threads $THREADS --memory-mapping \
                --report "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_kraken2_report.txt" \
                --output "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_kraken2_output.txt" "$READS_R1"
            bracken -d "$KRAKEN2_DB" -i "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_kraken2_report.txt" \
                -o "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_bracken_species.txt" \
                -w "$KRAKEN_DIR/$SRA_ID/${SRA_ID}_bracken_report.txt" -r "$BRACKEN_READ_LEN" -l S -t 10
        ) 2>&1 | tee -a "$SAMPLE_LOG"
    else
        log "WARNING: KRAKEN2_DB not found, skipping Fase 1c." | tee -a "$SAMPLE_LOG"
    fi

    # ---- Fase 2: assembly (MEGAHIT, -r for single-end) ----
    log "[Fase 2] MEGAHIT..." | tee -a "$SAMPLE_LOG"
    rm -rf "$ASSEMBLY_DIR/$SRA_ID"
    megahit -r "$READS_R1" -t $THREADS --min-contig-len $MIN_CONTIG_LEN \
        -o "$ASSEMBLY_DIR/$SRA_ID" 2>&1 | tee -a "$SAMPLE_LOG"
    check_step "MEGAHIT" "$ASSEMBLY_DIR/$SRA_ID/final.contigs.fa" || continue

    # ---- Fase 3: Tiara ----
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

    # ---- Fase 4A: bacterial binning (MetaWRAP --single-end) ----
    log "[Fase 4A] Bacterial MAG recovery..." | tee -a "$SAMPLE_LOG"
    set +u; eval "$("$MICROMAMBA_BIN_PEDRO" shell hook -s bash)"; micromamba activate -p "$PROKARYOTE_ENV"; set -u

    METAWRAP_READS_DIR="$BIN_PRO_DIR/$SRA_ID/_reads_link"
    mkdir -p "$METAWRAP_READS_DIR"
    ln -sf "$(readlink -f "$READS_R1")" "$METAWRAP_READS_DIR/${SRA_ID}.fastq"

    metawrap binning -o "$BIN_PRO_DIR/$SRA_ID" -t $THREADS -a "$PROK_FASTA" \
        --single-end --metabat2 --maxbin2 --concoct \
        "$METAWRAP_READS_DIR/${SRA_ID}.fastq" 2>&1 | tee -a "$SAMPLE_LOG"
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
        log "WARNING: bin_refinement failed (known MetaWRAP bug). Falling back to CheckM directly." | tee -a "$SAMPLE_LOG"
        CHECKM_DIR="$BIN_PRO_DIR/$SRA_ID/${SRA_ID}_checkm_direto"
        mkdir -p "$CHECKM_DIR/bins_filtrados_corrigido"
        checkm lineage_wf -t $THREADS -x fa "$BIN_PRO_DIR/$SRA_ID/metabat2_bins" "$CHECKM_DIR" \
            --reduced_tree --tab_table -f "$CHECKM_DIR/checkm_results.txt" 2>&1 | tee -a "$SAMPLE_LOG"

        awk -F'\t' -v comp=90 -v cont=5 \
            'NR>1 && $12+0 >= comp && $13+0 <= cont {print $1".fa"}' \
            "$CHECKM_DIR/checkm_results.txt" > "$CHECKM_DIR/bins_aprovados.txt"

        while read -r binfile; do
            [ "$binfile" = "bin.unbinned.fa" ] && continue
            cp "$BIN_PRO_DIR/$SRA_ID/metabat2_bins/$binfile" "$CHECKM_DIR/bins_filtrados_corrigido/" 2>/dev/null
        done < "$CHECKM_DIR/bins_aprovados.txt"

        BAC_BINS_FINAL_DIR="$CHECKM_DIR/bins_filtrados_corrigido"
        if [ -z "$(ls -A "$BAC_BINS_FINAL_DIR" 2>/dev/null)" ]; then
            log "WARNING: no bacterial bin met the quality threshold for $SRA_ID." | tee -a "$SAMPLE_LOG"
        fi
    fi

    # ---- Fase 4B: fungal binning (bwa mem single-end) ----
    EUK_BINS_DIR=""
    if [ -s "$EUK_FASTA" ]; then
        log "[Fase 4B] Fungal MAG recovery..." | tee -a "$SAMPLE_LOG"
        mkdir -p "$BIN_EUK_DIR/$SRA_ID/bins"
        bwa index "$EUK_FASTA" 2>&1 | tee -a "$SAMPLE_LOG"
        bwa mem -t $THREADS "$EUK_FASTA" "$READS_R1" 2>>"$SAMPLE_LOG" | \
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
    # Not configured in this version.

    log "Sample $SRA_ID finished." | tee -a "$SAMPLE_LOG"
done

log "=========================================================="
log " PIPELINE FINISHED"
log " Check $LOG_DIR/*.log per sample for warnings/errors."
log "=========================================================="

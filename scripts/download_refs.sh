#!/bin/bash

set -e
set -o pipefail
set -x

echo "Starting reference preparation pipeline..."
echo "========================================"

REFS_DIR=$1
mkdir -p $REFS_DIR
rm -f "$REFS_DIR/*.tmp"

# ==========================================
# 1. CaVEMan Blacklist (ENCODE hg38)
# ==========================================
BLACKLIST_FILE="$REFS_DIR/hg38-blacklist.v2.bed"
BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"

if [ ! -f "$BLACKLIST_FILE" ]; then
    echo "Downloading CaVEMan Blacklist..."
    wget -qO- "$BLACKLIST_URL" | zcat > "$BLACKLIST_FILE.tmp" && mv "$BLACKLIST_FILE.tmp" "$BLACKLIST_FILE" 
else
    echo "Skipped: $BLACKLIST_FILE already exists."
fi

# ==========================================
# 2. CaVEMan Indels (Mills & 1000G Gold Standard hg38)
# ==========================================
INDELS_FILE="$REFS_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
INDELS_URL="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
INDELS_TBI="$REFS_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
INDELS_TBI_URL="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

if [ ! -f "$INDELS_FILE" ]; then
    echo "Downloading CaVEMan Indels VCF..."
    wget -qO "$INDELS_FILE.tmp" "$INDELS_URL" && mv "$INDELS_FILE.tmp" "$INDELS_FILE"
else
    echo "Skipped: $INDELS_FILE already exists."
fi

if [ ! -f "$INDELS_TBI" ]; then
    echo "Downloading CaVEMan Indels Index..."
    wget -qO "$INDELS_TBI.tmp" "$INDELS_TBI_URL" && mv "$INDELS_TBI.tmp" "$INDELS_TBI"
else
    echo "Skipped: $INDELS_TBI already exists."
fi

# ==========================================
# 3. ASCAT GC Correction
# ==========================================
ASCAT_GC_DIR="$REFS_DIR/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+"
ASCAT_GC_FILE="$ASCAT_GC_DIR/ascat/SnpGcCorrections.tsv"
ASCAT_GC_URL="ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz"

if [ ! -d "$ASCAT_GC_DIR" ]; then
    echo "Downloading ASCAT GC Correction file..."
    wget -qO "$ASCAT_GC_DIR.tar.gz" "$ASCAT_GC_URL"
    tar -zxvf "$ASCAT_GC_DIR.tar.gz" -C "$REFS_DIR"
    # Sort the chromosomes numerically rather than lexicographically;
    # otherwise, the conversion from bigWig to bedGraph will be incorrect.
    head -n 1 $ASCAT_GC_FILE > sorted_snp_gc.tsv
    tail -n +2 $ASCAT_GC_FILE | sort -k2,2V -k3,3n >> sorted_snp_gc.tsv
    cp -f sorted_snp_gc.tsv $ASCAT_GC_FILE
else
    echo "Skipped: $ASCAT_GC_DIR already exists."
fi

# ==========================================
# 4. UCSC Regions (Centromeres, Telomeres, SuperDups)
# ==========================================
CENTROMERES_FILE="$REFS_DIR/centromeres.bed"
if [ ! -f "$CENTROMERES_FILE" ]; then
    echo "Downloading Centromeres..."
    wget -qO- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz" | zcat | cut -f2,3,4 | awk '{print $0"\tcentromere"}'  > "$CENTROMERES_FILE.tmp" && mv "$CENTROMERES_FILE.tmp" "$CENTROMERES_FILE"
else
    echo "Skipped: $CENTROMERES_FILE already exists."
fi

TELOMERES_FILE="$REFS_DIR/telomeres.bed"
if [ ! -f "$TELOMERES_FILE" ]; then
    echo "Downloading Telomeres..."
    wget -qO- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz" | zcat | awk 'BEGIN{FS="\t"} $8=="telomere" {print $2"\t"$3"\t"$4"\ttelomere"}' > "$TELOMERES_FILE.tmp" && mv "$TELOMERES_FILE.tmp" "$TELOMERES_FILE"
else
    echo "Skipped: $TELOMERES_FILE already exists."
fi

SUPERDUPS_FILE="$REFS_DIR/super_dups.bed"
if [ ! -f "$SUPERDUPS_FILE" ]; then
    echo "Downloading Segmental Duplications (SuperDups)..."
    wget -qO- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz" | zcat | cut -f2,3,4 | awk '{print $0"\tsuperDup"}' > "$SUPERDUPS_FILE.tmp" && mv "$SUPERDUPS_FILE.tmp" "$SUPERDUPS_FILE"
else
    echo "Skipped: $SUPERDUPS_FILE already exists."
fi

# ==========================================
# 5. Merge Ultimate Blacklist
# ==========================================
ULTIMATE_BLACKLIST="$REFS_DIR/caveman_blacklist.bed"
if [ ! -f "$ULTIMATE_BLACKLIST" ]; then
    echo "Merging all problematic regions into $ULTIMATE_BLACKLIST..."
    cat "$BLACKLIST_FILE" "$CENTROMERES_FILE" "$TELOMERES_FILE" "$SUPERDUPS_FILE" | sort -k1,1V -k2,2n | bedtools merge > "$ULTIMATE_BLACKLIST.tmp" && mv "$ULTIMATE_BLACKLIST.tmp" "$ULTIMATE_BLACKLIST" 
    echo "Ultimate Blacklist created."
else
    echo "Skipped: $ULTIMATE_BLACKLIST already exists."
fi

echo "========================================"
echo "Pipeline initialization complete!"

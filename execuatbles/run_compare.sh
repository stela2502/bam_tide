#!/usr/bin/env bash
set -euo pipefail

# ---- user tools (override via env if needed) ----
: "${BAMCOVERAGE_PY:=bamCoverage}"    # deeptools bamCoverage
: "${BAMCOVERAGE_RS:=bam-coverage}"   # your rust tool
: "${BWCOMPARE:=bw-compare}"          # your comparator

# ---- args ----
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.bam> [bin_width]"
  exit 2
fi

BAM="$1"
BIN_WIDTH="${2:-50}"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi

# ---- derive names ----
bam_base="$(basename "$BAM")"
bam_stem="${bam_base%.bam}"  # ok even if not .bam; just keeps full name
PY_BW="/tmp/python_${bam_stem}.bw"
RS_BW="/tmp/rust_${bam_stem}.bw"

# Put a copy of the compare output next to bw-compare (your request: "rust infil's folder")
#BWCOMPARE_PATH="$(command -v "$BWCOMPARE")"
#BWCOMPARE_DIR="$(cd "$(dirname "$BWCOMPARE_PATH")" && pwd)"
#REPORT_DIR="$BWCOMPARE_DIR"
#REPORT_DIR="/tmp"
#REPORT_TXT="${REPORT_DIR}/bw_compare_${bam_stem}_w${BIN_WIDTH}.txt"

# Also keep a tmp copy for quick inspection
TMP_REPORT="/tmp/bw_compare_${bam_stem}_w${BIN_WIDTH}.txt"

echo "BAM:        $BAM"
echo "BIN_WIDTH:  $BIN_WIDTH"
echo "PY_BW:      $PY_BW"
echo "RS_BW:      $RS_BW"
echo "REPORT:     $TMP_REPORT"
echo

# ---- 1) create expected python bw (if missing) ----
if [[ -s "$PY_BW" ]]; then
  echo "[1/3] Python BW exists, skipping: $PY_BW"
else
  echo "[1/3] Creating python BW with deeptools bamCoverage..."
  "$BAMCOVERAGE_PY" -b "$BAM" -o "$PY_BW" -bs "$BIN_WIDTH" --minMappingQuality 0
fi

# ---- 2) create rust bw (always regenerate, safer while iterating) ----
echo "[2/3] Creating rust BW with $BAMCOVERAGE_RS..."
"$BAMCOVERAGE_RS" -b "$BAM" -o "$RS_BW" -w "$BIN_WIDTH" --min-mapping-quality 0

# ---- 3) compare ----
echo "[3/3] Comparing..."
# capture stdout+stderr to report files AND show on terminal
set +e
OUT="$("$BWCOMPARE" --python-bw "$PY_BW" --rust-bw "$RS_BW" --bin-width "$BIN_WIDTH" --eps 1e-4 2>&1)"
RC=$?
set -e

echo "$OUT" | tee "$TMP_REPORT" > "$REPORT_TXT"

if [[ $RC -ne 0 ]]; then
  echo "bw-compare exited with code $RC" >&2
  exit $RC
fi

echo
echo "Saved:"
echo "  $TMP_REPORT"
echo "  $REPORT_TXT"

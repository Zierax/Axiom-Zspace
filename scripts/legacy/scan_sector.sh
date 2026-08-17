#!/bin/bash
# ══════════════════════════════════════════════════════════════════
# Axiom-ZSpace C Engine — Sector Scanner (V2 — Local CSV Mode)
# ══════════════════════════════════════════════════════════════════
#
# TWO-PHASE ARCHITECTURE for maximum speed:
#   Phase 1 (Python):  python export_sector_csv.py 5
#                      → Bulk-downloads all FITS, exports to CSV
#   Phase 2 (C/This):  ./scan_sector.sh 5
#                      → Reads local CSVs, pure computation, blazing fast
#
# Usage:
#   ./scan_sector.sh 5           # Scan sector 5 (all targets)
#   ./scan_sector.sh 5 100       # Limit to 100 targets
# ══════════════════════════════════════════════════════════════════

SECTOR=${1:-5}
MAX_TARGETS=${2:-0}
BINARY="./axiom_zspace"
CSV_DIR=".cache/sector_${SECTOR}_csv"
MANIFEST="${CSV_DIR}/manifest.json"
OUTPUT_DIR="axiom_output/sector_${SECTOR}_c"
DISCOVERIES_FILE="${OUTPUT_DIR}/discoveries.json"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo ""
echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD}  AXIOM-ZSPACE C ENGINE — SECTOR ${SECTOR} SCAN${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}"
echo ""

# Check binary
if [ ! -f "$BINARY" ]; then
    echo -e "${RED}ERROR: $BINARY not found. Compile first:${NC}"
    echo "  gcc -O3 -fopenmp -o axiom_zspace axiom_zspace.c -lm"
    exit 1
fi

# Check CSV directory
if [ ! -d "$CSV_DIR" ]; then
    echo -e "${YELLOW}  CSV data not found. Running bulk export first...${NC}"
    echo ""
    if [ "$MAX_TARGETS" -gt 0 ]; then
        python3 export_sector_csv.py "$SECTOR" --max "$MAX_TARGETS"
    else
        python3 export_sector_csv.py "$SECTOR"
    fi
    echo ""
fi

# Get list of CSV files
CSV_FILES=$(find "$CSV_DIR" -name "TIC_*.csv" -size +100c 2>/dev/null | sort)
if [ -z "$CSV_FILES" ]; then
    echo -e "${RED}ERROR: No CSV files found in ${CSV_DIR}${NC}"
    echo "  Run first:  python3 export_sector_csv.py ${SECTOR}"
    exit 1
fi

TOTAL=$(echo "$CSV_FILES" | wc -l)
if [ "$MAX_TARGETS" -gt 0 ] && [ "$MAX_TARGETS" -lt "$TOTAL" ]; then
    CSV_FILES=$(echo "$CSV_FILES" | head -n "$MAX_TARGETS")
    TOTAL=$MAX_TARGETS
    echo -e "${YELLOW}  Limiting to ${MAX_TARGETS} targets${NC}"
fi

echo -e "${CYAN}  Mode:    LOCAL CSV (no network — pure computation)${NC}"
echo -e "${CYAN}  Targets: ${TOTAL}${NC}"
echo -e "${CYAN}  Source:  ${CSV_DIR}/${NC}"
echo ""

# Create output directory
mkdir -p "${OUTPUT_DIR}/discoveries" 2>/dev/null

# ── Processing Loop ──────────────────────────────────────────────
PROCESSED=0
DISCOVERIES=0
PLANET_CANDIDATES=0
LIKELY_CANDIDATES=0
FALSE_POSITIVES=0
FAILED=0
START_TIME=$(date +%s)

# Start discoveries JSON
echo '{' > "$DISCOVERIES_FILE"
echo "  \"sector\": ${SECTOR}," >> "$DISCOVERIES_FILE"
echo '  "engine": "axiom_zspace_c_v2.0",' >> "$DISCOVERIES_FILE"
echo '  "planets": [' >> "$DISCOVERIES_FILE"
FIRST_DISC=1

while IFS= read -r CSV_FILE; do
    [ -z "$CSV_FILE" ] && continue

    # Extract TIC ID from filename
    TIC_ID=$(basename "$CSV_FILE" .csv | sed 's/TIC_//')
    
    PROCESSED=$((PROCESSED + 1))
    
    # Timing
    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TIME))
    [ "$ELAPSED" -lt 1 ] && ELAPSED=1
    RATE=$(echo "scale=1; $PROCESSED * 60 / $ELAPSED" | bc 2>/dev/null || echo "?")
    ETA_MIN=$(echo "scale=1; ($TOTAL - $PROCESSED) * $ELAPSED / $PROCESSED / 60" | bc 2>/dev/null || echo "?")

    # Progress bar
    PCT=$((PROCESSED * 100 / TOTAL))
    FILLED=$((PROCESSED * 30 / TOTAL))
    EMPTY=$((30 - FILLED))
    BAR=""
    for ((i=0; i<FILLED; i++)); do BAR="${BAR}█"; done
    for ((i=0; i<EMPTY; i++)); do BAR="${BAR}░"; done

    printf "\r  [%s] %3d%%  (%d/%d)  TIC %-10s | ★:%d FP:%d F:%d | %s/min ETA:%sm   " \
        "$BAR" "$PCT" "$PROCESSED" "$TOTAL" "$TIC_ID" "$DISCOVERIES" "$FALSE_POSITIVES" "$FAILED" "$RATE" "$ETA_MIN"

    # ── Load stellar metadata if available ────────────────────────
    META_FILE="${CSV_DIR}/TIC_${TIC_ID}.meta"
    TIC_FLAG=""
    if [ -f "$META_FILE" ]; then
        TIC_FLAG="--tic ${TIC_ID}"
    fi

    # ── Run C engine on local CSV ──────────────────────────────────
    OUTPUT=$($BINARY --csv "$CSV_FILE" --tic "$TIC_ID" 2>&1)
    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        FAILED=$((FAILED + 1))
        continue
    fi

    # Parse results from output
    CVS=$(echo "$OUTPUT" | grep -oP 'CVS\s*=\s*\K[0-9.]+' | head -1 | tr -d '\r')
    PERIOD=$(echo "$OUTPUT" | grep -oP 'Period:\s+\K[0-9.]+' | head -1 | tr -d '\r')
    VERDICT=$(echo "$OUTPUT" | grep -oP 'CVS\s*=.*\|\s+\K[^|]+' | head -1 | sed 's/[[:space:]]*$//' | tr -d '\r')
    SNR=$(echo "$OUTPUT" | grep -oP 'SNR:\s+\K[0-9.]+' | head -1 | tr -d '\r')

    if [ -z "$CVS" ]; then
        FAILED=$((FAILED + 1))
        continue
    fi

    # Classify result
    IS_PLANET=$(awk -v cvs="$CVS" 'BEGIN {print (cvs >= 0.80) ? 1 : 0}')
    IS_LIKELY=$(awk -v cvs="$CVS" 'BEGIN {print (cvs >= 0.55) ? 1 : 0}')

    if [ "$IS_PLANET" = "1" ]; then
        DISCOVERIES=$((DISCOVERIES + 1))
        PLANET_CANDIDATES=$((PLANET_CANDIDATES + 1))
        ICON="★"
        LABEL="PLANET CANDIDATE"
    elif [ "$IS_LIKELY" = "1" ]; then
        DISCOVERIES=$((DISCOVERIES + 1))
        LIKELY_CANDIDATES=$((LIKELY_CANDIDATES + 1))
        ICON="◆"
        LABEL="LIKELY CANDIDATE"
    else
        FALSE_POSITIVES=$((FALSE_POSITIVES + 1))
        continue
    fi

    # Move discovery card
    CARD="discovery_card_ZS-T-${TIC_ID}-01.json"
    [ -f "$CARD" ] && mv "$CARD" "${OUTPUT_DIR}/discoveries/"

    # Append to discoveries JSON
    [ "$FIRST_DISC" -eq 0 ] && echo "," >> "$DISCOVERIES_FILE"
    FIRST_DISC=0
    cat >> "$DISCOVERIES_FILE" << EOD
    {
      "#": ${DISCOVERIES},
      "tic_id": "${TIC_ID}",
      "zspace_id": "ZS-T-${TIC_ID}-01",
      "period_days": ${PERIOD:-0},
      "snr": ${SNR:-0},
      "cvs": ${CVS},
      "verdict": "${LABEL}"
    }
EOD

    # Print discovery
    echo ""
    echo -e "  ${GREEN}${ICON} #${DISCOVERIES} TIC ${TIC_ID} | P=${PERIOD}d | SNR=${SNR} | CVS=${CVS} | ${LABEL}${NC}"

done <<< "$CSV_FILES"

# ── Finalize ────────────────────────────────────────────────────
echo ""
echo ""

NOW=$(date +%s)
TOTAL_SEC=$((NOW - START_TIME))
TOTAL_MIN=$(echo "scale=1; $TOTAL_SEC / 60" | bc 2>/dev/null || echo "$((TOTAL_SEC/60))")
FINAL_RATE=$(echo "scale=1; $PROCESSED * 60 / $TOTAL_SEC" | bc 2>/dev/null || echo "?")

# Close discoveries JSON
echo "" >> "$DISCOVERIES_FILE"
echo "  ]," >> "$DISCOVERIES_FILE"
echo "  \"total_discoveries\": ${DISCOVERIES}," >> "$DISCOVERIES_FILE"
echo "  \"planet_candidates\": ${PLANET_CANDIDATES}," >> "$DISCOVERIES_FILE"
echo "  \"likely_candidates\": ${LIKELY_CANDIDATES}," >> "$DISCOVERIES_FILE"
echo "  \"false_positives\": ${FALSE_POSITIVES}," >> "$DISCOVERIES_FILE"
echo "  \"failed\": ${FAILED}," >> "$DISCOVERIES_FILE"
echo "  \"total_scanned\": ${PROCESSED}," >> "$DISCOVERIES_FILE"
echo "  \"elapsed_minutes\": ${TOTAL_MIN}," >> "$DISCOVERIES_FILE"
echo "  \"rate_per_minute\": ${FINAL_RATE}," >> "$DISCOVERIES_FILE"
echo "  \"scan_date\": \"$(date -u +%Y-%m-%dT%H:%M:%SZ)\"" >> "$DISCOVERIES_FILE"
echo "}" >> "$DISCOVERIES_FILE"

# Clean up stray cards
rm -f discovery_card_ZS-T-*-01.json 2>/dev/null

echo -e "${BOLD}  ╔══════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}  ║  SECTOR ${SECTOR} C-SCAN COMPLETE                                ║${NC}"
printf  "  ║  Processed:  %5d / %-5d                           ║\n" "$PROCESSED" "$TOTAL"
printf  "  ║  ★ Planets:  %5d  |  ◆ Likely: %5d                ║\n" "$PLANET_CANDIDATES" "$LIKELY_CANDIDATES"
printf  "  ║  ✗ FP:       %5d  |  Failed:  %5d                ║\n" "$FALSE_POSITIVES" "$FAILED"
printf  "  ║  Time: %s min  |  Rate: %s/min                  ║\n" "$TOTAL_MIN" "$FINAL_RATE"
echo -e "${BOLD}  ╚══════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "  Results: ${CYAN}${DISCOVERIES_FILE}${NC}"
echo ""

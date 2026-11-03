#!/usr/bin/env bash
# Convenience wrapper: ./run_pipeline.sh <input.adat> [output_dir] [group_column]
# Runs the pipeline in Docker if available, otherwise in the active conda env.
set -euo pipefail
IN="${1:?usage: ./run_pipeline.sh <input.adat> [output_dir] [group_column]}"
OUT="${2:-somascan_qc_out}"
GRP="${3:-}"
IN_DIR="$(cd "$(dirname "$IN")" && pwd)"; IN_BASE="$(basename "$IN")"
mkdir -p "$OUT"; OUT_DIR="$(cd "$OUT" && pwd)"
ARGS=(run -i "/in/$IN_BASE" -o /out)
[ -n "$GRP" ] && ARGS+=(-g "$GRP")
if command -v docker >/dev/null 2>&1 && docker image inspect somascan-qc:1.1.0 >/dev/null 2>&1; then
  echo "[somascan-qc] running in Docker"
  docker run --rm -v "$IN_DIR":/in:ro -v "$OUT_DIR":/out somascan-qc:1.1.0 "${ARGS[@]}"
else
  echo "[somascan-qc] running locally (conda/pip install)"
  somascan-qc run -i "$IN" -o "$OUT_DIR" ${GRP:+-g "$GRP"}
fi
echo "[somascan-qc] open $OUT_DIR/qc_report.html"

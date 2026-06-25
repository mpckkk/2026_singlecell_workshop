#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PY="/Library/Developer/CommandLineTools/usr/bin/python3"
export OMP_NUM_THREADS=4

mkdir -p logs figures tables objects

echo "=== Phase 0: QC ==="
Rscript scripts/00_qc.R 2>&1 | tee logs/00_qc.log

echo "=== Phase 1: Seurat ==="
Rscript scripts/01_seurat_sc.R 2>&1 | tee logs/01_seurat.log

echo "=== Phase 3: DESeq2 bulk ==="
Rscript scripts/03_deseq2_bulk.R 2>&1 | tee logs/03_deseq2.log

echo "=== Phase 2: PAGA (Scanpy) ==="
"$PY" scripts/02_paga.py 2>&1 | tee logs/02_paga.log

echo "=== Phase 4: Co-expression ==="
Rscript scripts/04_coexpression.R 2>&1 | tee logs/04_coexpression.log

echo "=== Phase 5: BisqueRNA ==="
Rscript scripts/05_bisque.R 2>&1 | tee logs/05_bisque.log

echo "=== Phase 5: scVI ==="
"$PY" scripts/05_scvi.py 2>&1 | tee logs/05_scvi.log

echo "=== Phase 6: Report ==="
Rscript scripts/06_report.R 2>&1 | tee logs/06_report.log

"$PY" -m pip freeze > logs/pip_freeze.txt
Rscript -e 'writeLines(capture.output(sessionInfo()), "logs/sessionInfo_R.txt")'

echo "Pipeline complete."

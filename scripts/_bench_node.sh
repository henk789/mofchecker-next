#!/bin/bash
# Run throughput benchmarks on a dedicated Capella compute node.
# /projects is read-only on compute nodes; outputs go to node-local /tmp.
set -e
cd /projects/p2/p_fm_mofs/mofchecker-next
echo "=== node: $(hostname)  cores(cgroup): $(nproc) ==="
echo
echo "########## EVAL set (accelerated diagnostics, no has_oms) ##########"
.venv/bin/python scripts/benchmark_throughput.py \
    --structures-file bench_structures_150.json --workers 10 --set eval \
    --out /tmp/throughput_eval.json
echo
echo "########## FULL set (all geometric diagnostics, incl has_oms) ##########"
.venv/bin/python scripts/benchmark_throughput.py \
    --structures-file bench_structures_150.json --workers 10 --set full \
    --out /tmp/throughput_full.json
echo "=== BENCH_NODE_DONE ==="

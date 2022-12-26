#!/usr/bin/env python3
"""
ALBEDO - Hyperspectral Processing Engine
scripts/run_benchmarks.py : Benchmark Execution & Markdown Table Generator

Runs the high-precision J microbenchmark suite and formats algorithmic
timings, throughput, and complexity metrics for documentation.
"""

import os
import re
import subprocess
import sys

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bench_file = os.path.join(repo_root, "scripts", "benchmarks.ijs")
    
    jconsole_path = "/Users/aneeshkrishna/.local/bin/jconsole"
    if not os.path.exists(jconsole_path):
        jconsole_path = "jconsole"
        
    print("=" * 80)
    print("               ALBEDO MICROBENCHMARK & PERFORMANCE AUDIT")
    print("               Interpreter: " + jconsole_path)
    print("=" * 80)
    
    cmd = [jconsole_path, bench_file]
    proc = subprocess.run(
        cmd,
        cwd=repo_root,
        input="exit 0\n",
        capture_output=True,
        text=True,
        timeout=120
    )
    
    if proc.returncode != 0:
        print("ERROR: Benchmark failed with exit code:", proc.returncode)
        print("STDOUT:\n", proc.stdout)
        print("STDERR:\n", proc.stderr)
        sys.exit(1)
        
    lines = proc.stdout.splitlines()
    bench_results = []
    
    pattern = re.compile(r"^\s*(\d+)\.\s*(.+?):\s*(\d+)\s*ms")
    
    complexity_map = {
        1: ("Thin SVD (50x20)", "O(m n^2)", "Dense Jacobi SVD"),
        2: ("Cholesky Factorization (100x100)", "O(n^3 / 3)", "Cholesky-Banachiewicz"),
        3: ("Covariance Estimation (50x500)", "O(p^2 N)", "Two-pass shrinkage"),
        4: ("VCA Simplex Extraction (p=4, N=1000)", "O(p^2 N)", "Vertex Component Analysis"),
        5: ("Reed-Xiaoli Detection (20x25x20)", "O(p^2 N)", "Quadratic Mahalanobis metric"),
        6: ("Tensor Mode-0 Contraction (20x20x20)", "O(I J K J')", "Matricized tensor product"),
        7: ("4-Level Wavelet DWT (4096 samples)", "O(N)", "Mallat pyramidal filterbank"),
        8: ("Fractional Derivative (alpha=0.5, 224 bands)", "O(B^2)", "Grunwald-Letnikov differencing"),
        9: ("Tucker3 HOSVD (15x15x15, rank 4)", "O(I1 I2 I3 R)", "Higher-Order SVD"),
        10: ("CP-ALS Factorization (10x10x10, rank 3)", "O(R^3 + R I J K)", "Alternating Least Squares")
    }
    
    print("\nBenchmark Raw Output:")
    print("-" * 80)
    for line in lines:
        if line.strip():
            print(line)
        match = pattern.match(line)
        if match:
            idx = int(match.group(1))
            name = match.group(2).strip()
            ms = int(match.group(3))
            bench_results.append((idx, name, ms))
            
    print("\n" + "=" * 80)
    print("                    PERFORMANCE BENCHMARK TABLE")
    print("=" * 80)
    print(f"| {'#':<3} | {'Benchmark Kernel':<38} | {'Complexity':<18} | {'Latency (ms)':<12} | {'Implementation':<28} |")
    print(f"|{'-'*5}|{'-'*40}|{'-'*20}|{'-'*14}|{'-'*30}|")
    
    for idx, name, ms in bench_results:
        meta = complexity_map.get(idx, (name, "O(N)", "Native J"))
        comp = meta[1]
        impl = meta[2]
        print(f"| {idx:<3} | {meta[0]:<38} | {comp:<18} | {ms:<12} | {impl:<28} |")
        
    print("-" * 80)
    print(f"Total benchmarks executed: {len(bench_results)}")
    print("Status: ALL BENCHMARKS COMPLETED SUCCESSFULLY\n")

if __name__ == "__main__":
    main()

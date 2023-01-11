# ALBEDO: Hyperspectral Imaging and Satellite Spectroscopy Processing Engine in J
# Root Build, Benchmark, and Test Automation Makefile

PYTHON ?= python3
PROJECT_DIR = $(shell pwd)

.PHONY: all test check bench loc sim unmix detect viewer clean help

all: check test bench

test:
	@echo "==> Running full ALBEDO test suite (35 suites)..."
	@$(PYTHON) scripts/run_tests.py

check:
	@echo "==> Verifying zero em-dash compliance..."
	@$(PYTHON) scripts/verify_no_em_dash.py

bench:
	@echo "==> Running high-precision microbenchmark suite..."
	@$(PYTHON) scripts/run_benchmarks.py

loc:
	@echo "==> Computing codebase line statistics..."
	@$(PYTHON) scripts/check_loc.py

sim:
	@echo "==> Simulating synthetic 32x32 hyperspectral scene..."
	@./bin/albedo sim --rows 32 --cols 32

unmix:
	@echo "==> Running VCA endmember extraction and FCLS abundance estimation..."
	@./bin/albedo unmix --endmembers 4 --method vca

detect:
	@echo "==> Running Global Reed-Xiaoli (RX) anomaly detection..."
	@./bin/albedo detect --method rx

viewer:
	@echo "==> Starting Albedo WebGL Dashboard on http://localhost:8000..."
	@$(PYTHON) viewer/server.py --port 8000

clean:
	@rm -rf build/ __pycache__ */__pycache__ *.tmp
	@echo "==> Workspace cleaned."

help:
	@echo "ALBEDO Make Targets:"
	@echo "  make test      Run all 35 verification test suites"
	@echo "  make bench     Run microbenchmark performance audit"
	@echo "  make check     Enforce typography and styling invariants"
	@echo "  make loc       Compute lines of code breakdown"
	@echo "  make sim       Generate synthetic hyperspectral cube"
	@echo "  make unmix     Execute VCA/FCLS spectral unmixing demo"
	@echo "  make detect    Execute Global Reed-Xiaoli anomaly detector"
	@echo "  make viewer    Launch WebGL hyperspectral visualizer"

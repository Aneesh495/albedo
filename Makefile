# Albedo: Hyperspectral Imaging and Satellite Spectroscopy Processing Engine in J
# Root Build and Test Automation Makefile

PYTHON ?= python3
PROJECT_DIR = $(shell pwd)

.PHONY: all test check loc clean

all: check test

test:
	@echo "==> Running Albedo test suites..."
	@$(PYTHON) scripts/run_tests.py

check:
	@echo "==> Checking em-dash constraints..."
	@$(PYTHON) scripts/verify_no_em_dash.py

loc:
	@echo "==> Reporting lines of code..."
	@$(PYTHON) scripts/check_loc.py

clean:
	@rm -rf build/ __pycache__ scripts/__pycache__ tests/__pycache__ *.tmp

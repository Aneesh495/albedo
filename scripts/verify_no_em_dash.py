#!/usr/bin/env python3
"""
Verification script: Enforce zero em-dashes in comments, documentation, and source code.
Exits with code 1 if any em-dash (\\u2014) is discovered.
"""

import os
import sys

EM_DASH = "\u2014"

def scan_file(filepath):
    errors = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            for lineno, line in enumerate(f, start=1):
                if EM_DASH in line:
                    errors.append((lineno, line.strip()))
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
    return errors

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    violations = 0
    scanned_count = 0

    skip_dirs = {".git", "build", "vendor"}

    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d not in skip_dirs]
        for fname in filenames:
            ext = os.path.splitext(fname)[1].lower()
            if ext in [".ijs", ".md", ".txt", ".py", ".c", ".h", ".js", ".html", ".css", ".json", ".yaml", ".yml", ".sh"]:
                fpath = os.path.join(dirpath, fname)
                scanned_count += 1
                errs = scan_file(fpath)
                if errs:
                    rel_path = os.path.relpath(fpath, root_dir)
                    for lineno, line in errs:
                        print(f"EM DASH DETECTED: {rel_path}:{lineno}: {line}", file=sys.stderr)
                        violations += 1

    print(f"Scanned {scanned_count} files for em-dash violations.")
    if violations > 0:
        print(f"FAILED: Found {violations} em-dash violation(s).", file=sys.stderr)
        sys.exit(1)
    else:
        print("PASSED: Zero em-dashes found across entire codebase.")
        sys.exit(0)

if __name__ == "__main__":
    main()

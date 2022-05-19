#!/usr/bin/env python3
"""
Lines of code counter for Albedo.
Analyzes J (.ijs) code, supporting Python, C, Markdown, and test suites.
"""

import os
import sys

def count_lines(filepath):
    total = 0
    code = 0
    comments = 0
    blanks = 0

    ext = os.path.splitext(filepath)[1].lower()

    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                total += 1
                stripped = line.strip()
                if not stripped:
                    blanks += 1
                elif ext == ".ijs":
                    if stripped.startswith("NB."):
                        comments += 1
                    else:
                        code += 1
                elif ext in [".py", ".sh"]:
                    if stripped.startswith("#"):
                        comments += 1
                    else:
                        code += 1
                elif ext in [".c", ".h", ".js"]:
                    if stripped.startswith("//") or stripped.startswith("/*") or stripped.startswith("*"):
                        comments += 1
                    else:
                        code += 1
                elif ext == ".md":
                    comments += 1
                else:
                    code += 1
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)

    return total, code, comments, blanks

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    stats = {}
    lang_stats = {}

    skip_dirs = {".git", "build", "vendor"}

    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d not in skip_dirs]
        for fname in filenames:
            ext = os.path.splitext(fname)[1].lower()
            if ext in [".ijs", ".py", ".c", ".h", ".js", ".html", ".css", ".md", ".sh"]:
                fpath = os.path.join(dirpath, fname)
                tot, c, comm, bl = count_lines(fpath)
                rel_dir = os.path.relpath(dirpath, root_dir)
                top_dir = rel_dir.split(os.sep)[0] if rel_dir != "." else "root"

                if top_dir not in stats:
                    stats[top_dir] = [0, 0, 0, 0]
                stats[top_dir][0] += tot
                stats[top_dir][1] += c
                stats[top_dir][2] += comm
                stats[top_dir][3] += bl

                lang = "J (.ijs)" if ext == ".ijs" else ext
                if lang not in lang_stats:
                    lang_stats[lang] = [0, 0, 0, 0]
                lang_stats[lang][0] += tot
                lang_stats[lang][1] += c
                lang_stats[lang][2] += comm
                lang_stats[lang][3] += bl

    print("================================================================================")
    print("                      ALBEDO LINES OF CODE REPORT                               ")
    print("================================================================================")
    print(f"{'Language':<15} {'Total Lines':>12} {'Code Lines':>12} {'Comment Lines':>15} {'Blank Lines':>12}")
    print("--------------------------------------------------------------------------------")
    total_all = [0, 0, 0, 0]
    for lang, counts in sorted(lang_stats.items(), key=lambda x: -x[1][0]):
        print(f"{lang:<15} {counts[0]:>12} {counts[1]:>12} {counts[2]:>15} {counts[3]:>12}")
        for i in range(4):
            total_all[i] += counts[i]
    print("--------------------------------------------------------------------------------")
    print(f"{'TOTAL':<15} {total_all[0]:>12} {total_all[1]:>12} {total_all[2]:>15} {total_all[3]:>12}")
    print("\nBreakdown by Directory:")
    print(f"{'Directory':<20} {'Total':>10} {'Code':>10} {'Comments':>12} {'Blanks':>10}")
    print("--------------------------------------------------------------------------------")
    for d, counts in sorted(stats.items(), key=lambda x: -x[1][0]):
        print(f"{d:<20} {counts[0]:>10} {counts[1]:>10} {counts[2]:>12} {counts[3]:>10}")
    print("================================================================================")

if __name__ == "__main__":
    main()

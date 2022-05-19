#!/usr/bin/env python3
"""
Test runner for Albedo J test suites.
Discovers and executes all test_*.ijs files via jconsole.
"""

import os
import sys
import glob
import subprocess
import time

def find_jconsole():
    # Priority: repo bin/, user local bin, system PATH, /Users/aneeshkrishna/.local/bin/jconsole
    candidates = [
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin", "jconsole"),
        "/Users/aneeshkrishna/.local/bin/jconsole",
        "/Users/aneeshkrishna/.local/share/j9.6/bin/jconsole",
        "/tmp/j96/j9.6/bin/jconsole"
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    # Check PATH
    import shutil
    p = shutil.which("jconsole")
    if p:
        return p
    return None

def run_suite(jconsole_bin, test_file, project_root):
    # Execute test by feeding:
    # 1. cd to project_root
    # 2. load test_framework
    # 3. load test_file
    # 4. exit 0
    cmd_input = f"""
1!:44 < '{project_root}'
load '{test_file}'
exit 0
"""
    start_time = time.time()
    try:
        proc = subprocess.Popen(
            [jconsole_bin],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=project_root
        )
        stdout, stderr = proc.communicate(input=cmd_input, timeout=30)
        elapsed = time.time() - start_time
        return proc.returncode, stdout, stderr, elapsed
    except subprocess.TimeoutExpired:
        proc.kill()
        return -1, "", "TIMEOUT after 30 seconds", 30.0
    except Exception as e:
        return -1, "", str(e), 0.0

def main():
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tests_dir = os.path.join(project_root, "tests")

    jconsole = find_jconsole()
    if not jconsole:
        print("ERROR: jconsole executable not found. Please install J9.6 or configure PATH.", file=sys.stderr)
        sys.exit(1)

    # Collect all test_*.ijs files (excluding test_framework.ijs)
    test_files = sorted(glob.glob(os.path.join(tests_dir, "test_*.ijs")))
    test_files = [f for f in test_files if not f.endswith("test_framework.ijs")]

    if not test_files:
        print("No test files found in tests/")
        sys.exit(0)

    print("================================================================================")
    print("                     ALBEDO AUTOMATED TEST SUITE RUNNER                         ")
    print(f"                     J Interpreter: {jconsole}")
    print(f"                     Discovered {len(test_files)} Test Suites")
    print("================================================================================")

    passed_suites = 0
    failed_suites = 0
    total_time = 0.0

    for test_path in test_files:
        rel_name = os.path.relpath(test_path, project_root)
        print(f"--> Running {rel_name:<45} ... ", end="", flush=True)

        code, stdout, stderr, elapsed = run_suite(jconsole, test_path, project_root)
        total_time += elapsed

        # Check for errors in output or non-zero exit code
        has_error = False
        if code != 0:
            has_error = True
        if "|error" in stdout or "|value error" in stdout or "|syntax error" in stdout or "[FAIL]" in stdout:
            has_error = True
        if "FAILED" in stdout:
            has_error = True

        if not has_error:
            print(f"[\033[92mPASS\033[0m] ({elapsed:.3f}s)")
            passed_suites += 1
        else:
            print(f"[\033[91mFAIL\033[0m] ({elapsed:.3f}s)")
            failed_suites += 1
            print("--------------------------------------------------------------------------------")
            if stdout:
                print("STDOUT:")
                for line in stdout.strip().split("\n"):
                    print(f"  {line}")
            if stderr:
                print("STDERR:")
                for line in stderr.strip().split("\n"):
                    print(f"  {line}")
            print("--------------------------------------------------------------------------------")

    print("================================================================================")
    print(f"RESULTS: {passed_suites} PASSED, {failed_suites} FAILED out of {len(test_files)} suites in {total_time:.2f}s")
    print("================================================================================")

    if failed_suites > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()

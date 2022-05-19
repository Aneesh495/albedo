NB. ==============================================================================
NB. Albedo Test Framework for J
NB. Provides lightweight assertion primitives, numerical tolerance checks,
NB. and structured test reporting for test suites.
NB. ==============================================================================

cocurrent 'z'

NB. Global counters
TEST_PASS_COUNT=: 0
TEST_FAIL_COUNT=: 0
CURRENT_SUITE_NAME=: 'Unnamed Suite'

NB. Initialize suite
suite_begin=: 3 : 0
  CURRENT_SUITE_NAME=: y
  TEST_PASS_COUNT=: 0
  TEST_FAIL_COUNT=: 0
  smoutput '=== Running Suite: ', CURRENT_SUITE_NAME, ' ==='
)

NB. Assert condition is truthy (1)
assert_true=: 4 : 0
  msg=. x
  cond=. y
  if. cond do.
    TEST_PASS_COUNT=: TEST_PASS_COUNT + 1
  else.
    TEST_FAIL_COUNT=: TEST_FAIL_COUNT + 1
    smoutput '  [FAIL] Assertion failed: ', msg
  end.
)

NB. Assert condition is falsey (0)
assert_false=: 4 : 0
  msg=. x
  cond=. y
  msg assert_true 0 = cond
)

NB. Assert exact equality (match)
assert_equal=: 4 : 0
  msg=. x
  'actual expected'=. y
  if. actual -: expected do.
    TEST_PASS_COUNT=: TEST_PASS_COUNT + 1
  else.
    TEST_FAIL_COUNT=: TEST_FAIL_COUNT + 1
    smoutput '  [FAIL] Equality failed: ', msg
    smoutput '    Expected: ', (": expected)
    smoutput '    Actual:   ', (": actual)
  end.
)

NB. Assert floating point near within tolerance
assert_near=: 4 : 0
  msg=. x
  'actual expected tol'=. y
  diff=. | actual - expected
  max_diff=. >./ , diff
  if. max_diff <: tol do.
    TEST_PASS_COUNT=: TEST_PASS_COUNT + 1
  else.
    TEST_FAIL_COUNT=: TEST_FAIL_COUNT + 1
    smoutput '  [FAIL] Float near check failed: ', msg
    smoutput '    Max diff: ', (": max_diff), ' (tol: ', (": tol), ')'
  end.
)

NB. Assert array shape matches expected
assert_shape=: 4 : 0
  msg=. x
  'actual expected_shape'=. y
  act_shape=. $ actual
  if. act_shape -: expected_shape do.
    TEST_PASS_COUNT=: TEST_PASS_COUNT + 1
  else.
    TEST_FAIL_COUNT=: TEST_FAIL_COUNT + 1
    smoutput '  [FAIL] Shape check failed: ', msg
    smoutput '    Expected shape: ', (": expected_shape)
    smoutput '    Actual shape:   ', (": act_shape)
  end.
)

NB. End suite and report
suite_end=: 3 : 0
  total=. TEST_PASS_COUNT + TEST_FAIL_COUNT
  smoutput '--- Suite Summary: ', CURRENT_SUITE_NAME, ' ---'
  smoutput '    Passed: ', (": TEST_PASS_COUNT), '/', (": total)
  if. TEST_FAIL_COUNT > 0 do.
    smoutput '    FAILED: ', (": TEST_FAIL_COUNT), ' tests'
    exit 1
  else.
    smoutput '    ALL TESTS PASSED'
  end.
)

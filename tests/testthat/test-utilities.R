# Test utility functions

test_that("Utility functions work correctly", {
  # Test null coalescing operator
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(5 %||% 10, 5)
})

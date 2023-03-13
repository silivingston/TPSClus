
test_that("TPSfit fits splines correctly", {
  expect_equal(TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
                  ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)$GAMsfitted$centered_x[1], -0.636163231)
})

test_that("TPSfit coefs are correct", {
  expect_equal(TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
                      ID="SubjectID", knots_time=c(0, 91, 182, 273, 365), n_fit_times=10)$GAMscoef$GAMbasis[1], -0.146563048)
})

test_that("Error should be produced if no knots provided", {
  expect_error(TPSfit(TS.sim, vars=c("Var1", "Var2", "Var3"), time="Time",
                      ID="SubjectID", n_fit_times=10))
})


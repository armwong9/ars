#testing the overall function
library(numDeriv)

test_that("Invalid input", {
  x <- 3
  
  expect_error(ars(x))
  
  #have some non log concave input?
  
  x <- function(x) {
    dgamma(x, 2)
  } #?
  
  
})

test_that("Correct input",{
  x <- dnorm
  num_samples <- 23
  
  expect_length(ars(x, n = num_samples), num_samples)
  expect_type(ars(x), 'double')
  
  #how do we want to test the output?
  
  x <- function(x) {
    dgamma(x, shape=5)
  }
  
  num_samples <- 12
  
  expect_length(ars(x, n = num_samples, bounds=c(0, Inf), starting_points = c(1,10)), num_samples)
  expect_type(ars(x, n = num_samples, bounds=c(0, Inf), starting_points = c(1,10)), 'double')
  
})


test_that("get the right distribution", {
  
  # normal
  expect_true(ks.test(ars(dnorm, n=1000), rnorm(1000))$p.value > 0.05)

  # gamma
  gam <- function(x) dgamma(x, shape=5)
  expect_true(ks.test(ars(gam, n=1000, bounds = c(0, Inf)), rgamma(1000, shape = 5))$p.value > 0.05)

  # beta
  bet <- function(x) dbeta(x,2,2)
  expect_true(ks.test(ars(bet, n=1000, bounds = c(0, 1)), rbeta(1000, 2,2))$p.value > 0.05)
  
  
})



#testing the module that do anything complicated 

#These below are some ideas

#sample_s 
#calc_functions or the 3 find functions?



test_that("Correct input for starting points",{
  
  good <- c(-1,1)
  bad <- c(1,-1)
  
  bounds <- c(-Inf, Inf)
  
  density <- function(x){
    log(dnorm(x))
  }
  
  expect_equal(get_StartingPoints(good, density, bounds), good)
  expect_error(get_StartingPoints(bad, density, bounds), "starting points are not on both sides of max value of function, choose better start points
                  or leave blank")
  
})


test_that("Unit test for find_Z",{
  
  f_vals <- c(1,1)
  f_deriv_vals <- c(1,-1)
  T_vals <- c(-1,1)
  
  expect_error(find_Z('b', 'b', 'b'))
  expect_error(find_Z(f_vals, 'b', 'b'))
  expect_error(find_Z(f_vals, f_deriv_vals, 'b'))
  expect_error(find_Z(f_vals, f_deriv_vals, 1))
  
  expect_equal(find_Z(f_vals, f_deriv_vals, T_vals), 0)
  
})


test_that("Unit test for find_U",{
  
  x <- 2
  f_vals <- c(1,1)
  f_deriv_vals <- c(1,-1)
  T_vals <- c(-1,1)
  Z_vals <- 0
  bounds <- c(-Inf, Inf)
  
  expect_error(find_U('b', 'b', 'b', 'b', 'b'))
  expect_error(find_U(x, 'b', 'b', 'b', 'b'))
  expect_error(find_U(x, f_vals, 'b', 'b', 'b', 'b'))
  expect_error(find_U(x, f_vals, f_deriv_vals, 'b', 'b', 'b'))
  expect_error(find_U(x, f_vals, f_deriv_vals, Z_vals, 'b', 'b'))
  expect_error(find_U(x, f_vals, f_deriv_vals, Z_vals, T_vals, 'b'))
  expect_error(find_U(x, f_vals, f_deriv_vals, Z_vals, 5, bounds))
  
  expect_equal(find_U(x, f_vals, f_deriv_vals, Z_vals, T_vals, bounds), 0)
  
})

test_that("Unit test for find_L",{
  
  x <- 2
  f_vals <- c(1,1)
  T_vals <- c(-1,1)
  bounds <- c(-Inf, Inf)
  
  expect_error(find_L('b', 'b', 'b', 'b'))
  expect_error(find_L(x, f_vals, 'b', 'b'))
  expect_error(find_L(c, f_vals, f_deriv_vals, 'b'))
  expect_error(find_L(x, f_vals, f_deriv_vals, 5))
  
  expect_equal(find_L(0, f_vals, T_vals, bounds), 1)

  
})



test_that("Unit test for initialize",{
  

  f <- dnorm
  starting_points <- c(-1,1)
  bounds <- c(-Inf, Inf)
  
  expect_error(initialize('b', 'b', 'b'))
  expect_error(initialize(dnorm, 'b', 'b'))
  expect_error(initialize(dnorm, starting_points, 'b'))
  expect_error(initialize(dnorm, starting_points, 5))
  
})



test_that("Unit test for update_step",{
  
  x <- 2
  t_vals <- c(-1,1)
  f <- dnorm
  f_der <- function(x) grad(f, x)
  bounds <- c(-Inf, Inf)
  k <- 1
  
  expect_error(update_step('b', 'b', 'b', 'b', 'b', 'b'))
  expect_error(update_step(x, 'b', 'b', 'b', 'b', 'b'))
  expect_error(update_step(x, t_vals, 'b', 'b', 'b', 'b'))
  expect_error(update_step(x, t_vals, f, 'b', 'b', 'b'))
  expect_error(update_step(x, t_vals, f, f_der, 'b', 'b'))
  expect_error(update_step(x, t_vals, f, f_der, bounds, 'b'))
  
})




test_that("Unit test for update_T",{
  
  x <- 2
  t_vals <- c(-1,1)
  k <- 1
  
  expect_error(update_T('b', 'b', 'b'))
  expect_error(update_T(x, 'b', 'b'))
  expect_error(update_T(x, t_vals, 'b'))

  
})




test_that("Unit test for sampling_step",{
  
  T_vals <- c(-1,1)
  f <- dnorm
  f_der <- function(x) grad(f, x)
  bounds <- c(-Inf, Inf)
  Z_vals <- 0
  u_vals <- c(1,1)
  f_vals <- c(1,1)
  f_deriv_vals <- c(1,-1)
  
  expect_error(sampling_step('b', 'b', 'b', 'b', 'b', 'b'))
  expect_error(sampling_step(T_vals, f, f_der, bounds, Z_vals, u_vals, 5, f_vals))
  
})


test_that("Unit test for squeezing_test",{
  
  w <- 0.01
  U_star <- 6
  L_star <- 4
  
  expect_error(squeezing_test('b', 'b', 'b'))
  expect_equal(squeezing_test(w,U_star,L_star)$accepted, TRUE)
  expect_equal(squeezing_test(w,U_star,L_star)$h_eval, FALSE)
  
})

test_that("Unit test for rejection_test",{
  
  w <- 0.01
  U_star <- 6
  h_star <- 4
  
  expect_error(rejection_test('b', 'b', 'b'))
  expect_equal(rejection_test(w,h_star,U_star), TRUE)
  
})




test_that("Unit test for calc_functions",{
  
  x <- 2
  t_vals <- c(-1,1)
  f <- dnorm
  f_der <- function(x) grad(f, x)
  bounds <- c(-Inf, Inf)

  
  expect_error(update_step('b', 'b', 'b', 'b', 'b'))
  expect_error(update_step(x, 'b', 'b', 'b', 'b'))
  expect_error(update_step(x, t_vals, 'b', 'b', 'b'))
  expect_error(update_step(x, t_vals, f, 'b', 'b'))
  expect_error(update_step(x, t_vals, f, f_der, 'b'))
  expect_error(update_step(x, t_vals, f, f_der, 'b'))
  
})








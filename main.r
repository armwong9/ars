library(numDeriv)
library(assertthat)

# Main function that does initialization, sampling, and updating

ars <- function(old_density_function, starting_points=NA, bounds=c(-Inf, Inf), n=20){
  # The overall function that computes all the steps in the Adaptive Rejection 
  # Sampling. The function takes in a density function and has the option to 
  # take in starting points, bounds and number of points to be sampled.
  # If the starting points and bounds are not specified the function sets
  # starting points and bounds. The default sampling size is 20 points. 
  # The function returns an array with n sampled points from the given
  # density function. 
  
  # Initialize
  
  
  density_function <- function(x){
    return(log(old_density_function(x)))
  }
  
  starting_points <- get_StartingPoints(starting_points, old_density_function, bounds)
  
  
  initialized <- initialize(density_function, starting_points, bounds)
  
  f_deriv <- initialized$f_deriv
  T_values <- initialized$T_values
  calculated_f <- initialized$calculated_f
  calculated_f_deriv <- initialized$calculated_f_deriv
  Z_values <- initialized$Z_values
  U_values <- initialized$U_values
  L_values <- initialized$L_values
  
  
  # Sample and Update steps
  
  # Empty vector to store accepted samples. 
  return_values <- vector(mode = "numeric", length = n)
  count <- 1
  
  # Sample and update whilst we haven't got enough values
  while(count < n+1){ 
    # Sampling step
    
    sampled_step <- sampling_step(T_values, density_function, f_deriv, bounds, 
                                  Z_values, U_values, calculated_f_deriv, 
                                  calculated_f) 
    
    x_star <- sampled_step$x_star
    h_eval <- sampled_step$h_eval
    accepted <- sampled_step$accepted
    
    if(accepted == TRUE){
      # If the sample was accepted then add x_star to the return values
      
      return_values[count] <- x_star
      
      count <- count + 1
      
    }
    
    # Update step if h(x*) and h'(x*) were evaluated in sampling
    if(h_eval == TRUE){
      
      updated_step <- update_step(x_star, T_values, density_function,
                                  f_deriv, bounds, length(T_values))
      
      T_values <- updated_step$T_values
      
      Z_values <- updated_step$Z_values
      U_values <- updated_step$U_values
      L_values <- updated_step$L_values
      calculated_f <- updated_step$calculated_f
      calculated_f_deriv <- updated_step$calculated_f_deriv
      
      
    } 
  }
  
  assert_that(length(return_values)==n)
  
  return(return_values)
  
}

get_StartingPoints <- function(starting_points, density, bounds){
  # A function that considers the starting points and asserts that
  # they are on either side of the density maximum.
  # If the function gets NA for starting points, it assumes that no starting
  # points were inputted in the beginning and then the function evaluates 
  # starting points to be used. 
  
  # Check input
  assert_that(is.function(density))
  assert_that(length(bounds) == 2, is.numeric(bounds))
  assert_that(bounds[1] < bounds[2])
  
  # If starting points given, check if on either side of density max.
  
  if(is.numeric(starting_points)){
    
    max_value <- optim(starting_points[1], function(x) -1*density(x), method="BFGS")$par
    
    assert_that( ((starting_points[1] < max_value) && (starting_points[2] > max_value)) , msg= 
                   "starting points are not on both sides of max value of function, choose better start points
                  or leave blank" )
    
    # Check output
    assert_that(length(starting_points) == 2 )
    assert_that(is.double(starting_points))
    assert_that(noNA(starting_points))
    
    return(starting_points)
    
  }
  
  # Set initial value for the optim function based on the bounds.
  
  
  init <- 1
  
  if(bounds[1] == -Inf){
    
    if (bounds[2] == Inf){
      
      init <- 1
      
    } else {
      
      init <- bounds[2]-5
      
    }
  } else {
    
    if (bounds[2] == Inf){
      
      init <- bounds[1]+5
      
    } else {
      
      boundMean <- (bounds[2]-bounds[1])/2
      init <- boundMean
      
    }
    
  }
  
  # Evaluate the starting points with the maximum value of
  # of the density function.
  
  max_value <- optim(init, function(x) -1*density(x), method="BFGS")$par
  
  starting_points <- c(max_value-0.2, max_value+0.2)
  
  # Check output
  assert_that(length(starting_points) == 2 )
  assert_that(is.double(starting_points))
  assert_that(noNA(starting_points))
  
  return(starting_points)
  
}

find_Z <- function(f_vals, f_deriv_vals, T_values){
  # A function that evaluates the z values, the intersects of the tangents
  # at xj and xj+1 (the T values). The function takes in the T values (x_j's)
  # and the values of the density function and differentiated density function 
  # at these x_j values. 
  
  #Check input
  assert_that(is.numeric(f_vals), is.numeric(f_deriv_vals), 
              is.numeric(T_values))
  assert_that( min(c(length(f_vals), length(f_deriv_vals), 
                     length(T_values))) ==
                 max(c(length(f_vals), length(f_deriv_vals), 
                       length(T_values))))
  
  
  n <- length(T_values)-1
  
  xj1hj1 <- T_values[2:(n+1)] * f_deriv_vals[2:(n+1)]
  xjhj <- T_values[1:n] * f_deriv_vals[1:n]
  
  num <- f_vals[2:(n+1)] - f_vals[1:n] - xj1hj1 + xjhj

  den <- f_deriv_vals[1:n]-f_deriv_vals[2:(n+1)]
  
  # Check that we are not dividing with zero. 
  assert_that(all(den != 0))
  
  Z_values <- num/den
  
  #Check output
  assert_that(is.double(Z_values))
  assert_that(length(Z_values) == (length(T_values)-1))
  
  return(Z_values)
  
}

# note that at T_values L and U are just evaluations of f

find_U <- function(x, f_vals, f_deriv_vals, Z_values, T_values, bounds){
  
  # A function that calculates the upper bound (u) for given T values 
  # and returns the estimates of the upper bounds for given x value/-s
  # The function takes where the upper bound should be estiamted (x), 
  # the T values (x_j's) and the values of the density function and 
  # differentiated density function at these T values. It also takes in the
  # Z values (intersect points of the tangets) and the given bounds. 
  
  #Check input
  assert_that(is.double(x), is.double(f_vals), is.double(f_deriv_vals), 
              is.double(Z_values), is.double(T_values),
              is.double(bounds))
  
  assert_that( min(c(length(f_vals), length(f_deriv_vals), 
                     length(T_values))) ==
                 max(c(length(f_vals), length(f_deriv_vals), 
                       length(T_values))))
  
  assert_that(length(x)<=length(T_values))
  assert_that(length(Z_values) == (length(T_values)-1))
  assert_that(length(bounds) == 2, is.double(bounds))
  assert_that(bounds[1] < bounds[2])
  
  # If already evaluated.
  if(x %in% T_values){
    u <- f_vals[match(x,T_values)]
    
    # Check output
    assert_that(is.double(u))
    
    return(u)
  }
  
  ind <- findInterval(x, c(bounds[1], Z_values, bounds[2]))
  
  u <- f_vals[ind] + ((x-T_values[ind])*f_deriv_vals[ind])
  
  # Check output
  #assert_that(length(u) == length(T_values))
  assert_that(is.double(u))
  
  return(u)
  
}


find_L <- function(x, f_vals, T_values, bounds){
  
  # A function that calculates the lower bound (u) for given T values 
  # and returns the estimates of the upper bounds for given x value/-s
  # The function takes where the upper bound should be estiamted (x), 
  # the T values (x_j's) and the values of the density function 
  #  at these T values and the given bounds. 
  
  #Check input
  assert_that(is.double(x), is.double(f_vals),  
              is.double(T_values), is.double(bounds))
  assert_that( min(c(length(f_vals), length(T_values))) ==
                 max(c(length(f_vals), length(T_values))))
  
  assert_that(length(bounds) == 2, is.double(bounds))
  assert_that(bounds[1] < bounds[2])
  
  # If already evaluated.
  
  if(x %in% T_values){
    l <- f_vals[match(x,T_values)]
    
    # Check output
    assert_that(is.double(l))
    
    return(l)
  }
  
  ind <- findInterval(x, c(bounds[1], T_values, bounds[2]))
  
  
  if ( (ind==1) || (ind==(length(T_values)+1)) ){
    return(-Inf)
  }
  
  
  term_1 <- (T_values[ind]-x)*f_vals[ind-1]
  term_2 <- (x-T_values[ind-1])*f_vals[ind]
  term_3 <- (T_values[ind]-T_values[ind-1])
  
  assert_that(term_3 != 0)
  
  L_values <- (term_1+term_2)/term_3
  
  #Check output
  assert_that(is.double(L_values))
  
  return(L_values)
  
}


sample_s <- function(lb, ub, Z_values, U_values, T_values, calculated_f_deriv, 
                     calculated_f){
  # A function that samples a value from a piecewise exponential distribution
  # with given intersections of tangent lines (z) , upper bound values (u), 
  # T values and returns that sample. 
  # The function takes in these given values as well as the values of the 
  # density function and differentiated density function at the given T values. 
  # The function also takes in a lower bound (lb) and upper bound (ub). 
  
  # Check input
  
  assert_that(is.double(lb), is.double(ub), is.double(Z_values), 
              is.double(U_values), is.double(calculated_f_deriv), 
              is.double(calculated_f))
  assert_that(noNA(lb), noNA(ub), noNA(Z_values), noNA(U_values),
              noNA(calculated_f_deriv), noNA(calculated_f))
  assert_that(lb <= ub, length(lb) == 1, length(ub) == 1) 
  
  # Set bounds of integration.
  
  bounds <- c(lb, Z_values, ub)
  
  # Number of starting points
  pt <- length(bounds)-1
  
  
  # Vector to store each piecewise cdf evaluation
  cdfs <- rep(0, pt)
  
  for(i in seq(1, pt)){
    # Evaluate the integral
    
    up <- function(x){ 
      exp(calculated_f[i] + (x-T_values[i])*calculated_f_deriv[i]) 
    }
    integral <- tryCatch(integrate(up, lower=bounds[i], upper=bounds[i+1])$value, error=function(e) e, warning=function(w) w)
    if (!is.numeric(integral)){
      stop("Please choose appropriate bounds for your function", call. = FALSE)
    } else {
      cdfs[i] <- integral
    }
    
  }
  
  
  # Normalising constant
  normc <- sum(cdfs)
  
  # assert that not dividing by 0
  assert_that(normc != 0)
  
  cdfs <- cdfs/normc
  
  # Choose which cdf interval to sample from 
  cdf_idx <- sample(length(cdfs), prob = cdfs, size=1)
  U <- runif(1)*cdfs[cdf_idx]
  
  if (calculated_f_deriv[cdf_idx] == 0){
    
    # Sample
    y <- bounds[cdf_idx] + U*(bounds[cdf_idx+1] - bounds[cdf_idx])
  
  } else {
    
    num <- (calculated_f_deriv[cdf_idx]*U / exp(calculated_f[cdf_idx])) + exp((bounds[cdf_idx] - T_values[cdf_idx])*calculated_f_deriv[cdf_idx])
    # Using Taylor expansion of log(x) in case we are taking negative log is negative
    check <- ifelse((num < 0), 
                    (num-1) - (1/2)*(num-1)^2 + (1/3)*(num-1)^3, log(num))
    # Sample
    y <- T_values[cdf_idx] + check/calculated_f_deriv[cdf_idx]
  
  }

  # Check output
  assert_that(is.double(y), length(y) == 1)
  
  return(y)
  
}


initialize <- function(f, starting_points, bounds){
  
  # A initialization function that initializes T and calculates U, L, 
  # and S at the starting values for given bounds and density function. 
  # The function takes in a density function f, starting points and bounds.
  # The function returns the derivative of the density function, 
  # the T values (sorted starting points), the intersection points of
  # the tangent lines (z), upper bound values (u), and lower bound values (l).
  # Also returns the values of the density function and differentiated 
  # density function at the given T values.
  
  #Check input
  assert_that(is.function(f))
  assert_that(is.double(starting_points))
  assert_that(length(bounds) == 2, is.numeric(bounds))
  assert_that(bounds[1] < bounds[2])
  
  
  # Get function for calculating derivative at a certain point
  f_deriv <- function(x) {
    result <- tryCatch({
      der <- grad(f, x)
      return(der)
    }, error = function(err) {
      print(x)
      h <- 0.0001 #.Machine$double.eps^(1/4)
      num <- f(x+h) - f(x)
      print(num)
      den <- h
      der2 <- num/den
      return(der2)
    })
  }

  
  # Sort T to be in ascending order
  T_values <- sort(starting_points)
  
  # Construct, uk, lk, zk and evaluate f and f'
  
  calced_functions <- calc_functions(T_values, T_values, f, f_deriv, bounds)
  
  Z_values <- calced_functions$Z_values
  U_values <- calced_functions$U_values
  L_values <- calced_functions$L_values
  calculated_f <- calced_functions$calculated_f
  calculated_f_deriv <- calced_functions$calculated_f_deriv
  
  # Check outputs 
  
  assert_that(is.double(Z_values), is.double(U_values), is.double(L_values),
              is.double(T_values), is.double(bounds), 
              is.double(calculated_f_deriv), is.double(calculated_f))
  assert_that( min(c(length(T_values), length(U_values), 
                     length(L_values), length(calculated_f_deriv), 
                     length(calculated_f))) ==
                 max(c(length(T_values), length(U_values), 
                       length(L_values), length(calculated_f_deriv), 
                       length(calculated_f))))
  assert_that(length(Z_values) == (length(T_values)-1))
  assert_that(noNA(Z_values), noNA(U_values), noNA(L_values),
              noNA(T_values), noNA(bounds), 
              noNA(calculated_f_deriv), noNA(calculated_f))
  
  return(list("f_deriv" = f_deriv, 
              "T_values" = T_values, 
              "calculated_f" = calculated_f,
              "calculated_f_deriv" = calculated_f_deriv,
              "Z_values" = Z_values,
              "U_values" = U_values,
              "L_values" = L_values))
  
}



update_step <- function(x_star, T_values, f, f_deriv, bounds, k){
  # A function that performs the updating step. 
  # The function takes in a sampled x, the T values, density function,
  # the differentiated density function, bounds and k. 
  # The function returns the new k value and T values.
  
  # Check inputs
  assert_that(is.function(f), is.function(f_deriv))
  assert_that(length(bounds) == 2, is.double(bounds), noNA(bounds))
  assert_that(bounds[1] < bounds[2])
  assert_that(is.double(x_star), is.double(T_values), is.numeric(k))
  assert_that(noNA(T_values))
  assert_that(length(x_star) == 1, length(k) == 1, length(T_values) == k)
  # Ensure log concavity
  #assert_that(hessian(f, x_star)<=0)
  
  # Update T values
  
  T_values <- update_T(x_star, T_values, k)
  
  # Construct, uk+1, lk+1, zk+1 and evaluate f and f'
  
  calced_functions <- calc_functions(T_values, T_values, f, f_deriv, bounds)
  
  Z_values <- calced_functions$Z_values
  U_values <- calced_functions$U_values
  L_values <- calced_functions$L_values
  calculated_f <- calced_functions$calculated_f
  calculated_f_deriv <- calced_functions$calculated_f_deriv
  
  k <- k + 1 
  
  # Check outputs 
  assert_that(is.numeric(k), length(k) == 1)
  assert_that( min(c(length(T_values), length(U_values), 
                     length(L_values), length(calculated_f_deriv), 
                     length(calculated_f))) ==
                 max(c(length(T_values), length(U_values), 
                       length(L_values), length(calculated_f_deriv), 
                       length(calculated_f))))
  
  return(list("Z_values" = Z_values,
              "U_values" = U_values,
              "L_values" = L_values,
              "k" = k,
              "T_values" = T_values,
              "calculated_f" = calculated_f,
              "calculated_f_deriv" = calculated_f_deriv
  ))
  
}

update_T <- function(x_star, T_values , k){
  # A function that adds the accepted sample to
  # the T values and sorts the T values. 
  # The function takes in a sampled value, T values and k.
  # The function returns the updated T_k+1 values. 
  
  # Checking inputs
  
  assert_that(is.double(x_star), length(x_star) == 1)
  assert_that(is.integer(k), length(k) == 1)
  assert_that(is.double(T_values), length(T_values) == k)
  
  # Add sample to T values
  T_values[k+1] <- x_star
  
  # Relabel elements of Tk+1 in ascending order
  T_values <- sort(T_values)
  
  # Check output
  assert_that(is.numeric(T_values), length(T_values) == (k+1))
  
  return(T_values)
  
}



sampling_step <- function(T_values, f, f_deriv, bounds, Z_values, 
                          U_values, calculated_f_deriv, calculated_f){
  # A function that performs the sampling step by sampling
  # and then performs the squeezing test and rejection
  # test if relevant. The function returns the sample,
  # whether the sample was accepted and 
  # whether h, h' were evaluated. 
  
  # Check inputs
  assert_that(is.double(Z_values), is.double(U_values), 
              is.double(T_values), is.double(bounds), 
              is.double(calculated_f_deriv), is.double(calculated_f))
  assert_that(is.function(f), is.function(f_deriv))
  assert_that( min(c(length(T_values), length(U_values), 
                     length(calculated_f_deriv), length(calculated_f))) ==
                 max(c(length(T_values), length(U_values), 
                       length(calculated_f_deriv), length(calculated_f))))
  assert_that(length(Z_values) == (length(T_values)-1))
  assert_that(length(bounds) == 2)
  assert_that(bounds[1] < bounds[2])
  
  # Sample from s_k equation
  x_star <- sample_s(bounds[1], bounds[2], Z_values, U_values, T_values, 
                     calculated_f_deriv, calculated_f)
  
  # Sample from UNIF[0,1]
  w <- runif(1)
  
  # Evaluate the equations for T_k and the sampled value
  
  calced_functions <- calc_functions(x_star, T_values, f, f_deriv, bounds)
  
  U_star <- calced_functions$U_values
  L_star <- calced_functions$L_values
  h_star <- f(x_star)
  
  # Perform squeezing test
  
  squeezing_tested <- squeezing_test(w, U_star, L_star)
  h_eval <- squeezing_tested$h_eval
  accepted <- squeezing_tested$accepted
  
  if(accepted == FALSE){
    # If squeezing test not passed then do rejection_test
    
    accepted <- rejection_test(w, x_star, U_star)
    
  }
  
  # Check output
  assert_that(is.logical(accepted), length(accepted) == 1)
  assert_that(is.logical(h_eval), length(h_eval) == 1)
  assert_that(is.numeric(x_star), length(x_star) == 1)
  
  return(list("h_eval" = h_eval,
              "accepted" = accepted,
              "x_star" = x_star))
  
}

squeezing_test <- function(w, U_star, L_star){
  # A function that performs the squeezing test
  # for a given w, value of u_k equation for the sample
  # and given value of l_k equation for the sample
  
  # Check input 
  assert_that(is.numeric(w), is.numeric(U_star), is.numeric(L_star))
  assert_that(length(w) == 1, length(U_star) == 1, length(L_star) == 1)
  assert_that(U_star >= L_star)
  assert_that(w >= 0, w <= 1)
  
  tol <- exp(L_star - U_star)
  
  # Check if test is passed
  if( w <= tol ){
    #if test if passed then sample accepted and h not evaluated  
    
    accepted <- TRUE
    h_eval <- FALSE
    
  } else {
    # If test if not passed then sample not accepted and h will be evaluated 
    
    accepted <- FALSE
    h_eval <- TRUE
    
  }
  
  # Check output
  assert_that(is.logical(accepted), length(accepted) == 1)
  assert_that(is.logical(h_eval), length(h_eval) == 1)
  
  return(list("accepted" = accepted,
              "h_eval" = h_eval))
  
}

rejection_test <- function(w, h_star, U_star){
  # A function that performs the rejection test
  # for a given w, value of h for the sample
  # and given value of u_k equation for the sample
  
  # Check input 
  assert_that(is.numeric(w), is.numeric(h_star), is.numeric(U_star))
  assert_that(length(w) == 1, length(h_star) == 1, length(U_star) == 1)
  assert_that(w >= 0, w <= 1)
  
  tol <- exp(h_star - U_star)
  
  if(w <= tol){
    # If test is passed then accept sample 
    
    accepted <- TRUE
    
  }else{
    # Otherwise don't accept sample
    
    accepted <- FALSE
    
  }
  
  # Check output
  assert_that(is.logical(accepted), length(accepted) == 1)
  
  return(accepted)
  
}

calc_functions <- function(x, T_values, f, f_deriv, bounds){
  # A function that evaluates Z and equations for u_k, l_k and s_k
  # for given T_k values and returns the corresponding values of the
  # equations for given x values. 
  
  # Check inputs
  assert_that(is.double(x), is.double(T_values), is.double(bounds))
  assert_that(noNA(T_values), noNA(bounds))
  assert_that(is.function(f), is.function(f_deriv))
  assert_that(length(x) <= length(T_values), length(bounds) == 2)
  assert_that(bounds[1] < bounds[2])
  
  # Evaluate f and its derivative at the starting values
  calculated_f <- sapply(T_values, f)
  calculated_f_deriv <- sapply(T_values, f_deriv)
  
  # Find intersection points
  Z_values <- find_Z(calculated_f, calculated_f_deriv, T_values)
  
  # Evaluate U,L,S at starting values
  U_values <- sapply(x, find_U, calculated_f, calculated_f_deriv, Z_values, T_values, bounds)
  L_values <- sapply(x, find_L, calculated_f, T_values, bounds)  
  
  # Check output
  assert_that(length(U_values) == length(x), length(L_values) == length(x))
  assert_that(length(calculated_f_deriv) == length(calculated_f))
  assert_that(length(calculated_f_deriv) == length(T_values))
  #assert_that(length(U_values) == length(x), 
  # length(L_values) == length(x), length(calculated_f) == length(x), 
  #length(calculated_f_deriv) == length(x))
  #assert_that(all(U_values >= 0), all(L_values >= 0), all(Z_values >= 0))
  assert_that(is.double(Z_values), is.double(U_values), is.double(L_values),
              is.double(T_values), is.double(bounds), 
              is.double(calculated_f_deriv), is.double(calculated_f))
  assert_that(length(Z_values) == (length(T_values)-1))
  
  assert_that(noNA(Z_values), noNA(U_values), noNA(L_values),
              noNA(T_values), noNA(bounds), 
              noNA(calculated_f_deriv), noNA(calculated_f))
  
  return(list("Z_values" = Z_values,
              "U_values" = U_values,
              "L_values" = L_values,
              "calculated_f" = calculated_f,
              "calculated_f_deriv" = calculated_f_deriv))
  
}

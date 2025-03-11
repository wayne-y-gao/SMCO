# Run Comparisons of SMCO with Other Optimization Algorithms

# "Optimization via Strategic Law of Large Numbers"
# By: Xiaohong Chen, Zengjing Chen, Wayne Yuan Gao, Xiaodong Yan, Guodong Zhang, and Yu Zhang
# Date: March 10, 2025
# GitHub Repository Maintained by: Wayne Yuan Gao
# For quesions, comments, and bug reports, contact: waynegao@upenn.edu

# This Rscript:
# (i) runs comparisons of SMCO algorithm with a List of other optimization algorithms
# (ii) can be used to replicate the numerical experiments in Section 4 of the paper:
#      "Optimization via Strategic Law of Large Numbers" 

# Required files:
# (1) "SMCO.R":
# (2) "testfuncs.R": contains a range of deterministic test function configurations
# (3) "testfunc_NN1L.R": contains the random test function configurations based on neural networks
# (4) "GD.R, "SignGD.R", "SPSA.R": R implementations of these algorithms
# (5) R packages: nloptr, GenSA, DEoptim, GA, PSO, parma, nloptr, optimg, stats

# Core functions: 
# (1) "run_comparison": core function that runs the comparison among the algorithms and test configurations
#                       as specified at the end of this Rscript
# (2) "print_results": print summary statistics to txt file and save R objects to Rdata file
# (3)  code at the end of this Rscript that specifies and calls run_comparisons and print_results

# Helper functions:
# (5) "generate_start_points": generate multiple starting points using the Latin Hypercube Sampling procedure,
#                              to be used by SMCO and local optimizers GD, SignGD, SPSA, parma/ADAM, optim/L-BFGS, nloptr/BOBYQA 

##########################
# Generate Starting Points
generate_start_points <- function(n_starts, bounds_lower, bounds_upper, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  d <- length(bounds_lower)
  points <- matrix(0, nrow = n_starts, ncol = d)
  
  # Pre-calculate segments
  segments <- seq(0, 1, length.out = n_starts + 1)
  lower_segs <- segments[-length(segments)]
  upper_segs <- segments[-1]
  
  # Latin Hypercube Sampling
  for (i in 1:d) {
    # Vectorized sampling from segments
    points[, i] <- lower_segs + runif(n_starts) * (upper_segs - lower_segs)
    # Scale to bounds
    points[, i] <- bounds_lower[i] + points[, i] * (bounds_upper[i] - bounds_lower[i])
  }
  
  return(points)
}

# Main Function to Run Comparisons of Different Optimization Algorithms

run_comparison <- function(name_config, dim_config, opt_config, to_maximize, n_replications, algo_names, truth_known, SMCO_Rname, SMCO_options) {
  
  source("testfuncs.R")
  source("testfunc_NN1L.R")
  source(SMCO_Rname)
  #library(foreach)
  #library(doParallel)
  
  # # Detect the number of available cores
  # num_cores <- detectCores()
  # # Use n-1 cores to leave one core free for the operating system
  # cluster <- makeCluster(num_cores - 1)
  # # Register the cluster for parallel execution
  # registerDoParallel(cluster)
  # 
  test_config <- assign_config(name = name_config, d = dim_config, opt_config = opt_config)
  dim_config <- test_config$dim
  
  # Random shifting of the origin and the lower/upper bounds
  bounds_lower = test_config$bounds_lower
  bounds_upper = test_config$bounds_upper
  
  bounds_diff = bounds_upper  - bounds_lower
  origin_shift <- rep(0, test_config$d) # bounds_diff * rnorm(test_config$d, 0, 1) #rep(0,test_config$d)
  # bounds_lower = bounds_lower + origin_shift
  # bounds_upper = bounds_upper + origin_shift + runif(test_config$d, 0.025, 0.05) * bounds_diff 
  
  test_func <- function(x) test_config$f(x - origin_shift)        
  #test_func <- test_config$f        
  test_func_minus <- function(x) (-test_func(x))
  
  # Assign properly signed function to be used in optimization algo
  # Add function call tracker

  if (to_maximize) { 
    test_func_for_max <- function(x) {return(test_func(x))} # for algo that maximizes functions
    test_func_for_min <- function(x) {return(test_func_minus(x))} # for algo that minimizes functions
  } else {
    test_func_for_max <- function(x) {return(test_func_minus(x))}
    test_func_for_min <- function(x) {return(test_func(x))}
  } 
  
  if (truth_known) {
    true_opt <- if (to_maximize) test_config$f_max else test_config$f_min
  } else {
    true_opt <- NA
  }
  
  # bounds_lower = test_config$bounds_lower
  # bounds_upper = test_config$bounds_upper
  
  n_algo <- length(algo_names)
  
  time_algo <- matrix(NA, n_algo, n_replications)
  fopt_algo = matrix(NA, n_algo, n_replications)
  fcalls_algo = matrix(NA, n_algo, n_replications)
  
  for (i in 1:n_replications) {
    
    start_points <- generate_start_points(SMCO_options$n_starts, bounds_lower, bounds_upper, seed = NULL)
    
    if ("SMCO" %in% algo_names) {
      
      index_algo = match("SMCO", algo_names)
      
      start_time = Sys.time()
      results <- SMCO_multi(f = test_func_for_max,  
                            bounds_lower, bounds_upper, start_points = start_points, 
                           opt_control = list(
                             bounds_buffer = SMCO_options$bounds_buffer,
                             buffer_rand = SMCO_options$buffer_rand,
                             n_starts = SMCO_options$n_starts,
                             iter_max = SMCO_options$iter_max,
                             iter_nstart = 1,
                             iter_boost = 0,
                             tol_conv = SMCO_options$tol_conv,
                             refine_search = F,
                             refine_ratio = 0,
                             partial_option = SMCO_options$partial_option,
                             use_runmax = F,
                             use_parallel = SMCO_options$use_parallel))
      
      end_time = Sys.time()
      time_algo[index_algo , i] = difftime(end_time, start_time, units = 'secs')
      fopt = results$best_result$f_optimal
      fopt_algo[index_algo,i] = ifelse(to_maximize, fopt, -fopt)
    }
    
    if ("SMCO_R" %in% algo_names) {
      
      source(SMCO_Rname)
      index_algo = match("SMCO_R", algo_names)
      
      start_time = Sys.time()
      results <- SMCO_multi(f = test_func_for_max, 
                            bounds_lower, bounds_upper, start_points = start_points, 
                            opt_control = list(
                              bounds_buffer = SMCO_options$bounds_buffer,
                              buffer_rand = SMCO_options$buffer_rand,
                              n_starts = SMCO_options$n_starts,
                              iter_max = SMCO_options$iter_max,
                              iter_nstart = 1,
                              iter_boost = 0,
                              tol_conv = SMCO_options$tol_conv,
                              refine_search = T,
                              refine_ratio = 0.5,
                              partial_option = "center",
                              use_runmax = T,
                              use_parallel = SMCO_options$use_parallel))
      
      end_time = Sys.time()
      
      time_algo[index_algo , i] = difftime(end_time, start_time, units = 'secs')
      fopt = results$best_result$f_optimal
      fopt_algo[index_algo,i] = ifelse(to_maximize, fopt, -fopt)
    }
    
    if ("SMCO_BR" %in% algo_names) {
      
      
      index_algo = match("SMCO_BR", algo_names)
      
      start_time = Sys.time()
      results <- SMCO_multi(f = test_func_for_max, 
                            bounds_lower, bounds_upper, start_points = start_points, 
                            opt_control = list(
                              bounds_buffer = SMCO_options$bounds_buffer,
                              buffer_rand = SMCO_options$buffer_rand,
                              n_starts = SMCO_options$n_starts,
                              iter_max = SMCO_options$iter_max/2,
                              iter_nstart = 1,
                              iter_boost = 99,
                              tol_conv = SMCO_options$tol_conv,
                              refine_search = T,
                              refine_ratio = 0.5,
                              partial_option = SMCO_options$partial_option,
                              use_runmax = T,
                              use_parallel = SMCO_options$use_parallel))
      
      end_time = Sys.time()
      
      time_algo[index_algo , i] = difftime(end_time, start_time, units = 'secs')
      fopt = results$best_result$f_optimal
      fopt_algo[index_algo,i] = ifelse(to_maximize, fopt, -fopt)
      
    }
    
    if ("GenSA" %in% algo_names) {
      library(GenSA)
      index_algo = match("GenSA", algo_names)
      
      start_time = Sys.time()
      
      GenSA_result = GenSA(
        par = start_points[1,],
        fn = test_func_for_min,
        lower = bounds_lower,
        upper = bounds_upper
      )
      end_time = Sys.time()
      
      time_algo[index_algo, i] = difftime(end_time, start_time, units = 'secs')
      fopt = GenSA_result$value
      fopt_algo[index_algo, i] = ifelse(to_maximize, -fopt, fopt)
    }
    
    if ("DEoptim" %in% algo_names) {
      library(DEoptim)
      index_algo = match("DEoptim", algo_names)
      
      start_time = Sys.time()
      DEoptim_result = DEoptim(fn = test_func_for_min, bounds_lower, bounds_upper, 
                               control = DEoptim.control(trace = FALSE))
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = DEoptim_result$optim$bestval
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }

    if ("CMAES" %in% algo_names) {
      library(parma)
      index_algo = match("CMAES", algo_names)
      
      start_time = Sys.time()
      CMAES_result = cmaes(par = start_points[1,],
                            fun = test_func_for_min,
                            lower = bounds_lower, upper = bounds_upper)
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = CMAES_result$objective
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }
        
    if ("GA" %in% algo_names) {
      library(GA)
      index_algo = match("GA", algo_names)
      
      start_time = Sys.time()
      GA_result <- ga(type = "real-valued",
                     fitness = test_func_for_max,
                     lower = bounds_lower, upper = bounds_upper)
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = GA_result@fitnessValue
      fopt_algo[index_algo,i] = ifelse(to_maximize, fopt, -fopt)
    }    

    if ("stogo" %in% algo_names) {
      library(nloptr)
      index_algo= match("stogo", algo_names)
      
      start_time = Sys.time()
      stogo_result = stogo(x0 = start_points[1,],
                             fn = test_func_for_min,
                             lower = bounds_lower, upper = bounds_upper)
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = stogo_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }
    
    if ("PSO" %in% algo_names) {
      library(pso)
      index_algo = match("PSO", algo_names)
      
      start_time = Sys.time()
      PSO_result = psoptim(par = start_points[1,],
                           fn = test_func_for_min,
                           lower = bounds_lower, upper = bounds_upper)
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = PSO_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }   
 
    if ("GD" %in% algo_names) {
      source("GD.R")
      index_algo = match("GD", algo_names)
      
      start_time = Sys.time()
      results <- apply(start_points, 1, function(start_point) {
        gradient_descent(x0 = start_point,
               f = test_func_for_min,
               bounds_lower = bounds_lower,
               bounds_upper = bounds_upper)}
      )
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]      
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }

    if ("SignGD" %in% algo_names) {
      source("SignGD.R")
      index_algo = match("SignGD", algo_names)
      
      start_time = Sys.time()
      results <- apply(start_points, 1, function(start_point) {
        signGD(x0 = start_point,
               f = test_func_for_min,
               bounds_lower = bounds_lower,
               bounds_upper = bounds_upper)}
      )
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]      
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }  
    
    if  ((!to_maximize) && ("ADAM" %in% algo_names)) {
      library(optimg)
      index_algo = match("ADAM", algo_names)
      
      start_time = Sys.time()
      results <- apply(start_points, 1, function(start_point) {
        optimg(par = start_point,
               fn = test_func_for_min,
               method = "ADAM")}
      )
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]      
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }    
    
    if ("SPSA" %in% algo_names) {
      source("SPSA.R")
      index_algo = match("SPSA", algo_names)
      
      start_time = Sys.time()
      
      results <- apply(start_points, 1, function(start_point) {
        spsa(f = test_func_for_max,
             A = 50,
             a = 0.10,
             c = 1e-3,
             bounds_lower, bounds_upper, 
             startpoint = start_point
        )})
      f_values <- sapply(results, function(x) x$f_optimal)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$f_optimal
      fopt_algo[index_algo,i] = ifelse(to_maximize, fopt, -fopt)
      
    }
    
    
    if ("optimNM" %in% algo_names) {
      
      index_algo = match("optimNM", algo_names)
      start_time = Sys.time()
      
      results <- apply(start_points, 1, function(start_point) {
                       optim(par = start_point, 
                             fn = test_func_for_min,
                             method = "Nelder-Mead",
                             lower = bounds_lower,
                             upper = bounds_upper)}
      )
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)  
    }
    
    
    if ("optimLBFGS" %in% algo_names) {
      index_algo = match("optimLBFGS", algo_names)
      
      start_time = Sys.time()
      
      results <- apply(start_points, 1, function(start_point) {
        optim(par = start_point,
              fn = test_func_for_min,
              method = "L-BFGS-B",
              lower = bounds_lower,
              upper = bounds_upper
        )})
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]
      
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)  
    }
    
    
    if ("optimSANN" %in% algo_names) {
      index_algo = match("optimSANN", algo_names)
      
      start_time = Sys.time()
      results <- apply(start_points, 1, function(start_point) {
        optim(par = start_point, 
              fn = test_func_for_min,
              method = "SANN", 
              control = list(temp = 500, 
                             maxit = 0.4*1e+6, 
                             reltol = 1e-7))}
      )
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }
    
    
    
    if ("nloptr" %in% algo_names) {
      library(nloptr)
      index_algo = match("nloptr", algo_names)
      
      start_time = Sys.time()
      results <- apply(start_points, 1, function(start_point) {
        lbfgs(x0 = start_point,
               fn = test_func_for_min,
               lower = bounds_lower, upper = bounds_upper)})
      
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
      
    }
    
    if ("BOBYQA" %in% algo_names) {
      library(nloptr)
      index_algo = match("BOBYQA", algo_names)
      
      start_time = Sys.time()
      results <- apply(start_points, 1, function(start_point) {
        bobyqa(x0 = start_point,
               fn = test_func_for_min,
               lower = bounds_lower, upper = bounds_upper)})
      f_values <- sapply(results, function(x) x$value)
      best_idx <- which.max(f_values)
      best_result <- results[[best_idx]]      
      end_time = Sys.time()
      
      time_algo[index_algo,i] = difftime(end_time, start_time, units = 'secs')
      fopt = best_result$value
      fopt_algo[index_algo,i] = ifelse(to_maximize, -fopt, fopt)
    }
    
       
  }
  
  if (truth_known == F ||  is.na(true_opt)) {
    if (to_maximize) {
      fopt_algo[is.infinite(fopt_algo)] = NA
      best_opt = max(fopt_algo, na.rm = T)
    } else {
      fopt_algo[is.infinite(fopt_algo)] = NA
      best_opt = min(fopt_algo, na.rm = T)
    }
  } else {
    best_opt = true_opt
  }
  
  time_avg = rowMeans(time_algo)
  fcalls_avg = rowMeans(fcalls_algo)
  AE = abs(fopt_algo - best_opt)
  rMSE = sqrt(rowMeans(AE^2))
  MAE = rowMeans(AE)
  fopt_mean = rowMeans(fopt_algo)
  fopt_sd = apply(fopt_algo, 1, sd)
  if (best_opt == 0) { 
    rMSEpct = NA
    goodpct = rowMeans((AE < 1e-2)) 
  } else {
    rMSEpct = rMSE/abs(best_opt)
    goodpct = rowMeans((AE/abs(best_opt) < 1e-2))
  }
  
  sum_results <- data.frame(
    algo_names = algo_names,
    time_avg = time_avg, #round(time_avg, digits = 6),
    # fcalls_avg = fcalls_avg,
    rMSE = rMSE, #round(rMSE, digits = 6),
    rMSEpct = rMSEpct, #round(rMSEpct, digits = 6),
    MAE = MAE, #round(MAE, digits = 6),
    goodpct = goodpct, #round(goodpct, digits = 6),
    fopt_mean = fopt_mean, #round(fopt_mean, digits = 6),
    fopt_sd = fopt_sd  #round(fopt_sd, digits = 6)
  )
  
  return(list(name_config = name_config,
              dim_config  = dim_config,
              bounds_lower = bounds_lower,
              bounds_upper = bounds_upper,
              to_maximize = to_maximize,
              truth_known = truth_known,
              true_opt = true_opt,
              n_replications = n_replications,
              SMCO_options = SMCO_options,
              sum_results = sum_results,
              time_algo = time_algo,
              fopt_algo = fopt_algo))
  
}

print_results <- function(comp_results) {
  
  library(xtable)
  
  with(comp_results, {
    
    name_append = ""
    current_time <- format(Sys.time(), "%y%m%d_%H%M%S")
    name_to_maximize <- if (to_maximize) "Max" else "Min"
    
    filename <- paste0("Comp_", name_config, "_", dim_config, "d", name_to_maximize, name_append, ".txt")
    sink(filename)
    
    cat("Test Configuration:\n") 
    cat("\t", paste0(name_config, " ", dim_config, "d\n"))
    if (name_config == "ANN") {
      cat("\n")
      cat("\t", paste0("Input Dim: ", opt_config$input_dim, "\n"))
      cat("\t", paste0("Hidden Dim: ", opt_config$hidden_dim, "\n"))
      cat("\t", paste0("Noise Sigma: ", opt_config$noise_sigma, "\n"))
      cat("\t", paste0("Use True Reg Func for MSE: ", opt_config$use_truth, "\n"))
      cat("\t", paste0("# of Simulated Input Points: ", opt_config$n, "\n"))
      cat("\t", paste0("search range: on [", bounds_lower[1], ", ", bounds_upper[1], "] hypercube\n"))
    } else {
      cat("\t", paste0("on [", bounds_lower[1], ", ", bounds_upper[1], "] hypercube"), "\n")
      cat("\t", "True Optimum Known? ", paste0(truth_known), "\n")
    }
    cat("\n")
    cat("\t", "# Replications: ", n_replications, "\n")
    cat("\t", "Maximization? ", to_maximize, "\n")
    cat("\n")
    
    cat("SMCO Configuration:\n")
    cat("\t", "Rscript: ", paste(SMCO_Rname), "\n")
    cat("\t", "# Starting Points: ", SMCO_options$n_starts, "\n") # " ~ sqrt(dim_config) * 10\n")
    cat("\t", "iter_max: ", paste(SMCO_options$iter_max), "\n")
    cat("\t", "iter_boost: ", paste(SMCO_options$iter_boost), "\n")
    cat("\t", "bounds_buffer: ", paste(SMCO_options$bounds_buffer), "\n")
    cat("\t", "buffer_rand: ", paste(SMCO_options$buffer_rand), "\n")
    cat("\t", "tol_conv: ", paste(SMCO_options$tol_conv), "\n")
    cat("\t", "refine_search: ", paste(SMCO_options$refine_search), "\n")
    cat("\t", "refine_ratio: ", paste(SMCO_options$refine_ratio), "\n")
    cat("\t", "use_runmax: ", paste(SMCO_options$use_runmax), "\n")
    cat("\t", "use_parallel: ", paste(SMCO_options$use_parallel), "\n")
    cat("\n")
    
    if (truth_known == T && !is.na(true_opt)) {
      cat("True Optimal Value: ", true_opt, "\n")
      cat("\n")
    } else {
      cat("True Optimal Value: ", "Unknown", "\n")
      cat("\n")
    }
    
    print(sum_results)
    
    cat("\n")
    sink()
    # 
    # texname = paste0("Comp_", name_config, "_", dim_config, "d", name_to_maximize, name_append, "_tex.txt")
    # 
    # print_to_tex <- data.frame(rMSE = sum_results$rMSE, 
    #                            rMSEpct = sum_results$rMSEpct, 
    #                            time_avg = sum_results$time_avg)
    # 
    # row.names(print_to_tex) <- sum_results$algo_names 
    # 
    # sink(texname)
    # print(xtable(print_to_tex, display = c("s","e", "e", "f"), digits = 2))
    # sink()
    
    Rdataname <- paste0("Comp_", name_config, "_", dim_config, "d", name_to_maximize, name_append, ".RData")
    save.image(file = Rdataname)
    
  })
}

# Specify and run comparisons
if (T) { # Set to T to run comparisons

  # SMCO Configuration
  SMCO_options <- list(
    seed = NULL,
    bounds_buffer = 0.05,
    buffer_rand = F,
    tol_conv = 1e-8,
    #refine_search = T,         # "run_comparison"
    #refine_ratio = 0.5,        # "run_comparison"
    partial_option = "center",
    iter_max = 200,
    #iter_nstart = 100,         # Overridden in "run_comparison"
    #iter_boost = 0,            # Overridden in "run_comparison"
    #n_starts = 100,            # Overridden in dimension-adapted formula below
    #use_runmax = T,            # Overridden in "run_comparison"
    use_parallel = F
  )
  
  SMCO_Rname = "SMCO.R"
  
  # Specification for neural network test case only
  name_config <- "ANN"
  opt_config <- list(
    input_dim = 3,
    hidden_dim = 5,
    noise_sigma = 1,
    n = 1000,
    use_truth = T,
    noise_sigma = 1,
    seed = NULL
  )
  
  # Main comparison specifications:
  
  # Set number of replications
  n_replications = 100
  
  # Choose list of test functions
  name_config_list <- c("Rastrigin", "Qing") 
  # See testfuncs.R for complete list of test functions incorporated
  # c("ANN", Rastrigin", Ackley", "Griewank", "Michalewicz", "DixonPrice", #"Dropwave", "Shubert", "McCormick", "Zakharov", "DixonPrice", "Rosenbrock" , "EggHolder",...)
  
  # Set list of dimensions for test functions (if applicable)
  dim_config_list <- 10 #c(2, 10, 20, 50, 100, 200, 1000)
  
  # Set list of maximization/minimization 
  to_maximize_list <- F  # c(T,F)
  
  # Choose whether to calculate RMSE relative to known true optimal value (if applicable)
  truth_known = F
  
  # Set list of optimization algorithms to be included in the comparison
  algo_names =  c("SMCO", "SMCO_R", "SMCO_BR") #, "GD", "SignGD", "ADAM", "SPSA", "optimLBFGS", "BOBYQA", "GenSA", "DEoptim", "CMAES", "stogo", "GA", "PSO") 
    # List of algorithms  
    # c("SMCO", "SMCO_R", "SMCO_BR", "GD", "SignGD", "ADAM", "SPSA", "optimLBFGS", "BOBYQA", "GenSA", "DEoptim", "CMAES", "stogo", "GA", "PSO") 
  
  config_grid <- expand.grid(to_maximize = to_maximize_list, name_config = name_config_list, dim_config = dim_config_list)
  
  for (i in 1:nrow(config_grid)) {
    
    set.seed(123)
    
    name_config = config_grid$name_config[i]
    dim_config = config_grid$dim_config[i]
    to_maximize = config_grid$to_maximize[i]
    SMCO_options$n_starts = min(100, round(sqrt(dim_config) * 10, digits = 0))
    
    comp_results <- run_comparison(name_config = name_config, dim_config = dim_config, opt_config = opt_config, to_maximize = to_maximize,
                                 truth_known = truth_known, algo_names = algo_names, n_replications = n_replications,  
                                 SMCO_Rname =SMCO_Rname, SMCO_options = SMCO_options)
  
    print_results(comp_results)
  }
}

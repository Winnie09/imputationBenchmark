## Call python code from R
##
## Copyright Lingxue Zhu (lzhu@cmu.edu).
## All Rights Reserved.

## Note:
## If your python is installed under /usr/local/bin, 
## you may need to update the system path first by 
# Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"),sep=":"))


data_to_csv <- function(data, data_dir, file_prefix, varname) {
  res = list()
  if (!is.null(data)){
    filename = paste0(data_dir, "/", file_prefix, "_", varname, ".csv")
    write.table(data, file=filename, 
                quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
    res[varname] = filename
  }
  return(res)
}

PyGEM <- function(
                  ## Change this to your local path
                  py_script="/Users/lingxue/Documents/Thesis/SingleCell/scUnif/scUnif.py",
                  
                  BKexpr=NULL, ## sample-by-gene
                  K=3, 
                  SCexpr=NULL, ## sample-by-gene
                  G=NULL, ## cell-type info for single cell
                  data_dir = "data/", ## directory to hold data files
                  data_prefix = "", ## prefix added to data files
                  
                  ## model parameters
                  iMarkers = NULL,
                  init_A=NULL,  min_A=1e-6,
                  init_alpha=NULL, est_alpha=TRUE,
                  init_pkappa=NULL, init_ptau=NULL, ## mean and variance 
                  burnin=20, sample=20, thin=1, ## for SC Gibbs sampling
                  bk_mean_approx=TRUE, ## whether use mean-gibbs updates for BK e-step
                  MLE_CONV=1e-3, EM_CONV=1e-3, 
                  MLE_maxiter=1, EM_maxiter=2,
                  verbose=1,
                  out_dir="out/", ## output directory
                  output_prefix="out_", ## output prefix
                  log_dir="log/" ## logging directory
                  ) {
  
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive=TRUE)
  }
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive=TRUE)
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
  }
  
  #########################
  ## write expression data to file
  ## and record file names
  #########################
  arguments = list()

  ## note that cell type and gene index in python are 0-indexed
  if (!is.null(G)) {
    G = G-1
  }
  if (!is.null(iMarkers)) {
    iMarkers = iMarkers - 1
  }
  
  arguments = c(arguments,
                data_to_csv(BKexpr, data_dir, data_prefix, "bulk_expr_file"),
                data_to_csv(SCexpr, data_dir, data_prefix, "single_cell_expr_file"),
                ## note that cell type and gene index in python are 0-indexed
                data_to_csv(G, data_dir, data_prefix, "single_cell_type_file"),
                data_to_csv(iMarkers, data_dir, data_prefix, "iMarkers_file"),
                data_to_csv(init_A, data_dir, data_prefix, "initial_A_file"),
                data_to_csv(init_alpha, data_dir, data_prefix, "initial_alpha_file"))
  
  ################################
  ## other algorithm parameters
  #################################
  arguments = c(arguments,
              list(output_directory=out_dir, 
                   output_prefix=output_prefix, 
                   logging_file=paste0(log_dir, "/", output_prefix, ".log"),
                   EM_maxiter=EM_maxiter, Mstep_maxiter=MLE_maxiter,
                   EM_convergence_tol=EM_CONV, Mstep_convergence_tol=MLE_CONV,
                   gibbs_thinning=thin, gibbs_sample_number=sample, burn_in_length=burnin,
                   number_of_cell_types=K, 
                   mininimal_A=min_A,
                   verbose_level=verbose)
  )
  
  ## several parameters need special handling to have the right format
  if (!est_alpha){
    arguments$no_est_alpha = ""
  }
  if (bk_mean_approx) {
    arguments$mean_approx = ""
  } else {
    arguments$no_mean_approx = ""
  }
  if (!is.null(init_pkappa)) {
    arguments$initial_kappa_mean_var=paste(init_pkappa, collapse = " ")
  }
  if (!is.null(init_ptau)) {
    arguments$initial_tau_mean_var=paste(init_ptau, collapse = " ")
  }
  
  
  #########################
  ## run python
  #########################
  ## command line args
  cmdargs = ""
  for (i in 1:length(arguments)) {
    cmdargs = paste0(cmdargs, " --", names(arguments)[i], " ", arguments[i])
  }
  
  ## run
  arguments$R_runlog <- system(paste0("python ", py_script, cmdargs), intern=TRUE)
  return (arguments)
}
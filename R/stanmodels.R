# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("alt_g", "alt_gl", "marg_dr_npp_g4", "marg_dr_npp_gl4", "marg_dr_npp_glpknown4", "marg_dr_pp_g4", "marg_dr_pp_gl4", "marg_dr_pp_glpknown4", "marg_ndr_pp_g4", "marg_ndr_pp_gl4", "marg_ndr_pp_glpknown4")

# load each stan module
Rcpp::loadModule("stan_fit4alt_g_mod", what = TRUE)
Rcpp::loadModule("stan_fit4alt_gl_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_dr_npp_g4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_dr_npp_gl4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_dr_npp_glpknown4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_dr_pp_g4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_dr_pp_gl4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_dr_pp_glpknown4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_ndr_pp_g4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_ndr_pp_gl4_mod", what = TRUE)
Rcpp::loadModule("stan_fit4marg_ndr_pp_glpknown4_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})

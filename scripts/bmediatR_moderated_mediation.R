#' Definions and encodings for models included in the moderated mediation analysis through Bayesian model selection
#'
#' This function prints out a table describing the models considered in the moderated mediation analysis.
#'
#' @export
#' @examples moderation_model_info()
moderation_model_info <- function(){
  writeLines(c("likelihood models for all hypothesis",
               "hypotheses encoded by presence (1) or absence (0) of 'X->m, y->m, X->y' edges on the DAG",
               "(*) denotes reverse causation 'm<-y', i denotes moderation",
               "H1: '-,0,0' / y does not depend on X or m",
               "H2: '-,1,0' / y depends on m but not X",
               "H3: '-,0,1' / y depends on X but not m",
               "H4: '-,1,1' / y depends on X and m",
               "H5: '0,-,-' / m does not depend on X",
               "H6: '1,-,-' / m depends on X",
               "H7: '0,*,-' / m depends on y but not X",
               "H8: '1,*,-' / m depends on X and y",
               "H9: 'i,-,-' / m depends on X and is moderated by t",
               "H10: 'i,*,-' / m depends on X and y and relation between X and m is moderated by t",
               "H11: '0,*i,-' / m depends on y but not X and is moderated by t",
               "H12: '1,*i,-' / m depends on X and y and y -> m is moderated by t",
               "H13: '-,0,i' / y depends on X but not m and is moderated by t",
               "H14: '-,1,i' / y depends on X and m and X->y is moderated by t",
               "H15: '-,i,0' / y depends on m but not X and is moderated by t",
               "H16: '-,i,1' / y depends on X and m and m -> y is moderated by t",
               "all include covariates Z and t",
               "",
               "combinations of hypotheses for all cases",
               "cases encoded by presence (1) or absence (0) of 'X->m, m->y, X->y' edges on the DAG",
               "(*) denotes reverse causation 'm<-y' and (i) denotes moderation",
               "c1:  '0,0,0' / H1 and H5",
               "c2:  '0,1,0' / H2 and H5",
               "c3:  '1,0,0' / H1 and H6",
               "c4:  '1,1,0' / H2 and H6 - complete mediation",
               "c5:  '0,0,1' / H3 and H5",
               "c6:  '0,1,1' / H4 and H5",
               "c7:  '1,0,1' / H3 and H6 - colocalization",
               "c8:  '1,1,1' / H4 and H6 - partial mediation",
               "c9:  '0,*,0' / H1 and H7",
               "c10: '1,*,0' / H1 and H8",
               "c11: '0,*,1' / H3 and H7 - complete med (reactive)",
               "c12: '1,*,1' / H3 and H8",
               "c13: '0,i,0' / H5 and H15",
               "c14: 'i,0,0' / H1 and H9",
               "c15: 'i,1,0' / H2 and H9 - moderated mediation A",
               "c16: '1,i,0' / H15 and H6 - moderated mediation B",
               "c17: '0,0,i' / H13 and H5",
               "c18: '0,1,i' / H14 and H5",
               "c19: '0,i,1' / H16 and H5",
               "c20: '1,0,i' / H13 and H6",
               "c21: 'i,0,1' / H3 and H9",
               "c22: '1,1,i' / H14 and H6",
               "c23: 'i,1,1' / H4 and H9",
               "c24: '1,i,1' / H16 and H6",
               "c25: '0,*i,0' / H1 and H11",
               "c26: 'i,*,0' / H1 and H10",
               "c27: '1,*i,0' / H1 and H12",
               "c28: '0,*,i' / H13 and H7",
               "c29: '0,*i,1' / H3 and H11",
               "c30: '1,*,i' / H13 and H8",
               "c31: 'i,*,1' / H3 and H10",
               "c32: '1,*i,1' / H3 and H12"))
}


#' Summaries of posterior model probabilities function
#'
#' This function takes the log posterior probability of the data (posterior likelihood) for the various models, the log prior model probabilities, and
#' returns log posterior odds
#'
#' @param ln_prob_data Log posterior likelihoods under the various models, returned by bmediatR().
#' @param ln_prior_c Log prior model probabilities. If posterior_summary() is being used for a non-default posterior odds
#' summary, the log prior model probabilities used with bmediatR() are stored in its output.
#' @param c_numerator The index of models to be summed in the numerator of the posterior odds. Models and their order provided with
#' model_summary().
#' @export
#' @examples posterior_summary()
posterior_summary_moderation <- function(ln_prob_data,
                                        ln_prior_c,
                                        c_numerator){
  #function to compute log odds from log probabilities
  ln_odds <- function(ln_p, numerator){
    ln_odds_numerator <- apply(ln_p[,numerator,drop=F], 1, matrixStats::logSumExp)
    ln_odds_denominator <- apply(ln_p[,-numerator,drop=F], 1, matrixStats::logSumExp)
    ln_odds <- ln_odds_numerator -ln_odds_denominator
  }

  #ensure c_numerator is a list
  if (!is.list(c_numerator)){
    c_numerator <- list(c_numerator)
  }

  #presets for ln_prior_c;
  if (ln_prior_c[1]=="complete"){
    ln_prior_c <- c(rep(0,8), rep(-Inf,4), rep(0, 12), rep(-Inf, 8))
  } else if (ln_prior_c[1]=="partial"){
    ln_prior_c <- c(rep(-Inf,4), rep(0,4), rep(-Inf,4), rep(-Inf, 4), rep(0, 8), rep(-Inf,8))
  } else if (ln_prior_c[1]=="reactive"){
    ln_prior_c <- rep(0,32)
  }

  #ensure ln_prior_c sum to 1 on probability scale and that it is a matrix
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
    ln_prior_c <- matrix(ln_prior_c, nrow(ln_prob_data), length(ln_prior_c), byrow=T)
  }

  #compute posterior probabilities for all cases
  #cases encoded by presence (1) or absence (0) of 'X->m, m->y, X->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y' and (i) denotes moderation
  #c1:  '0,0,0' / H1 and H5
  #c2:  '0,1,0' / H2 and H5
  #c3:  '1,0,0' / H1 and H6
  #c4:  '1,1,0' / H2 and H6 - complete mediation
  #c5:  '0,0,1' / H3 and H5
  #c6:  '0,1,1' / H4 and H5
  #c7:  '1,0,1' / H3 and H6 - colocalization
  #c8:  '1,1,1' / H4 and H6 - partial mediation
  #c9:  '0,*,0' / H1 and H7
  #c10: '1,*,0' / H1 and H8
  #c11: '0,*,1' / H3 and H7 - Reactive
  #c12: '1,*,1' / H3 and H8
  #c13: '0,i,0' / H5 and H15
  #c14: 'i,0,0' / H1 and H9
  #c15: 'i,1,0' / H2 and H9
  #c16: '1,i,0' / H15 and H6
  #c17: '0,0,i' / H13 and H5
  #c18: '0,1,i' / H14 and H5
  #c19: '0,i,1' / H16 and H5
  #c20: '1,0,i' / H13 and H6
  #c21: 'i,0,1' / H3 and H9
  #c22: '1,1,i' / H14 and H6
  #c23: 'i,1,1' / H4 and H9
  #c24: '1,i,1' / H16 and H6
  #c25: '0,*i,0' / H1 and H11
  #c26: 'i,*,0' / H1 and H10
  #c27: '1,*i,0' / H1 and H12
  #c28: '0,*,i' / H13 and H7
  #c29: '0,*i,1' / H3 and H11
  #c30: '1,*,i' / H13 and H8
  #c31: 'i,*,1' / H3 and H10
  #c32: '1,*i,1' / H3 and H12

  ln_post_c <- cbind(ln_prob_data[,1] + ln_prob_data[,5] + ln_prior_c[,1],
                     ln_prob_data[,2] + ln_prob_data[,5] + ln_prior_c[,2],
                     ln_prob_data[,1] + ln_prob_data[,6] + ln_prior_c[,3],
                     ln_prob_data[,2] + ln_prob_data[,6] + ln_prior_c[,4],
                     ln_prob_data[,3] + ln_prob_data[,5] + ln_prior_c[,5],
                     ln_prob_data[,4] + ln_prob_data[,5] + ln_prior_c[,6],
                     ln_prob_data[,3] + ln_prob_data[,6] + ln_prior_c[,7],
                     ln_prob_data[,4] + ln_prob_data[,6] + ln_prior_c[,8],
                     ln_prob_data[,1] + ln_prob_data[,7] + ln_prior_c[,9],
                     ln_prob_data[,1] + ln_prob_data[,8] + ln_prior_c[,10],
                     ln_prob_data[,3] + ln_prob_data[,7] + ln_prior_c[,11],
                     ln_prob_data[,3] + ln_prob_data[,8] + ln_prior_c[,12],
                     ln_prob_data[,5] + ln_prob_data[,15] + ln_prior_c[,13],
                     ln_prob_data[,1] + ln_prob_data[,9] + ln_prior_c[,14],
                     ln_prob_data[,2] + ln_prob_data[,9] + ln_prior_c[,15],
                     ln_prob_data[,15] + ln_prob_data[,6] + ln_prior_c[,16],
                     ln_prob_data[,13] + ln_prob_data[,5] + ln_prior_c[,17],
                     ln_prob_data[,14] + ln_prob_data[,5] + ln_prior_c[,18],
                     ln_prob_data[,16] + ln_prob_data[,5] + ln_prior_c[,19],
                     ln_prob_data[,13] + ln_prob_data[,6] + ln_prior_c[,20],
                     ln_prob_data[,3] + ln_prob_data[,9] + ln_prior_c[,21],
                     ln_prob_data[,14] + ln_prob_data[,6] + ln_prior_c[,22],
                     ln_prob_data[,4] + ln_prob_data[,9] + ln_prior_c[,23],
                     ln_prob_data[,16] + ln_prob_data[,6] + ln_prior_c[,24],
                     ln_prob_data[,1] + ln_prob_data[,11] + ln_prior_c[,25],
                     ln_prob_data[,1] + ln_prob_data[,10] + ln_prior_c[,26],
                     ln_prob_data[,1] + ln_prob_data[,12] + ln_prior_c[,27],
                     ln_prob_data[,13] + ln_prob_data[,7] + ln_prior_c[,28],
                     ln_prob_data[,3] + ln_prob_data[,11] + ln_prior_c[,29],
                     ln_prob_data[,13] + ln_prob_data[,8] + ln_prior_c[,30],
                     ln_prob_data[,3] + ln_prob_data[,10] + ln_prior_c[,31],
                     ln_prob_data[,3] + ln_prob_data[,12] + ln_prior_c[,32]
  )

  ln_ml <- apply(ln_post_c, 1, matrixStats::logSumExp)
  ln_post_c <- ln_post_c - ln_ml

  colnames(ln_post_c) <- c("0,0,0",
                           "0,1,0",
                           "1,0,0",
                           "1,1,0",
                           "0,0,1",
                           "0,1,1",
                           "1,0,1",
                           "1,1,1",
                           "0,*,0",
                           "1,*,0",
                           "0,*,1",
                           "1,*,1",
                           "0,i,0",
                           "i,0,0",
                           "i,1,0",
                           "1,i,0",
                           "0,0,i",
                           "0,1,i",
                           "0,i,1",
                           "1,0,i",
                           "i,0,1",
                           "1,1,i",
                           "i,1,1",
                           "1,i,1",
                           "0,*i,0",
                           "i,*,0",
                           "1,*i,0",
                           "0,*,i",
                           "0,*i,1",
                           "1,*,i",
                           "i,*,1",
                           "1,*i,1")
  rownames(ln_post_c) <- rownames(ln_prob_data)

  #compute prior odds for each combination of cases
  ln_prior_odds <- sapply(c_numerator, ln_odds, ln_p=ln_prior_c)
  ln_prior_odds <- matrix(ln_prior_odds, ncol=length(c_numerator))
  rownames(ln_prior_odds) <- rownames(ln_post_c)

  #compute posterior odds for each combination of cases
  ln_post_odds <- sapply(c_numerator, ln_odds, ln_p=ln_post_c)
  ln_post_odds <- matrix(ln_post_odds, ncol=length(c_numerator))
  rownames(ln_post_odds) <- rownames(ln_post_c)

  if (is.null(c_numerator)) {
    colnames(ln_post_odds) <- colnames(ln_prior_odds) <- c_numerator
  } else {
    colnames(ln_post_odds) <- colnames(ln_prior_odds) <- names(c_numerator)
  }

  #return results
  list(ln_post_c=ln_post_c, ln_post_odds=ln_post_odds, ln_prior_odds=ln_prior_odds, ln_ml=ln_ml)
}

#' Column indeces for commonly used posterior odds
#'
#' This helper function returns the columns of the log posterior model probabilities to be summed for
#' commonly desired log posterior odds summaries.
#'
#' @param odds_type The desired posterior odds.
#' @export
#' @examples return_preset_odds_index()
return_preset_odds_index_moderation <- function(odds_type = c("mediation",
                                                              "partial",
                                                              "complete",
                                                              "colocal",
                                                              "y_depends_x",
                                                              "reactive",
                                                              "moderation",
                                                              "moderation_a",
                                                              "moderation_b",
                                                              "moderation_b_reactive",
                                                              "moderation_c")) {

  presets <- list("mediation" = c(4, 8),
                  "partial" = 8,
                  "complete" = 4,
                  "colocal" = 7,
                  "y_depends_x" = c(4:8, 11, 12),
                  "reactive" = 9:12,
                  "moderation" = 13:32,
                  "moderation_a" = c(14, 15, 21, 23, 26, 31),
                  "moderation_b" = c(13, 16, 19, 24),
                  "moderation_b_reactive" = c(25, 27, 29, 32),
                  "moderation_c" = c(17, 18, 20, 22, 28, 30)
                  )

  index_list <- presets[odds_type]
  index_list
}

## Function to process data and optionally align them
process_data_moderation <- function(y, M, X, t,
                         Z = NULL, Z_y = NULL, Z_M = NULL,
                         w = NULL, w_y = NULL, w_M = NULL,
                         align_data = TRUE,
                         verbose = TRUE) {

  # Ensure y is a vector
  if (is.matrix(y)) { y <- y[,1] }

  # Ensure X, M, t, Z, Z_y, and Z_M are matrices
  X <- as.matrix(X)
  M <- as.matrix(M)
  t <- as.matrix(t)
  if (!is.null(Z)) { Z <- as.matrix(Z) }
  if (!is.null(Z_y)) { Z_y <- as.matrix(Z_y) }
  if (!is.null(Z_M)) { Z_M <- as.matrix(Z_M) }

  # Process covariate matrices
  if (is.null(Z_y)) { Z_y <- matrix(NA, length(y), 0); rownames(Z_y) <- names(y) }
  if (is.null(Z_M)) { Z_M <- matrix(NA, nrow(M), 0); rownames(Z_M) <- rownames(M) }

  if (!is.null(Z)) {
    if (align_data) {
      Z_y <- cbind(Z, Z_y[rownames(Z),])
      Z_M <- cbind(Z, Z_M[rownames(Z),])
    } else {
      Z_y <- cbind(Z, Z_y)
      Z_M <- cbind(Z, Z_M)
    }
  }

  # Process weight vectors
  if (is.null(w)) {
    if (is.null(w_y)) { w_y <- rep(1, length(y)); names(w_y) <- names(y) }
    if (is.null(w_M)) { w_M <- rep(1, nrow(M)); names(w_M) <- rownames(M) }
  } else {
    w_y <- w_M <- w
  }

  if (align_data) {
    # M and Z_M can have NAs
    overlapping_samples <- Reduce(f = intersect, x = list(names(y),
                                                          rownames(X),
                                                          rownames(Z_y),
                                                          names(w_y),
                                                          rownames(t)))

    if (length(overlapping_samples) == 0 | !any(overlapping_samples %in% unique(c(rownames(M), rownames(Z_M), names(w_M))))) {
      stop("No samples overlap. Check rownames of M, X, t, Z (or Z_y and Z_M) and names of y and w (or w_y and w_M).", call. = FALSE)
    } else if (verbose) {
      writeLines(text = c("Number of overlapping samples:", length(overlapping_samples)))
    }

    # Ordering
    y <- y[overlapping_samples]
    M <- M[overlapping_samples,, drop = FALSE]
    X <- X[overlapping_samples,, drop = FALSE]
    t <- t[overlapping_samples,, drop = FALSE]
    Z_y <- Z_y[overlapping_samples,, drop = FALSE]
    Z_M <- Z_M[overlapping_samples,, drop = FALSE]
    w_y <- w_y[overlapping_samples]
    w_M <- w_M[overlapping_samples]
  }

  # Drop observations with missing y, X, or t and update n
  complete_y <- !is.na(y)
  complete_X <- !apply(is.na(X), 1, any)
  complete_t <- !apply(is.na(X), 1, any)

  complete <- complete_y & complete_X & complete_t

  y <- y[complete]
  M <- M[complete,, drop = FALSE]
  X <- X[complete,, drop = FALSE]
  t <- t[complete,, drop = FALSE]
  Z_y <- Z_y[complete,, drop = FALSE]
  Z_M <- Z_M[complete,, drop = FALSE]
  w_y <- w_y[complete]
  w_M <- w_M[complete]

  # Drop columns of Z_y and Z_M that are invariant
  Z_y_drop <- which(apply(Z_y, 2, function(x) var(x)) == 0)
  Z_M_drop <- which(apply(Z_M, 2, function(x) var(x)) == 0)
  if (length(Z_y_drop) > 0) {
    if (verbose) {
      writeLines(paste("Dropping invariants columns from Z_y:", colnames(Z_y)[Z_y_drop]))
    }
    Z_y <- Z_y[,-Z_y_drop, drop = FALSE]
  }
  if (length(Z_M_drop) > 0) {
    if (verbose) {
      writeLines(paste("Dropping invariants columns from Z_M:", colnames(Z_M)[Z_M_drop]))
    }
    Z_M <- Z_M[,-Z_M_drop, drop = FALSE]
  }

  # Return processed data
  list(y = y,
       M = M,
       X = X,
       t = t,
       Z_y = Z_y, Z_M = Z_M,
       w_y = w_y, w_M = w_M)
}

#' Convert preset options to log prior case probabilities function
#'
#' This function takes the log prior case probabilities, and if a preset is provided, converts it into the formal log prior case
#' probability.
#'
#' @param ln_prior_c Log prior case probabilities. If posterior_summary_moderation() is being used for a non-default posterior odds
#' summary, the log prior case probabilities used with bmediatR() are stored in its output.
#' @export
#' @examples return_ln_prior_c_from_presets_moderation()
return_ln_prior_c_from_presets_moderation <- function(ln_prior_c) {

  #presets for ln_prior_c;
  if (ln_prior_c[1]=="complete"){
    ln_prior_c <- c(rep(0,8), rep(-Inf,4), rep(0, 12), rep(-Inf, 8))
  } else if (ln_prior_c[1]=="partial"){
    ln_prior_c <- c(rep(-Inf,4), rep(0,4), rep(-Inf,4), rep(-Inf, 4), rep(0, 8), rep(-Inf,8))
  } else if (ln_prior_c[1]=="reactive"){
    ln_prior_c <- rep(0,32)
  } else {
    ln_prior_c <- ln_prior_c
  }
  ln_prior_c
}



#' Bayesian model selection for moderated mediation analysis function
#'
#' This function takes an outcome (y), candidate mediators (M), a candidate moderator (t), and a driver as a design matrix (X) to perform a
#' Bayesian model selection analysis for mediation.
#'
#' @param y Vector or single column matrix of an outcome variable. Single outcome variable expected.
#' Names or rownames must match across M, X, t, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported.
#' Names or rownames must match across y, X, t, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param X Design matrix of the driver. Names or rownames must match across y, M, t, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE.
#' If align_data = FALSE, dimensions and order must match across inputs. One common application is for X to represent genetic information at a QTL,
#' either as founder strain haplotypes or variant genotypes, though X is generalizable to other types of variables.
#' @param t Design matrix of the moderator for the outcome and mediator. Names or rownames must match across y, M, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE.
#' If align_data = FALSE, dimensions and order must match across inputs.
#' @param Z DEFAULT: NULL. Design matrix of covariates that influence the outcome and mediator variables.
#' Names or rownames must match to those of y, M, X, t, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data=FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param Z_y DEFAULT: NULL. Design matrix of covariates that influence the outcome variable.
#' Names or rownames must match to those of y, M, X, t, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param Z_M DEFAULT: NULL. Design matrix of covariates that influence the mediator variables.
#' Names or rownames must match across y, M, X, t, Z_y, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param w DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis that applies to both
#' y and M. Names must match across y, M, X, t, Z, Z_y, and Z_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where w
#' is a vector of the number of individuals per strain. If no w, w_y, or w_M is given, observations are equally weighted as 1s for y and M.
#' If w is provided, it supercedes w_y and w_M.
#' @param w_y DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of y. Names must match across y, M, X, t, Z, Z_y, Z_M, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_y is a vector of the number of individuals per strain used to
#' measure y. If no w_y (or w) is given, observations are equally weighted as 1s for y.
#' @param w_M DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of M. Names must match across y, M, X, t, Z, Z_y, Z_M, and w_y (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_M is a vector of the number of individuals per strain use to
#' measure M. If no w_M (or w) is given, observations are equally weighted as 1s for M.
#' @param tau_sq_mu DEFAULT: 1000. Variance component for the intercept. The DEFAULT represents a diffuse prior, analagous to
#' a fixed effect term.
#' @param tau_sq_Z DEFAULT: 1000. Variance component for the covariates encoded in Z. The DEFAULT represents a diffuse prior, analagous
#' to fixed effect terms.
#' @param phi_sq DEFAULT: c(1, 1, 1). Each element of (a, b, c) represents one of the relationships being evaluated for mediation,
#' specifically the ratio of signal to noise. a is the effect of X on M, b is the effect of M on y, and c is the effect of X on y.
#' The DEFAULT represents relationships that explain 50% of the variation in the outcome variable.
#' @param ln_prior_c DEFAULT: "complete". The prior log case probabilities. See model_info() for description of likelihoods and their
#' combinations into cases. Simplified pre-set options are available, including "complete", "partial", and "reactive".
#' @param align_data DEFAULT: TRUE. If TRUE, expect vector and matrix inputes to have names and rownames, respectively. The overlapping data
#' will then be aligned, allowing the user to not have to reduce data to overlapping samples and order them.
#' @export
#' @examples bmediatR()
bmediatR_moderation <- function(y, M, X, t,
                                Z = NULL, Z_y = NULL, Z_M = NULL,
                                w = NULL, w_y = NULL, w_M = NULL,
                                kappa = 0.001, lambda = 0.001,
                                tau_sq_mu = 1000,
                                tau_sq_Z = 1000,
                                tau_sq_t = 1000, # Prior for the marginal effect of t
                                phi_sq = c(1, 1, 1, 1, 1, 1), # prior for edges a, b, c, and the interaction term for edge a, b, and c
                                ln_prior_c = "complete",
                                options_X = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                                options_t = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                                align_data = TRUE,
                                verbose = TRUE) {

  #presets for ln_prior_c;
  ln_prior_c <- return_ln_prior_c_from_presets_moderation(ln_prior_c = ln_prior_c)

  #ensure ln_prior_c sum to 1 on probability scale
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
  }

  #optionally align data
  processed_data <- process_data_moderation(y = y, M = M, X = X, t= t,
                                 Z_y = Z_y, Z_M = Z_M,
                                 w_y = w_y, w_M = w_M,
                                 align_data = align_data,
                                 verbose = verbose)
  y <- processed_data$y
  M <- processed_data$M
  X <- processed_data$X
  t <- processed_data$t
  Z_y <- processed_data$Z_y
  Z_M <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_M <- processed_data$w_M

  #dimension of y
  n <- length(y)

  #dimension of Z's
  p_y <- ncol(Z_y)
  p_M <- ncol(Z_M)

  #scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p_y > 0) { Z_y <- apply(Z_y, 2, scale) }
  if (p_M > 0) { Z_M <- apply(Z_M, 2, scale) }

  #optionally use sum-to-zero contrast for X
  #recommended when X is a matrix of factors, with a column for every factor level
  if (options_X$sum_to_zero == TRUE) {
    C <- sumtozero_contrast(ncol(X))
    X <- X%*%C
  }

  #optionally use sum-to-zero contrast for t
  #recommended when t is a matrix of factors, with a column for every factor level
  if (options_t$sum_to_zero == TRUE) {
    C <- sumtozero_contrast(ncol(t))
    t <- t%*%C
  }

  #optionally center and scale X
  X <- apply(X, 2, scale, center = options_X$center, scale = options_X$scale)

  #optionally center and scale t
  t <- apply(t, 2, scale, center = options_t$center, scale = options_t$scale)

  #dimension of X and t
  d <- ncol(X)
  k <- ncol(t)

  # Interaction between X and moderator
  tX <- model.matrix(~ 0 + X:t)
  # Dimension of interaction
  ix <- ncol(tX)

  # Interaction between y and moderator (for reactive models)
  ty <- model.matrix(~ 0 + y:t)
  # Dimension of interaction
  iy <- ncol(ty)

  #column design matrix for mu
  ones <- matrix(1, n)

  #begin Bayesian calculations
  if (verbose) { print("Initializing", quote = FALSE) }

  #reformat priors
  kappa = rep(kappa, 16)
  lambda = rep(lambda, 16)
  tau_sq_mu = rep(tau_sq_mu, 16)
  tau_sq_Z = rep(tau_sq_Z, 16)
  tau_sq_t = rep(tau_sq_t, 16)
  phi_sq_X = c(NA, NA, phi_sq[3], phi_sq[3], NA, phi_sq[1], NA, phi_sq[1],
               phi_sq[1], phi_sq[1], NA, phi_sq[1], phi_sq[3], phi_sq[3], NA, phi_sq[3])
  phi_sq_m = c(NA,phi_sq[2],NA,phi_sq[2],NA,NA,NA,NA,
               NA, NA, NA, NA, NA, phi_sq[2], phi_sq[2], phi_sq[2])
  phi_sq_y = c(NA,NA,NA,NA,NA,NA,phi_sq[2],phi_sq[2],
               NA, phi_sq[2], phi_sq[2], phi_sq[2], NA, NA, NA, NA)
  phi_sq_int = c(rep(NA, 8), # effect size prior for interaction terms
                 phi_sq[5], phi_sq[5], phi_sq[6], phi_sq[6], phi_sq[4], phi_sq[4], phi_sq[6], phi_sq[6])

  #identify likelihoods that are not supported by the prior
  #will not compute cholesky or likelihood for these
  calc_ln_prob_data <- rep(NA, 16)
  calc_ln_prob_data[1] <- any(!is.infinite(ln_prior_c[c(1,3,9,10)]))
  calc_ln_prob_data[2] <- any(!is.infinite(ln_prior_c[c(2,4)]))
  calc_ln_prob_data[3] <- any(!is.infinite(ln_prior_c[c(5,7,11,12)]))
  calc_ln_prob_data[4] <- any(!is.infinite(ln_prior_c[c(6,8)]))
  calc_ln_prob_data[5] <- any(!is.infinite(ln_prior_c[c(1,2,5,6)]))
  calc_ln_prob_data[6] <- any(!is.infinite(ln_prior_c[c(3,4,7,8)]))
  calc_ln_prob_data[7] <- any(!is.infinite(ln_prior_c[c(9,11)]))
  calc_ln_prob_data[8] <- any(!is.infinite(ln_prior_c[c(10,12)]))
  calc_ln_prob_data[9] <- any(!is.infinite(ln_prior_c[c(14,15,21,23)]))
  calc_ln_prob_data[10] <- any(!is.infinite(ln_prior_c[c(26, 31)]))
  calc_ln_prob_data[11] <- any(!is.infinite(ln_prior_c[c(25,29)]))
  calc_ln_prob_data[12] <- any(!is.infinite(ln_prior_c[c(27,32)]))
  calc_ln_prob_data[13] <- any(!is.infinite(ln_prior_c[c(28,30)]))
  calc_ln_prob_data[14] <- any(!is.infinite(ln_prior_c[c(18,22)]))
  calc_ln_prob_data[15] <- any(!is.infinite(ln_prior_c[c(13,16)]))
  calc_ln_prob_data[16] <- any(!is.infinite(ln_prior_c[c(19,24)]))

  #likelihood models for all hypothesis
  #hypotheses encoded by presence (1) or absence (0) of 'X->m, y->m, X->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y', i denotes moderation
  #H1: '-,0,0' / y does not depend on X or m
  #H2: '-,1,0' / y depends on m but not X
  #H3: '-,0,1' / y depends on X but not m
  #H4: '-,1,1' / y depends on X and m
  #H5: '0,-,-' / m does not depend on X
  #H6: '1,-,-' / m depends on X
  #H7: '0,*,-' / m depends on y but not X
  #H8: '1,*,-' / m depends on X and y
  #H9: 'i,-,-' / m depends on X and is moderated by t
  #H10: 'i,*,-' / m depends on X and y and relation between X and M is moderated by t
  #H11: '0,*i,-' / m depends on y but not X and is moderated by t
  #H12: '1,*i,-' / m depends on X and y and y -> M is moderated by t
  #H13: '-,0,i' / y depends on X but not m and is moderated by t
  #H14: '-,1,i' / y depends on X and m and X->y is moderated by t
  #H15: '-,i,0' / y depends on m but not X and is moderated by t
  #H16: '-,i,1' / y depends on X and m and m -> y is moderated by t
  #all include covariates Z and t

  #design matrices for H1,H3,H5-H13 complete cases (do not depend on m)
  X1 <- cbind(ones, Z_y, t)
  X3 <- cbind(ones, X, Z_y, t)
  X5 <- cbind(ones, Z_M, t)
  X6 <- cbind(ones, X, Z_M, t)
  X7 <- cbind(X5, y)
  X8 <- cbind(X6, y)
  X9 <- cbind(X6, tX)
  X10 <- cbind(X8, tX)
  X11 <- cbind(X7, ty)
  X12 <- cbind(X8, ty)
  X13 <- cbind(X3, tX)

  #check if all scale hyperparameters are identical for H1 and H5
  #implies sigma1 and sigma5 identical, used to reduce computations
  sigma5_equal_sigma1 <- all(lambda[1]==lambda[5],
                             tau_sq_mu[1] == tau_sq_mu[5],
                             tau_sq_Z[1] == tau_sq_Z[5],
                             identical(Z_y, Z_M),
                             identical(w_y, w_M))

  #check if all scale hyperparameters are identical for H3 and H6
  #implies sigma3 and sigma6 identical, used to reduce computations
  sigma6_equal_sigma3 <- all(lambda[3]==lambda[6],
                             tau_sq_mu[3] == tau_sq_mu[6],
                             tau_sq_Z[3] == tau_sq_Z[6],
                             identical(Z_y, Z_M),
                             identical(w_y, w_M))

  #check if all scale hyperparameters are identical for H9 and H13
  #implies sigma3 and sigma6 identical, used to reduce computations
  sigma13_equal_sigma9 <- all(lambda[9]==lambda[13],
                              tau_sq_mu[9] == tau_sq_mu[13],
                              tau_sq_Z[9] == tau_sq_Z[13],
                              identical(Z_y, Z_M),
                              identical(w_y, w_M),
                              tau_sq_t[9] == tau_sq_t[13])


  #prior variance matrices (diagonal) for H1-H16
  v1 <- c(tau_sq_mu[1], rep(tau_sq_Z[1], p_y), rep(tau_sq_t[1], k))
  v2 <- c(tau_sq_mu[2], rep(tau_sq_Z[2], p_y), rep(tau_sq_t[2], k), phi_sq_m[2])
  v3 <- c(tau_sq_mu[3], rep(phi_sq_X[3], d), rep(tau_sq_Z[3], p_y), rep(tau_sq_t[3], k))
  v4 <- c(tau_sq_mu[4], rep(phi_sq_X[4], d), rep(tau_sq_Z[4], p_y), rep(tau_sq_t[4], k), phi_sq_m[4])
  v7 <- c(tau_sq_mu[7], rep(tau_sq_Z[7], p_M), rep(tau_sq_t[7], k), phi_sq_y[7])
  v8 <- c(tau_sq_mu[8], rep(phi_sq_X[8], d), rep(tau_sq_Z[8], p_M), rep(tau_sq_t[8], k), phi_sq_y[8])
  v9 <- c(tau_sq_mu[9], rep(phi_sq_X[9], d), rep(tau_sq_Z[9], p_M), rep(tau_sq_t[9], k), rep(phi_sq_int[9], ix))
  v10 <- c(tau_sq_mu[10], rep(phi_sq_X[10], d), rep(tau_sq_Z[10], p_M), rep(tau_sq_t[10], k), rep(phi_sq_int[10], ix))
  v11 <- c(tau_sq_mu[11], rep(tau_sq_Z[11], p_M), rep(tau_sq_t[11], k), phi_sq_y[11], rep(phi_sq_int[11], iy))
  v12 <- c(tau_sq_mu[12], rep(phi_sq_X[12], d), rep(tau_sq_Z[12], p_M), rep(tau_sq_t[12], k), rep(phi_sq_int[12], iy))
  v14 <- c(tau_sq_mu[14], rep(phi_sq_X[14], d), rep(tau_sq_Z[14], p_y), rep(tau_sq_t[14], k), phi_sq_m[14], rep(phi_sq_int[14], ix))


  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]){
    v5 <- c(tau_sq_mu[5], rep(tau_sq_Z[5], p_M), rep(tau_sq_t[5], k))
  }

  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]){
    v6 <- c(tau_sq_mu[6], rep(phi_sq_X[6], d), rep(tau_sq_Z[6], p_M), rep(tau_sq_t[5], k))
  }

  if (!sigma13_equal_sigma9 | !calc_ln_prob_data[9]){
    v13 <- c(tau_sq_mu[13], rep(phi_sq_X[13], d), rep(tau_sq_Z[13], p_y), rep(tau_sq_t[13], k), rep(phi_sq_int[13], ix))
  }

  #scale matrices for H1,H3,H5-H8 complete cases, H9-H13 (do not depend on m)
  sigma1 <- crossprod(sqrt(lambda[1]*v1)*t(X1))
  sigma3 <- crossprod(sqrt(lambda[3]*v3)*t(X3))
  sigma7 <- crossprod(sqrt(lambda[7]*v7)*t(X7))
  sigma8 <- crossprod(sqrt(lambda[8]*v8)*t(X8))
  sigma9 <- crossprod(sqrt(lambda[9]*v9)*t(X9))
  sigma10 <- crossprod(sqrt(lambda[10]*v10)*t(X10))
  sigma11 <- crossprod(sqrt(lambda[11]*v11)*t(X11))
  sigma12 <- crossprod(sqrt(lambda[12]*v12)*t(X12))

  diag(sigma1) <- diag(sigma1) + lambda[1]/w_y
  diag(sigma3) <- diag(sigma3) + lambda[3]/w_y
  diag(sigma7) <- diag(sigma7) + lambda[7]/w_M
  diag(sigma8) <- diag(sigma8) + lambda[8]/w_M
  diag(sigma9) <- diag(sigma9) + lambda[9]/w_M
  diag(sigma10) <- diag(sigma10) + lambda[10]/w_M
  diag(sigma11) <- diag(sigma11) + lambda[11]/w_M
  diag(sigma12) <- diag(sigma12) + lambda[12]/w_M

  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]){
    sigma5 <- crossprod(sqrt(lambda[5]*v5)*t(X5))
    diag(sigma5) <- diag(sigma5) + lambda[5]/w_M
  }

  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]){
    sigma6 <- crossprod(sqrt(lambda[6]*v6)*t(X6))
    diag(sigma6) <- diag(sigma6) + lambda[6]/w_M
  }

  if (!sigma13_equal_sigma9 | !calc_ln_prob_data[9]){
    sigma13 <- crossprod(sqrt(lambda[13]*v13)*t(X13))
    diag(sigma13) <- diag(sigma13) + lambda[13]/w_y
  }

  #object to store likelihoods
  ln_prob_data=matrix(-Inf, ncol(M), 16)
  rownames(ln_prob_data) <- colnames(M)
  colnames(ln_prob_data) <- c("-,0,0",
                              "-,1,0",
                              "-,0,1",
                              "-,1,1",
                              "0,-,-",
                              "1,-,-",
                              "0,*,-",
                              "1,*,-",
                              "i,-,-",
                              "i,*,-",
                              "0,*i,-",
                              "1,*i,-",
                              "-,0,i",
                              "-,1,i",
                              "-,i,0",
                              "-,i,1")

  #identify batches of M that have the same pattern of missing values
  missing_m <- bmediatR:::batch_cols(M)

  #iterate over batches of M with same pattern of missing values
  if (verbose) { print("Iterating", quote = FALSE) }
  counter <- 0

  for (b in 1:length(missing_m)){
    #subset to non-missing observations
    index <- rep(T, length(y))
    index[missing_m[[b]]$omit] <- FALSE

    if (any(index)){
      y_subset <- y[index]
      w_y_subset <- w_y[index]
      w_M_subset <- w_M[index]
      t_subset <- t[index,,drop = F]

      #cholesky matrices for H1,H3,H5-H8 non-missing observations (do not depend on m)
      if (calc_ln_prob_data[1]){sigma1_chol_subset <- chol(sigma1[index,index])}
      if (calc_ln_prob_data[3]){sigma3_chol_subset <- chol(sigma3[index,index])}
      if (calc_ln_prob_data[7]){sigma7_chol_subset <- chol(sigma7[index,index])}
      if (calc_ln_prob_data[8]){sigma8_chol_subset <- chol(sigma8[index,index])}
      if (calc_ln_prob_data[9]){sigma9_chol_subset <- chol(sigma9[index,index])}
      if (calc_ln_prob_data[10]){sigma10_chol_subset <- chol(sigma10[index,index])}
      if (calc_ln_prob_data[11]){sigma11_chol_subset <- chol(sigma11[index,index])}
      if (calc_ln_prob_data[12]){sigma12_chol_subset <- chol(sigma12[index,index])}

      if (sigma5_equal_sigma1 & calc_ln_prob_data[1]) {
        sigma5_chol_subset <- sigma1_chol_subset
      } else if (calc_ln_prob_data[5]) {
        sigma5_chol_subset <- chol(sigma5[index,index])
      }

      if (sigma6_equal_sigma3 & calc_ln_prob_data[3]) {
        sigma6_chol_subset <- sigma3_chol_subset
      } else if (calc_ln_prob_data[6]) {
        sigma6_chol_subset <- chol(sigma6[index,index])
      }

      if (sigma13_equal_sigma9 & calc_ln_prob_data[9]){
        sigma13_chol_subset <- sigma9_chol_subset
      } else if (calc_ln_prob_data[13]){
        sigma13_chol_subset <- chol(sigma13[index,index])
      }

      #compute H1 and H3 outside of the mediator loop (invariant)
      if (calc_ln_prob_data[1]) { ln_prob_data1 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma1_chol_subset, df = kappa[1]) }
      if (calc_ln_prob_data[3]) { ln_prob_data3 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma3_chol_subset, df = kappa[3]) }
      if (calc_ln_prob_data[13]){ln_prob_data13 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma13_chol_subset, df = kappa[13])}

      #iterate over mediators
      for (i in missing_m[[b]]$cols) {
        counter <- counter + 1
        if (counter%%1000==0 & verbose) { print(paste(counter, "of", ncol(M)), quote=F) }

        #set current mediator non-missing observations
        m_subset <- M[index,i]

        # Calculate design matrix for interaction of M and t
        tM <- model.matrix(~ 0 + m_subset:t_subset)
        im <- ncol(tM)

        #design matrix for H2 and H4 non-missing observations
        X2_subset <- cbind(X1[index,,drop = FALSE], m_subset)
        X4_subset <- cbind(X3[index,,drop = FALSE], m_subset)
        X14_subset <- cbind(X4_subset, tX[index,,drop=FALSE])
        X15_subset <- cbind(X2_subset, tM)
        X16_subset <- cbind(X4_subset, tM)

        v15 <- c(tau_sq_mu[15], rep(tau_sq_Z[15], p_y), rep(tau_sq_t[15], k), phi_sq_m[15], rep(phi_sq_int[15], im))
        v16 <- c(tau_sq_mu[16], rep(phi_sq_X[16], d), rep(tau_sq_Z[16], p_y), rep(tau_sq_t[16], k), phi_sq_m[16], rep(phi_sq_int[16], im))

        #scale and cholesky matrices for H2 and H4 non-missing observations
        sigma2_subset <- crossprod(sqrt(lambda[2]*v2)*t(X2_subset))
        sigma4_subset <- crossprod(sqrt(lambda[4]*v4)*t(X4_subset))
        sigma14_subset <- crossprod(sqrt(lambda[14]*v14)*t(X14_subset))
        sigma15_subset <- crossprod(sqrt(lambda[15]*v15)*t(X15_subset))
        sigma16_subset <- crossprod(sqrt(lambda[16]*v16)*t(X16_subset))

        diag(sigma2_subset) <- diag(sigma2_subset) + lambda[2]/w_y_subset
        diag(sigma4_subset) <- diag(sigma4_subset) + lambda[4]/w_y_subset
        diag(sigma14_subset) <- diag(sigma14_subset) + lambda[14]/w_y_subset
        diag(sigma15_subset) <- diag(sigma15_subset) + lambda[15]/w_y_subset
        diag(sigma16_subset) <- diag(sigma16_subset) + lambda[16]/w_y_subset

        if (calc_ln_prob_data[2]){sigma2_chol_subset <- chol(sigma2_subset)}
        if (calc_ln_prob_data[4]){sigma4_chol_subset <- chol(sigma4_subset)}
        if (calc_ln_prob_data[14]){sigma14_chol_subset <- chol(sigma14_subset)}
        if (calc_ln_prob_data[15]){sigma15_chol_subset <- chol(sigma15_subset)}
        if (calc_ln_prob_data[16]){sigma16_chol_subset <- chol(sigma16_subset)}

        #compute likelihoods for H1-H16
        if (calc_ln_prob_data[1]){ln_prob_data[i,1] <- ln_prob_data1}
        if (calc_ln_prob_data[3]){ln_prob_data[i,3] <- ln_prob_data3}
        if (calc_ln_prob_data[13]){ln_prob_data[i,13] <- ln_prob_data13}

        if (calc_ln_prob_data[2]){ln_prob_data[i,2] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma2_chol_subset, df = kappa[2])}
        if (calc_ln_prob_data[4]){ln_prob_data[i,4] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma4_chol_subset, df = kappa[4])}
        if (calc_ln_prob_data[5]){ln_prob_data[i,5] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma5_chol_subset, df = kappa[5])}
        if (calc_ln_prob_data[6]){ln_prob_data[i,6] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma6_chol_subset, df = kappa[6])}
        if (calc_ln_prob_data[7]){ln_prob_data[i,7] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma7_chol_subset, df = kappa[7])}
        if (calc_ln_prob_data[8]){ln_prob_data[i,8] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma8_chol_subset, df = kappa[8])}
        if (calc_ln_prob_data[9]){ln_prob_data[i,9] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma9_chol_subset, df = kappa[9])}
        if (calc_ln_prob_data[10]){ln_prob_data[i,10] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma10_chol_subset, df = kappa[10])}
        if (calc_ln_prob_data[11]){ln_prob_data[i,11] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma11_chol_subset, df = kappa[11])}
        if (calc_ln_prob_data[12]){ln_prob_data[i,12] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma12_chol_subset, df = kappa[12])}
        if (calc_ln_prob_data[14]){ln_prob_data[i,14] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma14_chol_subset, df = kappa[14])}
        if (calc_ln_prob_data[15]){ln_prob_data[i,15] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma15_chol_subset, df = kappa[15])}
        if (calc_ln_prob_data[16]){ln_prob_data[i,16] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma16_chol_subset, df = kappa[16])}
      }
    }
  }

  #compute posterior probabilities for all cases
  #cases encoded by presence (1) or absence (0) of 'X->m, m->y, X->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y' and (i) denotes moderation
  #c1:  '0,0,0' / H1 and H5
  #c2:  '0,1,0' / H2 and H5
  #c3:  '1,0,0' / H1 and H6
  #c4:  '1,1,0' / H2 and H6 - complete mediation
  #c5:  '0,0,1' / H3 and H5
  #c6:  '0,1,1' / H4 and H5
  #c7:  '1,0,1' / H3 and H6 - colocalization
  #c8:  '1,1,1' / H4 and H6 - partial mediation
  #c9:  '0,*,0' / H1 and H7
  #c10: '1,*,0' / H1 and H8
  #c11: '0,*,1' / H3 and H7 - Reactive
  #c12: '1,*,1' / H3 and H8
  #c13: '0,i,0' / H5 and H15
  #c14: 'i,0,0' / H1 and H9
  #c15: 'i,1,0' / H2 and H9
  #c16: '1,i,0' / H15 and H6
  #c17: '0,0,i' / H13 and H5
  #c18: '0,1,i' / H14 and H5
  #c19: '0,i,1' / H16 and H5
  #c20: '1,0,i' / H13 and H6
  #c21: 'i,0,1' / H3 and H9
  #c22: '1,1,i' / H14 and H6
  #c23: '1,i,1' / H4 and H9
  #c24: '1,i,1' / H16 and H6
  #c25: '0,*i,0' / H1 and H11
  #c26: 'i,*,0' / H1 and H10
  #c27: '1,*i,0' / H1 and H12
  #c28: '0,*,i' / H13 and H7
  #c29: '0,*i,1' / H3 and H11
  #c30: '1,*,i' / H13 and H8
  #c31: 'i,*,1' / H3 and H10
  #c32: '1,*i,1' / H3 and H12

  preset_odds_index <- return_preset_odds_index_moderation()
  output <- posterior_summary_moderation(ln_prob_data, ln_prior_c, preset_odds_index)
  colnames(output$ln_post_odds) <- colnames(output$ln_prior_odds) <- colnames(output$ln_post_odds)

  #return results
  output$ln_prior_c <- matrix(ln_prior_c, nrow = 1)
  colnames(output$ln_prior_c) <- colnames(output$ln_post_c)

  output$ln_prob_data <- ln_prob_data
  output <- output[c("ln_prob_data", "ln_post_c", "ln_post_odds", "ln_prior_c", "ln_prior_odds", "ln_ml")]

  if (verbose) { print("Done", quote = FALSE) }
  output
}


#' Bayesian model selection for moderated mediation analysis function, version 0
#'
#' This function takes an outcome (y), candidate mediators (M), a candidate moderator (t), and a driver as a design matrix (X) to perform a
#' Bayesian model selection analysis for mediation.
#'
#' @param y Vector or single column matrix of an outcome variable. Single outcome variable expected.
#' Names or rownames must match across M, X, t, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported.
#' Names or rownames must match across y, X, t, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param X Design matrix of the driver. Names or rownames must match across y, M, t, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE.
#' If align_data = FALSE, dimensions and order must match across inputs. One common application is for X to represent genetic information at a QTL,
#' either as founder strain haplotypes or variant genotypes, though X is generalizable to other types of variables.
#' @param t Design matrix of the moderator for the outcome vand mediator. Names or rownames must match across y, M, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE.
#' If align_data = FALSE, dimensions and order must match across inputs.
#' @param Z DEFAULT: NULL. Design matrix of covariates that influence the outcome and mediator variables.
#' Names or rownames must match to those of y, M, X, t, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data=FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param Z_y DEFAULT: NULL. Design matrix of covariates that influence the outcome variable.
#' Names or rownames must match to those of y, M, X, t, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param Z_M DEFAULT: NULL. Design matrix of covariates that influence the mediator variables.
#' Names or rownames must match across y, M, X, t, Z_y, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param w DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis that applies to both
#' y and M. Names must match across y, M, X, t, Z, Z_y, and Z_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where w
#' is a vector of the number of individuals per strain. If no w, w_y, or w_M is given, observations are equally weighted as 1s for y and M.
#' If w is provided, it supercedes w_y and w_M.
#' @param w_y DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of y. Names must match across y, M, X, t, Z, Z_y, Z_M, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_y is a vector of the number of individuals per strain used to
#' measure y. If no w_y (or w) is given, observations are equally weighted as 1s for y.
#' @param w_M DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of M. Names must match across y, M, X, t, Z, Z_y, Z_M, and w_y (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_M is a vector of the number of individuals per strain use to
#' measure M. If no w_M (or w) is given, observations are equally weighted as 1s for M.
#' @param tau_sq_mu DEFAULT: 1000. Variance component for the intercept. The DEFAULT represents a diffuse prior, analagous to
#' a fixed effect term.
#' @param tau_sq_Z DEFAULT: 1000. Variance component for the covariates encoded in Z. The DEFAULT represents a diffuse prior, analagous
#' to fixed effect terms.
#' @param phi_sq DEFAULT: c(1, 1, 1). Each element of (a, b, c) represents one of the relationships being evaluated for mediation,
#' specifically the ratio of signal to noise. a is the effect of X on M, b is the effect of M on y, and c is the effect of X on y.
#' The DEFAULT represents relationships that explain 50% of the variation in the outcome variable.
#' @param ln_prior_c DEFAULT: "complete". The prior log case probabilities. See model_info() for description of likelihoods and their
#' combinations into cases. Simplified pre-set options are available, including "complete", "partial", and "reactive".
#' @param align_data DEFAULT: TRUE. If TRUE, expect vector and matrix inputes to have names and rownames, respectively. The overlapping data
#' will then be aligned, allowing the user to not have to reduce data to overlapping samples and order them.
#' @export
#' @examples bmediatR()
bmediatR_v0_moderation <- function(y, M, X, t,
                                Z = NULL, Z_y = NULL, Z_M = NULL,
                                w = NULL, w_y = NULL, w_M = NULL,
                                kappa = 0.001, lambda = 0.001,
                                tau_sq_mu = 1000,
                                tau_sq_Z = 1000,
                                tau_sq_t = 1000, # Prior for the marginal effect of t
                                phi_sq_X = c(NA, NA, 1, 1, NA, 1, NA, 1,
                                             1, 1, NA, 1, 1, 1, NA, 1),
                                phi_sq_m = c(NA, 1, NA, 1, NA, NA, NA, NA,
                                             NA, NA, NA, NA, NA, 1, 1, 1),
                                phi_sq_y = c(NA, NA, NA, NA, NA, NA, 1, 1,
                                             NA, 1, 1, 1, NA, NA, NA, NA),
                                phi_sq_int = c(rep(NA, 8), # effect size prior for interaction term
                                               1, 1, 1, 1, 1, 1, 1, 1),
                                ln_prior_c = "complete",
                                options_X = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                                options_t = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                                align_data = TRUE,
                                verbose = TRUE) {

  #presets for ln_prior_c;
  ln_prior_c <- return_ln_prior_c_from_presets_moderation(ln_prior_c = ln_prior_c)

  #ensure ln_prior_c sum to 1 on probability scale
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
  }

  #optionally align data
  processed_data <- process_data_moderation(y = y, M = M, X = X, t= t,
                                            Z_y = Z_y, Z_M = Z_M,
                                            w_y = w_y, w_M = w_M,
                                            align_data = align_data,
                                            verbose = verbose)
  y <- processed_data$y
  M <- processed_data$M
  X <- processed_data$X
  t <- processed_data$t
  Z_y <- processed_data$Z_y
  Z_M <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_M <- processed_data$w_M

  #dimension of y
  n <- length(y)

  #dimension of Z's
  p_y <- ncol(Z_y)
  p_M <- ncol(Z_M)

  #scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p_y > 0) { Z_y <- apply(Z_y, 2, scale) }
  if (p_M > 0) { Z_M <- apply(Z_M, 2, scale) }

  #optionally use sum-to-zero contrast for X
  #recommended when X is a matrix of factors, with a column for every factor level
  if (options_X$sum_to_zero == TRUE) {
    C <- sumtozero_contrast(ncol(X))
    X <- X%*%C
  }

  #optionally use sum-to-zero contrast for t
  #recommended when t is a matrix of factors, with a column for every factor level
  if (options_t$sum_to_zero == TRUE) {
    C <- sumtozero_contrast(ncol(t))
    t <- t%*%C
  }

  #optionally center and scale X
  X <- apply(X, 2, scale, center = options_X$center, scale = options_X$scale)

  #optionally center and scale t
  t <- apply(t, 2, scale, center = options_t$center, scale = options_t$scale)

  #dimension of X and t
  d <- ncol(X)
  k <- ncol(t)

  # Interaction between X and moderator
  tX <- model.matrix(~ 0 + X:t)
  # Dimension of interaction
  ix <- ncol(tX)

  # Interaction between y and moderator
  ty <- model.matrix(~ 0 + y:t)
  # Dimension of interaction
  iy <- ncol(ty)

  #column design matrix for mu
  ones <- matrix(1, n)

  #begin Bayesian calculations
  if (verbose) { print("Initializing", quote = FALSE) }

  #reformat priors
  kappa = rep(kappa, 16)
  lambda = rep(lambda, 16)
  tau_sq_mu = rep(tau_sq_mu, 16)
  tau_sq_Z = rep(tau_sq_Z, 16)
  tau_sq_t = rep(tau_sq_t, 16)

  #identify likelihoods that are not supported by the prior
  #will not compute cholesky or likelihood for these
  calc_ln_prob_data <- rep(NA, 16)
  calc_ln_prob_data[1] <- any(!is.infinite(ln_prior_c[c(1,3,9,10)]))
  calc_ln_prob_data[2] <- any(!is.infinite(ln_prior_c[c(2,4)]))
  calc_ln_prob_data[3] <- any(!is.infinite(ln_prior_c[c(5,7,11,12)]))
  calc_ln_prob_data[4] <- any(!is.infinite(ln_prior_c[c(6,8)]))
  calc_ln_prob_data[5] <- any(!is.infinite(ln_prior_c[c(1,2,5,6)]))
  calc_ln_prob_data[6] <- any(!is.infinite(ln_prior_c[c(3,4,7,8)]))
  calc_ln_prob_data[7] <- any(!is.infinite(ln_prior_c[c(9,11)]))
  calc_ln_prob_data[8] <- any(!is.infinite(ln_prior_c[c(10,12)]))
  calc_ln_prob_data[9] <- any(!is.infinite(ln_prior_c[c(14,15,21,23)]))
  calc_ln_prob_data[10] <- any(!is.infinite(ln_prior_c[c(26, 31)]))
  calc_ln_prob_data[11] <- any(!is.infinite(ln_prior_c[c(25,29)]))
  calc_ln_prob_data[12] <- any(!is.infinite(ln_prior_c[c(27,32)]))
  calc_ln_prob_data[13] <- any(!is.infinite(ln_prior_c[c(28,30)]))
  calc_ln_prob_data[14] <- any(!is.infinite(ln_prior_c[c(18,22)]))
  calc_ln_prob_data[15] <- any(!is.infinite(ln_prior_c[c(13,16)]))
  calc_ln_prob_data[16] <- any(!is.infinite(ln_prior_c[c(19,24)]))

  #likelihood models for all hypothesis
  #hypotheses encoded by presence (1) or absence (0) of 'X->m, y->m, X->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y', i denotes moderation
  #H1: '-,0,0' / y does not depend on X or m
  #H2: '-,1,0' / y depends on m but not X
  #H3: '-,0,1' / y depends on X but not m
  #H4: '-,1,1' / y depends on X and m
  #H5: '0,-,-' / m does not depend on X
  #H6: '1,-,-' / m depends on X
  #H7: '0,*,-' / m depends on y but not X
  #H8: '1,*,-' / m depends on X and y
  #H9: 'i,-,-' / m depends on X and is moderated by t
  #H10: 'i,*,-' / m depends on X and y and relation between X and M is moderated by t
  #H11: '0,*i,-' / m depends on y but not X and is moderated by t
  #H12: '1,*i,-' / m depends on X and y and y -> M is moderated by t
  #H13: '-,0,i' / y depends on X but not m and is moderated by t
  #H14: '-,1,i' / y depends on X and m and X->y is moderated by t
  #H15: '-,i,0' / y depends on m but not X and is moderated by t
  #H16: '-,i,1' / y depends on X and m and m -> y is moderated by t
  #all include covariates Z and t

  #design matrices for H1,H3,H5-H13 complete cases (do not depend on m)
  X1 <- cbind(ones, Z_y, t)
  X3 <- cbind(ones, X, Z_y, t)
  X5 <- cbind(ones, Z_M, t)
  X6 <- cbind(ones, X, Z_M, t)
  X7 <- cbind(X5, y)
  X8 <- cbind(X6, y)
  X9 <- cbind(X6, tX)
  X10 <- cbind(X8, tX)
  X11 <- cbind(X7, ty)
  X12 <- cbind(X8, ty)
  X13 <- cbind(X3, tX)

  #check if all scale hyperparameters are identical for H1 and H5
  #implies sigma1 and sigma5 identical, used to reduce computations
  sigma5_equal_sigma1 <- all(lambda[1]==lambda[5],
                             tau_sq_mu[1] == tau_sq_mu[5],
                             tau_sq_Z[1] == tau_sq_Z[5],
                             identical(Z_y, Z_M),
                             identical(w_y, w_M))

  #check if all scale hyperparameters are identical for H3 and H6
  #implies sigma3 and sigma6 identical, used to reduce computations
  sigma6_equal_sigma3 <- all(lambda[3]==lambda[6],
                             tau_sq_mu[3] == tau_sq_mu[6],
                             tau_sq_Z[3] == tau_sq_Z[6],
                             identical(Z_y, Z_M),
                             identical(w_y, w_M))

  #check if all scale hyperparameters are identical for H9 and H13
  #implies sigma3 and sigma6 identical, used to reduce computations
  sigma13_equal_sigma9 <- all(lambda[9]==lambda[13],
                              tau_sq_mu[9] == tau_sq_mu[13],
                              tau_sq_Z[9] == tau_sq_Z[13],
                              identical(Z_y, Z_M),
                              identical(w_y, w_M),
                              tau_sq_t[9] == tau_sq_t[13])


  #prior variance matrices (diagonal) for H1-H16
  v1 <- c(tau_sq_mu[1], rep(tau_sq_Z[1], p_y), rep(tau_sq_t[1], k))
  v2 <- c(tau_sq_mu[2], rep(tau_sq_Z[2], p_y), rep(tau_sq_t[2], k), phi_sq_m[2])
  v3 <- c(tau_sq_mu[3], rep(phi_sq_X[3], d), rep(tau_sq_Z[3], p_y), rep(tau_sq_t[3], k))
  v4 <- c(tau_sq_mu[4], rep(phi_sq_X[4], d), rep(tau_sq_Z[4], p_y), rep(tau_sq_t[4], k), phi_sq_m[4])
  v7 <- c(tau_sq_mu[7], rep(tau_sq_Z[7], p_M), rep(tau_sq_t[7], k), phi_sq_y[7])
  v8 <- c(tau_sq_mu[8], rep(phi_sq_X[8], d), rep(tau_sq_Z[8], p_M), rep(tau_sq_t[8], k), phi_sq_y[8])
  v9 <- c(tau_sq_mu[9], rep(phi_sq_X[9], d), rep(tau_sq_Z[9], p_M), rep(tau_sq_t[9], k), rep(phi_sq_int[9], ix))
  v10 <- c(tau_sq_mu[10], rep(phi_sq_X[10], d), rep(tau_sq_Z[10], p_M), rep(tau_sq_t[10], k), rep(phi_sq_int[10], ix))
  v11 <- c(tau_sq_mu[11], rep(tau_sq_Z[11], p_M), rep(tau_sq_t[11], k), phi_sq_y[11], rep(phi_sq_int[11], iy))
  v12 <- c(tau_sq_mu[12], rep(phi_sq_X[12], d), rep(tau_sq_Z[12], p_M), rep(tau_sq_t[12], k), rep(phi_sq_int[12], iy))
  v14 <- c(tau_sq_mu[14], rep(phi_sq_X[14], d), rep(tau_sq_Z[14], p_y), rep(tau_sq_t[14], k), phi_sq_m[14], rep(phi_sq_int[14], ix))


  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]){
    v5 <- c(tau_sq_mu[5], rep(tau_sq_Z[5], p_M), rep(tau_sq_t[5], k))
  }

  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]){
    v6 <- c(tau_sq_mu[6], rep(phi_sq_X[6], d), rep(tau_sq_Z[6], p_M), rep(tau_sq_t[5], k))
  }

  if (!sigma13_equal_sigma9 | !calc_ln_prob_data[9]){
    v13 <- c(tau_sq_mu[13], rep(phi_sq_X[13], d), rep(tau_sq_Z[13], p_y), rep(tau_sq_t[13], k), rep(phi_sq_int[13], ix))
  }

  #scale matrices for H1,H3,H5-H8 complete cases, H9-H13 (do not depend on m)
  sigma1 <- crossprod(sqrt(lambda[1]*v1)*t(X1))
  sigma3 <- crossprod(sqrt(lambda[3]*v3)*t(X3))
  sigma7 <- crossprod(sqrt(lambda[7]*v7)*t(X7))
  sigma8 <- crossprod(sqrt(lambda[8]*v8)*t(X8))
  sigma9 <- crossprod(sqrt(lambda[9]*v9)*t(X9))
  sigma10 <- crossprod(sqrt(lambda[10]*v10)*t(X10))
  sigma11 <- crossprod(sqrt(lambda[11]*v11)*t(X11))
  sigma12 <- crossprod(sqrt(lambda[12]*v12)*t(X12))

  diag(sigma1) <- diag(sigma1) + lambda[1]/w_y
  diag(sigma3) <- diag(sigma3) + lambda[3]/w_y
  diag(sigma7) <- diag(sigma7) + lambda[7]/w_M
  diag(sigma8) <- diag(sigma8) + lambda[8]/w_M
  diag(sigma9) <- diag(sigma9) + lambda[9]/w_M
  diag(sigma10) <- diag(sigma10) + lambda[10]/w_M
  diag(sigma11) <- diag(sigma11) + lambda[11]/w_M
  diag(sigma12) <- diag(sigma12) + lambda[12]/w_M

  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]){
    sigma5 <- crossprod(sqrt(lambda[5]*v5)*t(X5))
    diag(sigma5) <- diag(sigma5) + lambda[5]/w_M
  }

  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]){
    sigma6 <- crossprod(sqrt(lambda[6]*v6)*t(X6))
    diag(sigma6) <- diag(sigma6) + lambda[6]/w_M
  }

  if (!sigma13_equal_sigma9 | !calc_ln_prob_data[9]){
    sigma13 <- crossprod(sqrt(lambda[13]*v13)*t(X13))
    diag(sigma13) <- diag(sigma13) + lambda[13]/w_y
  }

  #object to store likelihoods
  ln_prob_data=matrix(-Inf, ncol(M), 16)
  rownames(ln_prob_data) <- colnames(M)
  colnames(ln_prob_data) <- c("-,0,0",
                              "-,1,0",
                              "-,0,1",
                              "-,1,1",
                              "0,-,-",
                              "1,-,-",
                              "0,*,-",
                              "1,*,-",
                              "i,-,-",
                              "i,*,-",
                              "0,*i,-",
                              "1,*i,-",
                              "-,0,i",
                              "-,1,i",
                              "-,i,0",
                              "-,i,1")

  #identify batches of M that have the same pattern of missing values
  missing_m <- bmediatR:::batch_cols(M)

  #iterate over batches of M with same pattern of missing values
  if (verbose) { print("Iterating", quote = FALSE) }
  counter <- 0

  for (b in 1:length(missing_m)){
    #subset to non-missing observations
    index <- rep(T, length(y))
    index[missing_m[[b]]$omit] <- FALSE

    if (any(index)){
      y_subset <- y[index]
      w_y_subset <- w_y[index]
      w_M_subset <- w_M[index]

      #cholesky matrices for H1,H3,H5-H8 non-missing observations (do not depend on m)
      if (calc_ln_prob_data[1]){sigma1_chol_subset <- chol(sigma1[index,index])}
      if (calc_ln_prob_data[3]){sigma3_chol_subset <- chol(sigma3[index,index])}
      if (calc_ln_prob_data[7]){sigma7_chol_subset <- chol(sigma7[index,index])}
      if (calc_ln_prob_data[8]){sigma8_chol_subset <- chol(sigma8[index,index])}
      if (calc_ln_prob_data[9]){sigma9_chol_subset <- chol(sigma9[index,index])}
      if (calc_ln_prob_data[10]){sigma10_chol_subset <- chol(sigma10[index,index])}
      if (calc_ln_prob_data[11]){sigma11_chol_subset <- chol(sigma11[index,index])}
      if (calc_ln_prob_data[12]){sigma12_chol_subset <- chol(sigma12[index,index])}

      if (sigma5_equal_sigma1 & calc_ln_prob_data[1]) {
        sigma5_chol_subset <- sigma1_chol_subset
      } else if (calc_ln_prob_data[5]) {
        sigma5_chol_subset <- chol(sigma5[index,index])
      }

      if (sigma6_equal_sigma3 & calc_ln_prob_data[3]) {
        sigma6_chol_subset <- sigma3_chol_subset
      } else if (calc_ln_prob_data[6]) {
        sigma6_chol_subset <- chol(sigma6[index,index])
      }

      if (sigma13_equal_sigma9 & calc_ln_prob_data[9]){
        sigma13_chol_subset <- sigma9_chol_subset
      } else if (calc_ln_prob_data[13]){
        sigma13_chol_subset <- chol(sigma13[index,index])
      }

      #compute H1 and H3 outside of the mediator loop (invariant)
      if (calc_ln_prob_data[1]) { ln_prob_data1 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma1_chol_subset, df = kappa[1]) }
      if (calc_ln_prob_data[3]) { ln_prob_data3 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma3_chol_subset, df = kappa[3]) }
      if (calc_ln_prob_data[13]){ln_prob_data13 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma13_chol_subset, df = kappa[13])}

      #iterate over mediators
      for (i in missing_m[[b]]$cols) {
        counter <- counter + 1
        if (counter%%1000==0 & verbose) { print(paste(counter, "of", ncol(M)), quote=F) }

        #set current mediator non-missing observations
        m_subset <- M[index,i]
        t_subset <- t[index,,drop = F]

        # Calculate design matrix for interaction of M and t
        tM <- model.matrix(~ 0 + m_subset:t_subset)
        im <- ncol(tM)

        #design matrix for H2 and H4 non-missing observations
        X2_subset <- cbind(X1[index,,drop = FALSE], m_subset)
        X4_subset <- cbind(X3[index,,drop = FALSE], m_subset)
        X14_subset <- cbind(X4_subset, tX[index,,drop=FALSE])
        X15_subset <- cbind(X2_subset, tM)
        X16_subset <- cbind(X4_subset, tM)

        v15 <- c(tau_sq_mu[15], rep(tau_sq_Z[15], p_y), rep(tau_sq_t[15], k), phi_sq_m[15], rep(phi_sq_int[15], im))
        v16 <- c(tau_sq_mu[16], rep(phi_sq_X[16], d), rep(tau_sq_Z[16], p_y), rep(tau_sq_t[16], k), phi_sq_m[16], rep(phi_sq_int[15], im))

        #scale and cholesky matrices for H2 and H4 non-missing observations
        sigma2_subset <- crossprod(sqrt(lambda[2]*v2)*t(X2_subset))
        sigma4_subset <- crossprod(sqrt(lambda[4]*v4)*t(X4_subset))
        sigma14_subset <- crossprod(sqrt(lambda[14]*v14)*t(X14_subset))
        sigma15_subset <- crossprod(sqrt(lambda[15]*v15)*t(X15_subset))
        sigma16_subset <- crossprod(sqrt(lambda[16]*v16)*t(X16_subset))

        diag(sigma2_subset) <- diag(sigma2_subset) + lambda[2]/w_y_subset
        diag(sigma4_subset) <- diag(sigma4_subset) + lambda[4]/w_y_subset
        diag(sigma14_subset) <- diag(sigma14_subset) + lambda[14]/w_y_subset
        diag(sigma15_subset) <- diag(sigma15_subset) + lambda[15]/w_y_subset
        diag(sigma16_subset) <- diag(sigma16_subset) + lambda[16]/w_y_subset

        if (calc_ln_prob_data[2]){sigma2_chol_subset <- chol(sigma2_subset)}
        if (calc_ln_prob_data[4]){sigma4_chol_subset <- chol(sigma4_subset)}
        if (calc_ln_prob_data[14]){sigma14_chol_subset <- chol(sigma14_subset)}
        if (calc_ln_prob_data[15]){sigma15_chol_subset <- chol(sigma15_subset)}
        if (calc_ln_prob_data[16]){sigma16_chol_subset <- chol(sigma16_subset)}

        #compute likelihoods for H1-H16
        if (calc_ln_prob_data[1]){ln_prob_data[i,1] <- ln_prob_data1}
        if (calc_ln_prob_data[3]){ln_prob_data[i,3] <- ln_prob_data3}
        if (calc_ln_prob_data[13]){ln_prob_data[i,13] <- ln_prob_data13}

        if (calc_ln_prob_data[2]){ln_prob_data[i,2] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma2_chol_subset, df = kappa[2])}
        if (calc_ln_prob_data[4]){ln_prob_data[i,4] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma4_chol_subset, df = kappa[4])}
        if (calc_ln_prob_data[5]){ln_prob_data[i,5] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma5_chol_subset, df = kappa[5])}
        if (calc_ln_prob_data[6]){ln_prob_data[i,6] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma6_chol_subset, df = kappa[6])}
        if (calc_ln_prob_data[7]){ln_prob_data[i,7] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma7_chol_subset, df = kappa[7])}
        if (calc_ln_prob_data[8]){ln_prob_data[i,8] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma8_chol_subset, df = kappa[8])}
        if (calc_ln_prob_data[9]){ln_prob_data[i,9] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma9_chol_subset, df = kappa[9])}
        if (calc_ln_prob_data[10]){ln_prob_data[i,10] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma10_chol_subset, df = kappa[10])}
        if (calc_ln_prob_data[11]){ln_prob_data[i,11] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma11_chol_subset, df = kappa[11])}
        if (calc_ln_prob_data[12]){ln_prob_data[i,12] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma12_chol_subset, df = kappa[12])}
        if (calc_ln_prob_data[14]){ln_prob_data[i,14] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma14_chol_subset, df = kappa[14])}
        if (calc_ln_prob_data[15]){ln_prob_data[i,15] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma15_chol_subset, df = kappa[15])}
        if (calc_ln_prob_data[16]){ln_prob_data[i,16] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma16_chol_subset, df = kappa[16])}
      }
    }
  }

  #compute posterior probabilities for all cases
  #cases encoded by presence (1) or absence (0) of 'X->m, m->y, X->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y' and (i) denotes moderation
  #c1:  '0,0,0' / H1 and H5
  #c2:  '0,1,0' / H2 and H5
  #c3:  '1,0,0' / H1 and H6
  #c4:  '1,1,0' / H2 and H6 - complete mediation
  #c5:  '0,0,1' / H3 and H5
  #c6:  '0,1,1' / H4 and H5
  #c7:  '1,0,1' / H3 and H6 - colocalization
  #c8:  '1,1,1' / H4 and H6 - partial mediation
  #c9:  '0,*,0' / H1 and H7
  #c10: '1,*,0' / H1 and H8
  #c11: '0,*,1' / H3 and H7 - Reactive
  #c12: '1,*,1' / H3 and H8
  #c13: '0,i,0' / H5 and H15
  #c14: 'i,0,0' / H1 and H9
  #c15: 'i,1,0' / H2 and H9
  #c16: '1,i,0' / H15 and H6
  #c17: '0,0,i' / H13 and H5
  #c18: '0,1,i' / H14 and H5
  #c19: '0,i,1' / H16 and H5
  #c20: '1,0,i' / H13 and H6
  #c21: 'i,0,1' / H3 and H9
  #c22: '1,1,i' / H14 and H6
  #c23: '1,i,1' / H4 and H9
  #c24: '1,i,1' / H16 and H6
  #c25: '0,*i,0' / H1 and H11
  #c26: 'i,*,0' / H1 and H10
  #c27: '1,*i,0' / H1 and H12
  #c28: '0,*,i' / H13 and H7
  #c29: '0,*i,1' / H3 and H11
  #c30: '1,*,i' / H13 and H8
  #c31: 'i,*,1' / H3 and H10
  #c32: '1,*i,1' / H3 and H12

  preset_odds_index <- return_preset_odds_index_moderation()
  output <- posterior_summary_moderation(ln_prob_data, ln_prior_c, preset_odds_index)
  colnames(output$ln_post_odds) <- colnames(output$ln_prior_odds) <- colnames(output$ln_post_odds)

  #return results
  output$ln_prior_c <- matrix(ln_prior_c, nrow = 1)
  colnames(output$ln_prior_c) <- colnames(output$ln_post_c)

  output$ln_prob_data <- ln_prob_data
  output <- output[c("ln_prob_data", "ln_post_c", "ln_post_odds", "ln_prior_c", "ln_prior_odds", "ln_ml")]

  if (verbose) { print("Done", quote = FALSE) }
  output
}


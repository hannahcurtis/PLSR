# Algorithm to run PLSR model and select variables based on VIP score/AIC

library(vip)
library(pls)

# function to automatically select the number of variables to use via vip
run_model <- function(X, Y, df, ncomp=5) {
  set.seed(0)
  
  # fit a PLS model
  formula <- as.formula(paste(paste(Y,collapse="+"),paste(X,collapse="+"),sep="~"))
  return(plsr(formula, data=df, ncomp=ncomp, scale=TRUE, validation="CV", segments=5))
}

compute_vi <- function(pls_model, ncomp=5) {
  # calculate the vi scores
  vi_df <- data.frame(vi(pls_model, ncomp=ncomp, scale=FALSE))
  
  # Arrange data frame from lowest to highest VIP scores
  vi_df <- arrange(vi_df, Importance)
  
  return(vi_df)
}

compute_aic <- function(model, ncomp) {
  n <- dim(model$residuals)[1]
  return(n * log(sum(model$residuals[,,ncomp]^2)/n) + 2 * ncomp)
}

get_ncomp_for_min_aic <- function(model) {
  min_aic <- 1000000
  min_aic_ncomp <- -1
  
  # for the number of components that we're selecting
  for (k in 1:dim(model$residuals)[3]) {
    aic = compute_aic(model, i) 
    
    if (aic < min_aic) {
      min_aic <- aic
      min_aic_ncomp <- i
    }
  }
  return(min_aic_ncomp)
}

get_ncomp_for_max_qscore <- function(q_score) {
  # find the max of q_score "w/o overfitting" by forcing it to always increase by 10%
  max_qscore <- q_score$val[2]
  max_ncomp <- 1
  for (i in 3:length(q_score$comps)) {
    if (q_score$val[i] > max_qscore && q_score$val[i] > max_qscore * 1.1) {
      max_qscore <- q_score$val[i]
      max_ncomp <- i-1
    } else {
      break
    }
  }
  return(max_ncomp)
}

get_blank_df <- function() {
  return(data.frame(max_q_score=numeric(), min_aic=numeric(), nvar=numeric(), ncomp=numeric(), X=character(), Y=character(), last_removed=character(), least_vi=character(), least_vi_val=numeric()))
}

eval_vars <- function(X, Y, df, ncomp, to_remove=NULL, method="q-score") {
  if (!is.null(to_remove)) {
    X <- X[X!=to_remove]
  }
  # fit a PLS model
  pls_model <- run_model(X, Y, df, ncomp)
  q_score <- pls::R2(pls_model)
  if (method == "aic") {
    best_ncomp <- get_ncomp_for_min_aic(pls_model)
  } else if (method == "q-score") {
    best_ncomp <- get_ncomp_for_max_qscore(q_score)
  } 
  min_aic <- compute_aic(pls_model, best_ncomp)
  max_q_score <- q_score$val[best_ncomp+1]
  
  vi <- compute_vi(pls_model, best_ncomp)
  # clean up
  if (is.null(to_remove)) {
   to_remove <- ""
  }
  return(data.frame(
    max_q_score=max_q_score, min_aic=min_aic, nvar=length(X), ncomp=best_ncomp, X=paste(X, collapse="+"), Y=paste(Y, collapse="+"), last_removed=to_remove, least_vi=vi$Variable[1], least_vi_val=vi$Importance[1]
  ))
}

iter_fun <- function(X, Y, df, ncomp, inner_func, comp_func, method) {
  sub_model_results <- get_blank_df()
  # try removing each variable individually and compute the q-score
  for (to_remove in X) {
    res_df <- eval_vars(X, Y, df, min(ncomp, length(X)-2), to_remove, method)
    sub_model_results <- rbind(sub_model_results, res_df)
  }
  # select the row with the highest q-score
  best_row <- sub_model_results[comp_func(inner_func(sub_model_results)),]
  
  return(best_row)
}

q_score_obj <- function(X, Y, df, ncomp) {
  return(iter_fun(X, Y, df, ncomp, function(x) x$max_q_score, which.max, "q-score"))
}

vi_objective <- function(X,Y,df,ncomp) {
  res_df <- eval_vars(X, Y, df, min(ncomp, length(X)-2), NULL, "q-score")
  res_df <- eval_vars(X, Y, df, min(ncomp, length(X)-2), res_df$least_vi, "q-score")
  return(res_df)
}

vi_aic_objective <- function(X,Y,df,ncomp) {
  res_df <- eval_vars(X, Y, df, min(ncomp, length(X)-2), NULL, "aic")
  res_df <- eval_vars(X, Y, df, min(ncomp, length(X)-2), res_df$least_vi, "aic")
  return(res_df)  
}

aic_objective <- function(X,Y,df,ncomp) {
  return(iter_fun(X, Y, df, ncomp, function(x) x$aic, which.min, "aic"))
}

get_best_vars <- function(X, Y, df, ncomp=5, objective="q-score") {
  if (objective == "q-score") { obj_fun <- q_score_obj }
  else if (objective == "vi") { obj_fun <- vi_objective }
  else if (objective == "vi_aic") { obj_fun <- vi_aic_objective }
  else if (objective == "aic") { obj_fun <- aic_objective }
  else { stop("objective must be 'q-score', 'aic', 'vi', or 'vi_aic") }
  
  # define a dataframe with columns X, Y, qscore, ncomp
  model_results <- get_blank_df()
  model_results <- rbind(model_results, eval_vars(X, Y, df, ncomp))
  while (length(X) > 3)
  {
    new_row <- obj_fun(X,Y,df,ncomp)
    if (nrow(new_row) != 1) {
      print(c("binding new row of length", length(new_row)))
      print(new_row)
      stop()
    }
    model_results <- rbind(model_results, new_row)
    X <- X[X!=new_row$last_removed]
  }
  return(model_results)
}

get_best_vars_both <- function(X, Y, df, ncomp=5)
{
  
  q_score_results <- get_best_vars(X, Y, df, ncomp, "q-score")
  vi_results <- get_best_vars(X, Y, df, ncomp, "vi")
  # get the row w/ the highest q-score
  best_q_score <- q_score_results[which.max(q_score_results$max_q_score),]
  # get the row w/ the highest vi score
  best_vi <- vi_results[which.max(vi_results$max_q_score),]
  
  # create the union of the variables
  q_score_x <- strsplit(best_q_score$X, "\\+")[[1]]
  vi_x <- strsplit(best_vi$X, "\\+")[[1]]
  best_vars <- unique(c(q_score_x, vi_x))
  # get the best model using the union of the variables
  best_model_vi_reduce <- get_best_vars(best_vars, Y, df, ncomp, "vi")
  best_both <- best_model_vi_reduce[which.max(best_model_vi_reduce$max_q_score),]
  
  # print useful outputs
  print(c("max q-score:", best_q_score$max_q_score, best_q_score$nvar))
  print(c("q-score X vars only:", paste(setdiff(q_score_x, vi_x), collapse="+")))
  print(c("max vi q2:", best_vi$max_q_score, best_vi$nvar))
  print(c("vi X vars only:", paste(setdiff(vi_x, q_score_x), collapse="+")))
  print(c("max both: ", best_both$max_q_score, best_both$nvar))
  print(c("max vi combined x vars:", paste(best_model_vi_reduce$X, collapse="+")))
  
  print("additively adding X vars from Q2 into results of VI...")
  additive_model_results <- get_blank_df()
  while (length(q_score_x) > 0) {
    sub_model_results <- get_blank_df()
    for (to_add in q_score_x) {
      # add variable to the vector X
      X_temp <- c(vi_x, to_add)
      res_df <- eval_vars(X_temp, Y, df, min(ncomp, length(X_temp)-2))
      res_df$last_removed = to_add
      sub_model_results <- rbind(sub_model_results, res_df)
    }
    best_row <- sub_model_results[which.max(sub_model_results$max_q_score),]
    additive_model_results <- rbind(additive_model_results, best_row)
    # remove the var we added from q_score_x
    q_score_x <- q_score_x[q_score_x != best_row$last_removed]
    vi_x <- c(vi_x, to_add)
  }
  # select the row with the highest q-score
  best_row <- additive_model_results[which.max(additive_model_results$max_q_score),]
  print(c("max q-score from additive model:", best_row$max_q_score))
  print(c("additive model vars unique from vi model:", paste(setdiff(strsplit(best_row$X, "\\+")[[1]], vi_x), collapse="+")))
  
  return(list(q_score_results, vi_results, best_model_vi_reduce, additive_model_results))
}


get_best_vars_generic <- function(X, Y, df, ncomp=5, reduce_model="vi", overall_obj="q-score")
{
  order <- c(overall_obj, reduce_model)
  
  
  print(c("Getting best vars for ", order[1]))
  results_1 <- get_best_vars(X, Y, df, ncomp, order[1])
  print(c("Getting best vars for ", order[2]))
  results_2 <- get_best_vars(X, Y, df, ncomp, order[2])
  print(c("Getting best vars for both using", overall_obj))
  if (overall_obj == "q-score") {
    best_1 <- results_1[which.max(results_1$max_q_score),]
    best_2 <- results_2[which.max(results_2$max_q_score),]
  } else if (overall_obj == "aic") {
    best_1 <- results_1[which.min(results_1$aic),]
    best_2 <- results_2[which.min(results_2$aic),]
  } else {
    stop("overall_obj must be 'q-score' or 'aic'")
  }
  
  # create the union of the variables
  print("taking union of variables")
  x_1 <- strsplit(best_1$X, "\\+")[[1]]
  x_2 <- strsplit(best_2$X, "\\+")[[1]]
  best_vars <- unique(c(x_1, x_2))
  
  # get the best model using the union of the variables
  print(c("reducing union using ", reduce_model))
  best_model_reduce <- get_best_vars(best_vars, Y, df, ncomp, reduce_model)
  if (overall_obj == "q-score") {
    best_both <- best_model_reduce[which.max(best_model_reduce$max_q_score),]
  } else if (overall_obj == "aic") {
    best_both <- best_model_reduce[which.min(best_model_reduce$aic),]
  } else {
    stop("overall_obj must be 'q-score' or 'aic'")
  }

  # print useful outputs
  print(c(order[1], "(q_score, aic, nvar):", best_1$max_q_score, best_1$min_aic, best_1$nvar))
  print(c(order[1], "X vars only:", paste(setdiff(x_1, x_2), collapse="+")))
  print(c(order[2], "(q_score, aic, nvar):", best_2$max_q_score, best_2$min_aic, best_2$nvar))
  print(c(order[2], "X vars only:", paste(setdiff(x_2, x_1), collapse="+")))
  print(c("max both (q_score, aic, nvar):", best_1$max_q_score, best_1$min_aic, best_1$nvar))
  print(c("max", reduce_model, "combined x vars:", paste(best_model_reduce$X, collapse="+")))
  
  print(c("additively adding X vars from ", order[1], " into results of ", reduce_model, "..."))
  additive_model_results <- get_blank_df()
  while (length(x_1) > 0) {
    sub_model_results <- get_blank_df()
    for (to_add in x_1) {
      # add variable to the vector X
      X_temp <- c(x_2, to_add)
      res_df <- eval_vars(X_temp, Y, df, min(ncomp, length(X_temp)-2), overall_obj)
      res_df$last_removed = to_add
      sub_model_results <- rbind(sub_model_results, res_df)
    }
    if (overall_obj == "q-score") {
      best_row <- sub_model_results[which.max(sub_model_results$max_q_score),]
    } else if (overall_obj == "aic") {
      best_row <- sub_model_results[which.min(sub_model_results$aic),]
    } else {
      stop("overall_obj must be 'q-score' or 'aic'")
    }
    additive_model_results <- rbind(additive_model_results, best_row)
    # remove the var we added from q_score_x
    x_1 <- x_1[x_1 != best_row$last_removed]
    x_2 <- c(x_2, to_add)
  }
  # select the row with the highest q-score
  if (overall_obj == "q-score") {
    best_row <- additive_model_results[which.max(additive_model_results$max_q_score),]
  } else if (overall_obj == "aic") {
    best_row <- additive_model_results[which.min(additive_model_results$aic),]
  } else {
    stop("overall_obj must be 'q-score' or 'aic'")
  }
  print(c("max q-score from additive model:", best_row$max_q_score))
  print(c("max aic from additive model:", best_row$aic))
  print(c("additive model vars unique from",reduce_model,"model:", paste(setdiff(strsplit(best_row$X, "\\+")[[1]], x_2), collapse="+")))
  
  return(list(results_1, results_2, best_model_reduce, additive_model_results))
}


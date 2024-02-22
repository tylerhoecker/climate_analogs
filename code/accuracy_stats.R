#' @description Calculate a set of accuracy metrics for AIMs
#' @param predicted 
#' @param expected 
#' @param fact_levels 
#' @return A table of analog and project environmental attributes

accuracy_stats <- function(aim, observed, fact_levels){

  predicted <- factor(
    as.character(aim),
    levels = as.character(lf_levels)
  )

  expected <- factor(
    as.character(observed),
    levels = as.character(lf_levels)
  )
  
  summary_stats <- caret::confusionMatrix(predicted, expected)[["byClass"]]|>
  as.data.frame() |>
  filter(Prevalence > 0, !is.na(F1)) |>
  summarise(across(
    c(Sensitivity,Specificity,Precision,Recall,F1,`Balanced Accuracy`),
    mean
  ))

  return(by_class_stats)
}

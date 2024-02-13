recordAttrition <- function(table, reason){
  n <- attr(table, "cohort_attrition") |> 
    filter(row_number() == max(row_number())) |> 
    pull(N)
  
  attr(table, "cohort_attrition") <- attr(table, "cohort_attrition") |> 
    union_all(
      tibble(
        "N"      = table |> nrow(),
        "Reason" = reason,
        "Excluded" = n -(table |> nrow())
      )
    ) |> 
    distinct()
  return(table)
}
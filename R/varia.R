#' Export list of dataframes as multiple sheets of xlsx file
#'
#' @name export_mass_xlsx
#' @param list_of_dfs list of dataframes to export
#' @param xlsx_name character string specyfying name of exported .xlsx
#'
#' @import xlsx
#' 
#' @export
#'

export_mass_xlsx <- function(list_of_dfs, xlsx_name) {

  write.xlsx(x = as.data.frame(list_of_dfs[[1]]),
             file = as.character(xlsx_name),
             sheetName = as.character(names(list_of_dfs)[1]))

  for(index in c(1:(length(list_of_dfs) - 1))) {

    write.xlsx(x = as.data.frame(list_of_dfs[[index+1]]),
               file = as.character(xlsx_name),
               sheetName = as.character(names(list_of_dfs)[index+1]),
               append = TRUE)

  }

}

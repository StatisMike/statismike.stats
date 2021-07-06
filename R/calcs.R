#### libraries ####

library(tidyverse)
library(psych)
library(lavaan)
library(officer)

#### functions ####
#' Calc Cronbach's alpha for multiple scales
#'
#' @name calc_alphas_for_scales
#'
#' @param data raw data as data.frame object
#' @param scales data.frame object describing scales and which items it consists
#' of. The structure is : \code{data.frame(scale = "item_1, -item_2, ... , item_n")},
#' where the negatives means that item is need to be reversed
#' @param values \code{c(min, max)} value of item
#' @param conf.level the confidence level for confidence intervals
#'
#' @return A data.frame containing number of items, raw Cronbachs alphas, standarized
#' alphas and lower and upper CI for alphas
#'
#' @import tidyverse
#' @import psych
#' @import cocron
#'

calc_alphas_for_scales <- function(data, scales, values = c(1,5), conf.level = 0.95) {

  dataframe <- data
  scalenames <- names(scales)
  datascale <- scales
  #  dataprefix <- prefix
  min_value <- values[1]
  max_value <- values[2]

  all_alphas <- data.frame(n = 0,
                           raw_alpha = 0,
                           std.alpha = 0,
                           lower_ci = 0,
                           upper_ci = 0)

  for (scale in scalenames) {

    items <- unlist(strsplit(as.character(datascale[scale]), split = ", ", fixed = TRUE))

    to_reverse <- ifelse(grepl(items, pattern = "-", fixed = TRUE) == TRUE, TRUE, FALSE)

    items_names <- str_remove(items, "-")

    to_reverse_names <- items_names[to_reverse]

    raw_values_straight <- select(dataframe, all_of(items_names), -all_of(to_reverse_names))

    raw_values_reverse <- select(dataframe, all_of(to_reverse_names)) %>%
      transmute_all(function(x){max_value + min_value - x})

    raw_values_all <- cbind(raw_values_straight, raw_values_reverse)

    y <- psych::alpha(raw_values_all)

    temp_alpha <- y$total[1:2]
    temp_ci <- cocron::cronbach.alpha.CI(alpha = y$total[1, 2], n = length(y$scores), items = nrow(y$item.stats), conf.level = conf.level)
    output1 <- data.frame(n = nrow(y$item.stats),
                          raw_alpha = temp_alpha[1],
                          std.alpha = temp_alpha[2],
                          lower_ci = temp_ci[1],
                          upper_ci = temp_ci[2])

    row.names(output1) <- scale

    all_alphas <- rbind(all_alphas, output1)

  }

  all_alphas <- all_alphas[-1,]
  return(all_alphas)

}

#' Calculate sum/or means for multiple scales of questionnaire
#'
#' @name calc_items_to_scale
#'
#' @param data raw data as a data.frame object
#' @param scales data.frame object describing scales and which items it consists
#' of. The structure is : \code{data.frame(scale = "item_1, -item_2, ... , item_n")},
#' where the negatives means that item is need to be reversed
#' @param values \code{c(min, max)} value of item
#' @param mean should the function calculate means of scales? If \code{FALSE}, then
#' the sums will be calculated
#' @param na.rm remove missing values?
#'
#' @return data.frame containing the summarized scales or mean of items in scales
#'
#' @import tidyverse
#'

calc_items_to_scale <- function(data, scales, values = c(1,5), mean = TRUE, na.rm = TRUE) {

  require(tidyverse)

  dataframe <- data
  scalenames <- names(scales)
  datascale <- scales
  #  dataprefix <- prefix
  min_value <- values[1]
  max_value <- values[2]

  scale_scores <- data.frame(id = c(1:nrow(dataframe)))

  for (scale in scalenames) {

    items <- unlist(strsplit(as.character(datascale[1, scale]), split = ", ", fixed = TRUE))

    to_reverse <- ifelse(grepl(items, pattern = "-", fixed = TRUE) == TRUE, TRUE, FALSE)

    items_names <- str_remove(items, "-")

    to_reverse_names <- items_names[to_reverse]

    raw_values_straight <- select(dataframe, all_of(items_names), -all_of(to_reverse_names))

    raw_values_reverse <- select(dataframe, all_of(to_reverse_names)) %>%
      transmute_all(function(x){max_value + min_value - x})

    raw_values_all <- cbind(raw_values_straight, raw_values_reverse)

    if(mean == TRUE) {
      scale_scores <- cbind(scale_scores, rowMeans(raw_values_all, na.rm = na.rm))
      names(scale_scores)[ncol(scale_scores)] <- scale
    } else {
      scale_scores <- cbind(scale_scores, rowSums(raw_values_all, na.rm = na.rm))
      names(scale_scores)[ncol(scale_scores)] <- scale
    }

  }

  return(scale_scores)

}

#' Summarize descriptives: mean and standard deviation for grouped variable
#'
#' @name summarize_descriptives
#'
#' @param data data.frame
#' @param var name of the variable for which to calculate descriptives
#' @param grouping.var names of the variables which should be used for grouping,
#' maximum 2
#'
#' @import tidyverse
#'

summarize_descriptives <- function(data, var, grouping.var) {
  require(tidyverse)

  if(length(grouping.var) == 1) {

    columns <- c(var, grouping.var)

    summ.data <- data %>%
      select(all_of(columns)) %>%
      rename(Variable = 1, G.Var = 2) %>%
      group_by(G.Var) %>%
      summarize(N = n(), percent = n()/nrow(data), M = mean(Variable, na.rm = T), SD = sd(Variable, na.rm = T)) %>%
      ungroup()

    colnames(summ.data)[1] <- paste(grouping.var)
    colnames(summ.data)[4] <- paste("M", var, sep = " ")
    colnames(summ.data)[5] <- paste("SD", var, sep = " ")

    return(summ.data)


  } else if(length(grouping.var) == 2){

    columns <- c(var, grouping.var)

    summ.data <- data %>%
      select(all_of(columns)) %>%
      rename(Variable = 1, G.Var.One = 2, G.Var.Two = 3) %>%
      group_by(G.Var.One, G.Var.Two) %>%
      summarize(N = n(), percent = n()/nrow(data), M = mean(Variable, na.rm = T), SD = sd(Variable, na.rm = T))

    colnames(summ.data)[1] <- paste(grouping.var[1])
    colnames(summ.data)[2] <- paste(grouping.var[2])
    colnames(summ.data)[5] <- paste("M", var, sep = " ")
    colnames(summ.data)[6] <- paste("SD", var, sep = " ")

    return(summ.data)

  } else  if(length(grouping.var) > 2){

    return("More grouping variables with with one variable descriptive statistics are not supported at this moment")

  } else {return("Error")}

}

#' Calculate scale descriptives for one or more variables with or without grouping
#'
#' @name calc_scale_desc
#'
#' @param data data.frame object
#' @param Vars.list character vector of variable names, for which the scale
#' descriptives will be calculated
#' @param grouping name of the grouping variable. defaults to "no grouping", where
#' there will not be any grouping
#' @param not.normal.list if \code{TRUE}, then output will be a character vector
#' with names of not normal distributed variables
#' @param normal.list if \code{TRUE}, then output will be a character vector with
#' names of normally distributed variable
#' @param na.rm remove missing values
#' @param normality_test which normality test to use for checking normality of
#' distribution. Defaults to \code{"shapiro"}, can also use \code{"lillie"} for
#' Lilliefors test.
#'
#' @import tidyverse
#' @import broom
#' @import moments
#' @import nortest
#'

calc_scale_desc <- function(data, Vars.list, grouping = "no grouping", not.normal.list = FALSE, normal.list = FALSE, na.rm = FALSE, normality_test = "shapiro") {

  # for calculating with grouping variable

  if(length(grouping) == 2) {#making long data

    data_long <- data %>%
      select(Vars.list, var2 = grouping[1], var3 = grouping[2]) %>%
      gather(key = "Variable", value = "val", -var2, -var3)

    #computing normality test

    moreThan3 <- data_long %>%
      group_by(var2, var3, Variable) %>%
      summarize(n = n()) %>%
      filter(n > 3) %>%
      select(Variable) %>%
      ungroup()

    if(normality_test == "shapiro") {

      normality_test <- data_long %>%
        filter(var2 %in% moreThan3$var2 & var3 %in% moreThan3$var3) %>%
        group_by(var2,var3, Variable) %>%
        do(tidy(shapiro.test(.$val))) %>%
        ungroup()

    }
    else if(normality_test == "lillie") {

      normality_test <- data_long %>%
        filter(var2 %in% moreThan3$var2 & var3 %in% moreThan3$var3) %>%
        group_by(var2, var3, Variable) %>%
        do(tidy(lillie.test(.$val))) %>%
        ungroup()

    }

    #computing descriptives: mean, median, sd, min, max, skewness and kurtosis

    skewness <- data_long %>%
      group_by(var2, var3, Variable) %>%
      summarize(n = n(),
                mean = mean(val, na.rm = na.rm),
                median = median(val, na.rm = na.rm),
                min = min(val, na.rm = na.rm),
                max = max(val, na.rm = na.rm),
                sd = sd(val, na.rm = na.rm),
                skewness = skewness(val, na.rm = na.rm),
                kurtosis = kurtosis(val, na.rm = na.rm)) %>%
      ungroup()

    #joining it together

    descriptives_full <- full_join(skewness, normality_test, by = c("var2", "var3", "Variable"))

    #changing back the column name of grouping variable to its original name

    colnames(descriptives_full)[1] <- grouping[1]
    colnames(descriptives_full)[2] <- grouping[2]

  }

  else if(grouping == "no grouping") {

    #calculating without grouping variable
    #creating long data

    data_long <- data %>%
      select(Vars.list) %>%
      gather(key = "Variable", value = "val")

    #Shapiro-Wilk

    if(normality_test == "shapiro") {
      normality_test <- data_long %>%
        group_by(Variable) %>%
        do(tidy(shapiro.test(.$val))) %>%
        ungroup() } else if(normality_test == "lillie") {

          normality_test <- data_long %>%
            group_by(Variable) %>%
            do(tidy(lillie.test(.$val))) %>%
            ungroup()
        }

    #descriptives

    skewness <- data_long %>%
      group_by(Variable) %>%
      summarize(n = n(),
                mean = mean(val, na.rm = na.rm),
                median = median(val, na.rm = na.rm),
                min = min(val, na.rm = na.rm),
                max = max(val, na.rm = na.rm),
                sd = sd(val, na.rm = na.rm),
                skewness = skewness(val, na.rm = na.rm),
                kurtosis = kurtosis(val, na.rm = na.rm)) %>%
      ungroup()

    #joining it all together

    descriptives_full <- inner_join(skewness, normality_test, by = "Variable")

  }

  else {

    #making long data

    data_long <- data %>%
      select(all_of(Vars.list), var2 = grouping) %>%
      gather(key = "Variable", value = "val", -var2)

    #computing normality test

    moreThan3 <- data_long %>%
      group_by(var2, Variable) %>%
      summarize(n = n()) %>%
      filter(n > 3) %>%
      select(Variable) %>%
      ungroup()

    if(normality_test == "shapiro") {

      normality_test <- data_long %>%
        filter(var2 %in% moreThan3$var2) %>%
        group_by(var2, Variable) %>%
        do(tidy(shapiro.test(.$val))) %>%
        ungroup()

    }
    else if(normality_test == "lillie") {

      normality_test <- data_long %>%
        filter(var2 %in% moreThan3$var2) %>%
        group_by(var2, Variable) %>%
        do(tidy(lillie.test(.$val))) %>%
        ungroup()

    }

    #computing descriptives: mean, median, sd, min, max, skewness and kurtosis

    skewness <- data_long %>%
      group_by(var2, Variable) %>%
      summarize(n = n(),
                mean = mean(val, na.rm = na.rm),
                median = median(val, na.rm = na.rm),
                min = min(val, na.rm = na.rm),
                max = max(val, na.rm = na.rm),
                sd = sd(val, na.rm = na.rm),
                skewness = skewness(val, na.rm = na.rm),
                kurtosis = kurtosis(val, na.rm = na.rm)) %>%
      ungroup()

    #joining it together

    descriptives_full <- full_join(skewness, normality_test, by = c("var2", "Variable"))

    #changing back the column name of grouping variable to its original name

    colnames(descriptives_full)[1] <- grouping}




  #if you don't specify not.normal.list, you get a table of descriptives and
  #and normality test for your variables (in groups if specified)

  if(not.normal.list == FALSE) {

    if(normal.list == TRUE) {

      list <- descriptives_full %>%
        filter(p.value >= 0.05) %>%
        select(Variable)

      return(list$Variable) } else {

        return(descriptives_full) }

  } else {

    #if you specify not.normal.list = TRUE, you get a vector instead
    #You can use it after to specify variables you DON'T want to use
    #in parametric tests (and WANT to use in non-parametrics)

    list <- descriptives_full %>%
      filter(p.value < 0.05) %>%
      select(Variable)

    return(list$Variable)

  }

}


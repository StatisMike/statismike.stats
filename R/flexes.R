#' Generate Wilcoxon test of differences between two groups and present results as a flextable
#'
#' @name flex_wilcox
#'
#' @param vars names of columns in data to test for differences
#' @param grouping name of column in data to distinguish groups
#' @param vars_names names of variables to print in flextable
#' @param grouping_names name of groups
#' @param group_name name of grouping in header
#' @param lang default to "english". Can also be "polish"
#' @param alternative either "two.sided" (default), "greater" or "less"
#' @param h_lines if there should be horizontal lines between variables
#'
#' @return flextable
#'
#' @import tidyverse
#' @import rstatix
#' @import flextable
#'

flex_wilcox <- function(data,
                        vars,
                        grouping,
                        vars_names = NULL,
                        grouping_names = NULL,
                        group_name = NULL,
                        lang = "english",
                        alternative = "two.sided",
                        h_lines = FALSE) {

  require(tidyverse)
  require(rstatix)
  require(flextable)

  sig_check <- function(x) {
    ifelse(x >= 0.1, " ",
           ifelse(x >= 0.05, ".",
                  ifelse(x >= 0.01, "*",
                         ifelse(x >= 0.001, "**", "***"))))
  }

  data <- data %>%
    select(all_of(grouping), all_of(vars))

  names(data)[1] <- "group"

  data <- data %>%
    gather(key = "Variable", value = "Value", -group)

  if(is.null(vars_names) == T){

    data$Variable <- factor(data$Variable, levels = vars)

  } else {

    data$Variable <- factor(data$Variable, levels = vars)
    data$Variable <- factor(data$Variable, labels = vars_names)

    vars <- vars_names

  }

  if(is.null(grouping_names) == F) {

    data$group <- factor(data$group,
                         labels = grouping_names)

  } else {

    data$group <- as.factor(data$group)

  }

  test_results <- data %>%
    group_by(Variable) %>%
    wilcox_test(Value ~ group, alternative = alternative) %>%
    ungroup()

  test_effsize <- data %>%
    group_by(Variable) %>%
    wilcox_effsize(Value ~ group, alternative = alternative) %>%
    ungroup()

  test_desc <- data %>%
    group_by(Variable, group) %>%
    summarize(mean = mean(Value, na.rm = T), median = median(Value, na.rm = T), .groups = "keep") %>%
    ungroup()

  test_desc["diff_M"] <- rep(0, times = nrow(test_desc))
  test_desc["diff_Me"] <- rep(0, times = nrow(test_desc))

  n_iter <- nrow(test_desc)

  for(var in 1:(n_iter/2)){

    mean_diff <- abs(test_desc$mean[1 + 2*(var - 1)] - test_desc$mean[2 + 2*(var - 1)])
    median_diff <- abs(test_desc$median[1 + 2*(var - 1)] - test_desc$median[2 + 2*(var - 1)])

    test_desc[(1 + 2*(var - 1)):(2 + 2*(var - 1)), 5] <- mean_diff
    test_desc[(1 + 2*(var - 1)):(2 + 2*(var - 1)), 6] <- median_diff


  }

  test_desc <- test_desc[1:n_iter,]

  test_wilc <- full_join(test_results, test_effsize, by = "Variable") %>%
    mutate(signif = sig_check(p))

  test_renamed <- test_desc %>%
    left_join(test_wilc[,c(1, 7:8, 16, 12)], by = "Variable")

  test_renamed <- mutate_at(test_renamed,
                            .vars = c(3:7, 10),
                            .funs = formatC,
                            digits = 2,
                            format = "f") %>%
    mutate(p = formatC(p, digits = 3, format = "f"))

  out <- flextable(test_renamed)

  header_names <- if(lang == "english") {

    data_frame(col_keys = c("Variable", "group", "mean", "median",
                            "diff_M", "diff_Me", "statistic", "p",
                            "signif", "effsize"),
               names = c("Variable",
                         if_else(is.null(group_name), "Group", group_name),
                         "Mean", "Median", "Difference of means",
                         "Difference of medians", "U", "p", "p", "r"))


  } else {

    data_frame(col_keys = c("Variable", "group", "mean", "median",
                            "diff_M", "diff_Me", "statistic", "p",
                            "signif", "effsize"),
               names = c("Zmienna",
                         if_else(is.null(group_name), "Grupa", group_name),
                         "Średnia", "Mediana", "Różnica średnich",
                         "Różnica median", "U", "p", "p", "r"))

  }

  out <- out %>%
    set_header_df(mapping = header_names, key = "col_keys") %>%
    merge_at(j = 8:9, i = 1, part = "header") %>%
    theme_booktabs() %>%
    align(align = "center", part = "all")


  for(var in c(1:length(vars))){

    join_vert <- c(1+2*(var-1), 2+2*(var-1))

    out <- merge_at(out, i = join_vert, j = 1, part = "body") %>%
      merge_at(i = join_vert, j = 5, part = "body") %>%
      merge_at(i = join_vert, j = 6, part = "body") %>%
      merge_at(i = join_vert, j = 7, part = "body") %>%
      merge_at(i = join_vert, j = 8, part = "body") %>%
      merge_at(i = join_vert, j = 9, part = "body") %>%
      merge_at(i = join_vert, j = 10, part = "body")

    if(h_lines == TRUE){

      out <- hline(out, i = join_vert[2], border = fp_border(width = 1))

    } else { }


  }

  out <- out %>%
    width(j = 8, width = 0.80) %>%
    width(j = 9, width = 0.30) %>%
    align(j = 8, align = "right", part = "body") %>%
    align(j = 9, align = "left", part = "body") %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "all") %>%
    hline(border = fp_border(width = 2), i = length(vars)*2, part = "body") %>%
    fix_border_issues(part = "all") %>%
    padding(j = 8:9, padding = 1)

  if(lang == "english"){
    out <- width(out, j = 6, width = 1)

  }


  return(out)


}

#' Calculate correlations and present in flextable. Automatically tests normality to choose between Pearson or Spearman
#'
#' @name flex_corr
#'
#' @param vars_hor names of columns in data on horizontal
#' @param vars_ver names of columns in data on vertical
#' @param vars_hor_names names vars_hor to print
#' @param vars_ver_names names of vars_ver to print
#' @param normality_test which normality test to use? Either "shapiro" or "lillie"
#' @param adjust adjustement for multiple tests
#' @param lang default to "english". Can also be "polish"
#' @param p_values include p values?
#' @param col_widths how width should the cols be? If not default, then specify
#' list with numeric values for \code{r_width}, \code{m_width}, \code{p_width} and \code{s_width}
#'
#' @return flextable
#'
#' @import tidyverse
#' @import broom
#' @import flextable
#' @import psych
#' @import nortest
#'

flex_corr <- function(data,
                      vars_hor,
                      vars_ver,
                      vars_hor_names = vars_hor,
                      vars_ver_names = vars_ver,
                      normality_test = "shapiro",
                      adjust = "none",
                      lang = "english",
                      p_values = "include",
                      col_widths = "default") {

  #### CORR CALCULATION ####

  if(col_widths == "default"){

    r_width = 0.75
    m_width = 0.25
    p_width = 1
    s_width = 0.25

  } else {

    if(c("r_width", "m_width", "p_width", "s_width") %in% names(col_widths)){

      r_width = col_widths[["r_width"]]
      m_width = col_widths[["m_width"]]
      p_width = col_widths[["p_width"]]
      s_width = col_widths[["s_width"]]

    } else {

      warning("If the col_widths isn't 'default', it should be a list of four elements: r_width, m_width, p_width and s_width, each a numeric variable.")

    }

  }

  data1 <- select(data, all_of(vars_hor))
  data2 <- select(data, all_of(vars_ver))

  data_long1 <- data1 %>% gather(key = "Variable", value = "val")
  if (normality_test == "shapiro") {
    normality_test1 <- data_long1 %>% group_by(Variable) %>%
      do(tidy(shapiro.test(.$val))) %>% ungroup()
  }
  else if (normality_test == "lillie") {
    normality_test1 <- data_long1 %>% group_by(Variable) %>%
      do(tidy(lillie.test(.$val))) %>% ungroup()
  }
  normal_list1 <- normality_test1 %>% filter(p.value >= 0.05) %>%
    select(Variable) %>% unlist()
  data_long2 <- data2 %>% gather(key = "Variable", value = "val")
  if (normality_test == "shapiro") {
    normality_test2 <- data_long2 %>% group_by(Variable) %>%
      do(tidy(shapiro.test(.$val))) %>% ungroup()
  }
  else if (normality_test == "lillie") {
    normality_test2 <- data_long2 %>% group_by(Variable) %>%
      do(tidy(lillie.test(.$val))) %>% ungroup()
  }
  normal_list2 <- normality_test2 %>% filter(p.value >= 0.05) %>%
    select(Variable) %>% unlist()
  pearson_corr <- corr.test(data1, data2, method = "pearson",
                            adjust = adjust)
  pearson_corr_r <- as.data.frame(pearson_corr$r)
  pearson_corr_p <- as.data.frame(pearson_corr$p)
  spearman_corr <- corr.test(data1, data2, method = "spearman",
                             adjust = adjust)
  spearman_corr_r <- as.data.frame(spearman_corr$r)
  spearman_corr_p <- as.data.frame(spearman_corr$p)
  corr_all <- data.frame(id = c(1:nrow(pearson_corr_r)))
  sig_check <- function(x, p) {
    ifelse(x[[p]] >= 0.1, " ", ifelse(x[[p]] >= 0.05,
                                      ".", ifelse(x[[p]] >= 0.01, "*", ifelse(x[[p]] >=
                                                                                0.001, "**", "***"))))
  }

  #### Creating whole table ####

  for (col in c(1:ncol(pearson_corr_r))) {
    if (names(pearson_corr_r)[col] %in% normal_list2) {
      row_temp <- data.frame(r = 0, m = "", p = 0, s = "")

      names(row_temp)[1] <- paste("r", names(pearson_corr_r)[col],
                                  sep = "_")
      names(row_temp)[3] <- paste("p", names(pearson_corr_p)[col],
                                  sep = "_")
      names(row_temp)[4] <- paste("s", names(pearson_corr_p)[col],
                                  sep = "_")
      names(row_temp)[2] <- paste("m", names(pearson_corr_p)[col],
                                  sep = "_")

      for (row in c(1:nrow(pearson_corr_r))) {
        if (row.names(pearson_corr_r)[row] %in% normal_list1) {
          pearson_corr_temp <- data.frame(r = pearson_corr_r[row,
                                                             col], p = pearson_corr_p[row, col]) %>%
            mutate(s = sig_check(., 2), m = "P") %>%
            select(r, m, p, s)

          names(pearson_corr_temp)[(1:4)] <- c(paste("r",
                                                     names(pearson_corr_r)[col], sep = "_"),
                                               paste("m",
                                                     names(pearson_corr_p)[col], sep = "_"),
                                               paste("p",
                                                     names(pearson_corr_p)[col], sep = "_"),
                                               paste("s",
                                                     names(pearson_corr_p)[col], sep = "_"))

          row.names(pearson_corr_temp)[nrow(pearson_corr_temp)] <- row.names(pearson_corr_r)[row]

          row_temp <- rbind(row_temp, pearson_corr_temp)
        }

        else {

          spearman_corr_temp <- data.frame(r = spearman_corr_r[row,
                                                               col], p = spearman_corr_p[row, col]) %>%
            mutate(s = sig_check(., 2), m = "S")%>%
            select(r, m, p, s)

          names(spearman_corr_temp)[(1:4)] <- c(paste("r",
                                                      names(spearman_corr_r)[col], sep = "_"),
                                                paste("m",
                                                      names(spearman_corr_p)[col], sep = "_"),
                                                paste("p",
                                                      names(spearman_corr_p)[col], sep = "_"),
                                                paste("s",
                                                      names(spearman_corr_p)[col], sep = "_"))

          row.names(spearman_corr_temp)[nrow(spearman_corr_temp)] <- row.names(spearman_corr_r)[row]

          row_temp <- rbind(row_temp, spearman_corr_temp)
        }
      }
      corr_all <- cbind(corr_all, row_temp[-1, ])
    }
    else {
      row_temp <- data.frame(r = 0, m = "", p = 0, s = "")

      names(row_temp)[1] <- paste("r", names(pearson_corr_r)[col],
                                  sep = "_")
      names(row_temp)[3] <- paste("p", names(pearson_corr_p)[col],
                                  sep = "_")
      names(row_temp)[4] <- paste("s", names(pearson_corr_p)[col],
                                  sep = "_")
      names(row_temp)[2] <- paste("m", names(pearson_corr_p)[col],
                                  sep = "_")
      for (row in c(1:nrow(pearson_corr_r))) {

        spearman_corr_temp <- data.frame(r = spearman_corr_r[row,
                                                             col], p = spearman_corr_p[row, col]) %>%
          mutate(s = sig_check(., 2), m = "S")%>%
          select(r, m, p, s)

        names(spearman_corr_temp)[(1:4)] <- c(paste("r",
                                                    names(spearman_corr_r)[col], sep = "_"),
                                              paste("m",
                                                    names(spearman_corr_p)[col], sep = "_"),
                                              paste("p",
                                                    names(spearman_corr_p)[col], sep = "_"),
                                              paste("s",
                                                    names(spearman_corr_p)[col], sep = "_"))

        row.names(spearman_corr_temp)[nrow(spearman_corr_temp)] <- row.names(spearman_corr_r)[row]

        row_temp <- rbind(row_temp, spearman_corr_temp)
      }
      corr_all <- cbind(corr_all, row_temp[-1, ])
    }
  }

  corr_all <- corr_all[,-1] %>%
    mutate(
      across(starts_with("r_"),
             round,
             digits = 2) ) %>%
    mutate(
      across(starts_with("r_"),
             sprintf,
             fmt = "%0.2f",
             how = "replace") ) %>%

    mutate(
      across(starts_with("p_"),
             round,
             digits = 3) )  %>%

    mutate(
      across(starts_with("p_"),
             sprintf,
             fmt = "%0.3f",
             how = "replace") )  %>%

    mutate(
      across(.fns = as.character)
    )

  if(lang == "english"){

    corr_all <- cbind(tibble(Variable = vars_hor_names),
                      corr_all)

  } else {

    corr_all <- cbind(tibble(Zmienna = vars_hor_names),
                      corr_all)

  }

  corr_all[corr_all == "0.000"] <- "< 0.001"


  main_header_vars <- data_frame(col_keys = names(corr_all),
                                 names = c(names(corr_all[1]), rep("", times = length(names(corr_all))-1)))

  for(var in c(1:length(vars_ver))){

    namenum <- c(2, 3, 4, 5) + (var - 1)*4

    main_header_vars[namenum, "names"] <- vars_ver_names[var]

  }

  normal_n <- tibble(S = sum(corr_all == "S"),
                     P = sum(corr_all == "P"))

  if(normal_n$S > normal_n$P){

    sub_header_vars <- c(names(corr_all)[1], rep(c("ρ", "ρ", "p", "p"), times = length(vars_ver)))

  } else {

    sub_header_vars <- c(names(corr_all)[1], rep(c("r", "r", "p", "p"), times = length(vars_ver)))

  }

  if(normal_n$S == 0){

    corr_all[corr_all == "P"] <- " "

  } else if(normal_n$P == 0 | normal_n$S > normal_n$P){

    corr_all[corr_all == "S"] <- " "

  } else {

  }

  if(p_values == "include"){

    out <- flextable(corr_all) %>%
      set_header_df(mapping = main_header_vars, key = "col_keys") %>%
      add_header_row(values = sub_header_vars, top = F) %>%
      merge_at(i = 1:2, j = 1, part = "header") %>%
      theme_booktabs() %>%
      hline_top(border = fp_border(width = 2), part = "all") %>%
      hline_bottom(border = fp_border(width = 2), part = "all")

    for(var in c(1:length(vars_ver))){

      namenum <- c(2, 3, 4, 5) + (var - 1)*4

      out <- merge_at(out, i = 1, j = namenum, part = "header")           # main header merge
      out <- merge_at(out, i = 2, j = namenum[1:2], part = "header")      # r / rho merge
      out <- merge_at(out, i = 2, j = namenum[3:4], part = "header")      # p merge
      out <- align(out, i = c(1:2), namenum, part = "header",
                   align = "center")
      out <- align(out, i = c(1:length(vars_hor)),
                   j = namenum[c(2, 4)], align = "left", part = "body")   # symbols alignment
      out <- align(out, i = c(1:length(vars_hor)),
                   j = namenum[c(1, 3)], align = "right", part = "body")  # values alignment
      out <- width(out, j = namenum[1], width = r_width)                        # r values width
      out <- width(out, j = namenum[2], width = m_width)                     # method symbols width
      out <- width(out, j = namenum[3], width = p_width)                        # p values width
      out <- width(out, j = namenum[4], width = s_width)                     # signif width
      out <- padding(out, j = namenum[1], padding.right = 1) %>%
        padding(j = namenum[3], padding.right = 1) %>%
        padding(j = namenum[2], padding.left = 1) %>%
        padding(j = namenum[4], padding.left = 1)


    }

  } else {

    main_header_vars <- data_frame(col_keys =c(if_else(lang == "english", "Variable", "Zmienna"),
                                               rep("", times = length(vars_ver) * 3 - 1)),
                                   names = c(if_else(lang == "english", "Variable", "Zmienna"),
                                             rep("", times = length(vars_ver) * 3 - 1)))

    var_n <- 0

    for(var in vars_ver_names){

      var_n <- var_n + 1

      namenum <- c(2, 5, 3) + (var_n - 1)*4

      namenum_x <- c(2, 3, 4) + (var_n - 1)*3

      main_header_vars[namenum_x, 1] <- names(corr_all)[namenum]

      main_header_vars[namenum_x, 2] <- vars_ver_names[var_n]


    }

    out <- flextable(corr_all,
                     col_keys = main_header_vars$col_keys) %>%
      set_header_df(mapping = main_header_vars, key = "col_keys") %>%
      theme_booktabs() %>%
      hline_top(border = fp_border(width = 2), part = "all") %>%
      hline_bottom(border = fp_border(width = 2), part = "all")

    for(var in c(1:length(vars_ver))){

      namenum <- c(2, 3, 4) + (var - 1)*3

      out <- merge_at(out, i = 1, j = namenum, part = "header")           # main header merge
      out <- align(out, i = 1, namenum, part = "header",
                   align = "center")
      out <- align(out, i = c(1:length(vars_hor)),
                   j = namenum[c(2, 3)], align = "left", part = "body")   # symbols alignment
      out <- align(out, i = c(1:length(vars_hor)),
                   j = namenum[1], align = "right", part = "body")        # values alignment
      out <- width(out, j = namenum[1], width = r_width)                        # r values width
      out <- width(out, j = namenum[3], width = m_width)                     # method symbols width
      out <- width(out, j = namenum[2], width = s_width)                     # signif width
      out <- padding(out, j = namenum[1], padding.right = 1) %>%
        padding(j = namenum[2], padding.right = 1, padding.left = 1) %>%
        padding(j = namenum[3], padding.left = 1)

    }

  }

  footer_value <- if(normal_n$S == 0 | normal_n$P == 0 ){

    if(lang == "english") {

      "* p < .05. ** p < .01. and *** p < .001. Values of p < 0.1 are marked by a comma."

    } else {

      "* p < .05. ** p < .01. a *** p < .001. Wartości p < 0.1 oznaczone są kropką."

    }

  } else {

    if(lang == "english") {

      foot1 <- "* p < .05. ** p < .01. and *** i p < .001. Values of p < 0.1 are marked by a comma."
      foot2 <- if_else(normal_n$S > normal_n$P,
                       "Letter 'P' next to the value of correlation estimate means usage of Pearson test. Otherwise, Spearman test was used.",
                       "Letter 'S' next to the value of correlation estimate means usage os Spearman test. Otherwise, Pearson test was used.")
      paste(foot1, foot2, sep = " ")

    } else {

      foot1 <- "* p < .05. ** p < .01. a *** p < .001. Wartości p < 0.1 oznaczone są kropką."
      foot2 <- if_else(normal_n$S > normal_n$P,
                       "Litera 'P' przy współczynniku korelacji oznacza wykorzystanie testu Pearsona. W pozostałych przypadkach, wykorzystano test Spearmana",
                       "Litera 'S' przy współczynniku korelacji oznacza wykorzystanie testu Spearmana. W pozostałych przypadkach, wykorzystano test Pearsona")
      paste(foot1, foot2, sep = " ")

    }




  }

  out <- add_footer_lines(out, values = footer_value)

  # corr_all

}

#' Generate t Student tests of differences between two groups and present results as a flextable. Automatically applies Welsh where Levene test is significant.
#'
#' @name flex_t_test
#'
#' @param vars names of columns in data to test for differences
#' @param grouping name of column in data to distinguish groups
#' @param vars_names names of variables to print in flextable
#' @param grouping_names name of groups
#' @param group_name name of grouping in header
#' @param lang default to "english". Can also be "polish"
#' @param alternative either "two.sided" (default), "greater" or "less"
#' @param h_lines if there should be horizontal lines between variables
#'
#' @import dplyr
#' @import car
#' @import broom
#' @import rstatix
#' @import flextable
#'

flex_t_test <- function(data,
                        vars,
                        grouping,
                        vars_names = NULL,
                        grouping_names = NULL,
                        group_name = NULL,
                        lang = "english",
                        alternative = "two.sided",
                        h_lines = FALSE) {

  sig_check <- function(x) {
    ifelse(x >= 0.1, " ",
           ifelse(x >= 0.05, ".",
                  ifelse(x >= 0.01, "*",
                         ifelse(x >= 0.001, "**", "***"))))
  }


  cond <- "Both"

  data_renamed <- data %>% select(all_of(vars), groups = all_of(grouping)) %>%
    gather(key = "var", value = "value", -groups)

  if(is.null(vars_names) == T){

    data_renamed$var <- factor(data_renamed$var, levels = vars)

  } else {

    data_renamed$var <- factor(data_renamed$var, levels = vars)
    data_renamed$var <- factor(data_renamed$var, labels = vars_names)

    vars <- vars_names

  }


  if(is.null(grouping_names) == F) {

    data_renamed$groups <- factor(data_renamed$groups,
                                  labels = grouping_names)

  } else {

    data_renamed$groups <- as.factor(data_renamed$groups)

  }

  Levene_test <- data_renamed %>%
    select(var, value, groups) %>%
    filter(var %in% vars) %>%
    group_by(var) %>%
    do(tidy(car::leveneTest(.$value ~ .$groups))) %>%
    ungroup()

  homo.variance <- Levene_test %>%
    filter(p.value > 0.05) %>%
    select(var)

  hetero.variance <- Levene_test %>%
    filter(p.value < 0.05) %>%
    select(var)

  if (nrow(homo.variance) != 0) {
    t_test <- data_renamed %>%
      select(var, value, groups) %>%
      filter(var %in% homo.variance$var) %>%
      group_by(var) %>%
      do(tidy(t.test(.$value ~ .$groups, var.equal = TRUE, alternative = alternative))) %>%
      mutate(estimate = estimate1 - estimate2) %>%
      ungroup()

    coh_d <- data_renamed %>%
      filter(var %in% homo.variance$var) %>%
      select(var, value, groups) %>%
      group_by(var) %>%
      do(rstatix::cohens_d(., value ~ groups,
                           paired = FALSE,
                           var.equal = TRUE)) %>%
      ungroup()

    t_test <- inner_join(t_test, coh_d, by = "var") %>%
      mutate(signif = sig_check(p.value)) %>%
      select(Variable = var,
             Difference = estimate,
             statistic,
             df = parameter,
             p.value,
             signif,
             CI_lower = conf.low,
             CI_upper = conf.high,
             method = method,
             d_Cohen = effsize) %>%
      mutate(method = if_else(lang == "english", "Student test", "Test Studenta"),
             Difference = abs(Difference))
  } else {
    cond <- "OnlyHetero"
  }


  if (nrow(hetero.variance) != 0) {
    t_test_Welch <- data_renamed %>%
      select(var, value, groups) %>%
      filter(var %in% hetero.variance$var) %>%
      group_by(var) %>%
      do(tidy(t.test(.$value ~ .$groups, var.equal = FALSE, alternative = alternative))) %>%
      ungroup()

    coh_d_Welch <- data_renamed %>%
      filter(var %in% hetero.variance$var) %>%
      select(var, value, groups) %>% group_by(var) %>%
      do(rstatix::cohens_d(., value ~ groups, paired = FALSE,
                           var.equal = FALSE)) %>%
      ungroup()

    t_test_Welch <- inner_join(t_test_Welch, coh_d_Welch,
                               by = "var") %>%
      mutate(signif = sig_check(p.value)) %>%
      select(Variable = var,
             Difference = estimate,
             statistic,
             df = parameter,
             p.value,
             signif,
             CI_lower = conf.low,
             CI_upper = conf.high,
             method = method,
             d_Cohen = effsize) %>%
      mutate(method = if_else(lang == "english", "Welch test", "Test Welcha"),
             Difference = abs(Difference))
  } else {
    cond <- "OnlyHomo"
  }

  if (cond == "OnlyHomo") {
    t_test_out <- t_test

  } else {

    if (cond == "OnlyHetero") {

      t_test_out <- t_test_Welch

    } else {

      t_test_out <- rbind(t_test, t_test_Welch)

    }
  }

  means <- data_renamed %>%
    group_by(var, groups) %>%
    summarize(.groups = "keep", means = mean(value, na.rm = T))

  out <- left_join(means, t_test_out, by = c("var" = "Variable")) %>%
    arrange(var, all_of(grouping))

  out <- out %>%

    mutate(across(.cols = c("means", "Difference", "statistic", "df", "CI_lower", "CI_upper", "d_Cohen"),
                  .fns = formatC,
                  digits = 2,
                  format = "f")) %>%
    mutate(p.value = formatC(p.value, digits = 3, format = "f"))

  if(alternative == "greater") {

    out$CI_upper <- rep("∞",
                        times = length(vars)*2)

  } else if(alternative == "less") {

    out$CI_lower <- rep("-∞",
                        times = length(vars)*2)

  } else {

  }

  if(lang == "english") {

    header_names <- data_frame(col_keys = c("var", "groups", "means", "Difference",
                                            "statistic", "df", "p.value", "signif",
                                            "CI_lower", "CI_upper", "method", "d_Cohen"),
                               names = c("Variable",
                                         if_else(is.null(group_name) == TRUE, "Group", group_name),
                                         "Mean", "Difference",
                                         "t", "df", "p", "p",
                                         "<", ">", "Test", "Cohen's d"))

    sub_header_names <- c("Variable",
                          if_else(is.null(group_name) == TRUE, "Group", group_name),
                          "Mean", "Difference",
                          "t", "df", "p", "p",
                          "CI 95%", "CI 95%", "Test", "Cohen's d")

  } else {

    header_names <- data_frame(col_keys = c("var", "groups", "means", "Difference",
                                            "statistic", "df", "p.value", "signif",
                                            "CI_lower", "CI_upper", "method", "d_Cohen"),
                               names = c("Zmienna",
                                         if_else(is.null(group_name) == TRUE, "Grupa", group_name),
                                         "Średnia", "Różnica",
                                         "t", "df", "p", "p",
                                         "<", ">", "Test", "d Cohena"))

    sub_header_names <- c("Zmienna",
                          if_else(is.null(group_name) == TRUE, "Grupa", group_name),
                          "Średnia", "Różnica",
                          "t", "df", "p", "p",
                          "CI 95%", "CI 95%", "Test", "d Cohena")

  }


  out <- flextable(out) %>%
    set_header_df(mapping = header_names, key = "col_keys") %>%
    add_header_row(values = sub_header_names) %>%
    merge_at(j = 1, i = 1:2, part = "header") %>%
    merge_at(j = 2, i = 1:2, part = "header") %>%
    merge_at(j = 3, i = 1:2, part = "header") %>%
    merge_at(j = 4, i = 1:2, part = "header") %>%
    merge_at(j = 5, i = 1:2, part = "header") %>%
    merge_at(j = 6, i = 1:2, part = "header") %>%
    merge_at(j = 7:8, i = 1:2, part = "header") %>%
    merge_at(j = 9:10, i = 1, part = "header") %>%
    merge_at(j = 11, i = 1:2, part = "header") %>%
    merge_at(j = 12, i = 1:2, part = "header") %>%
    theme_booktabs() %>%
    align(align = "center", part = "all")


  for(var in c(1:length(vars))){

    join_vert <- c(1+2*(var-1), 2+2*(var-1))

    out <- merge_at(out, i = join_vert, j = 1, part = "body") %>%
      merge_at(i = join_vert, j = 4, part = "body") %>%
      merge_at(i = join_vert, j = 5, part = "body") %>%
      merge_at(i = join_vert, j = 6, part = "body") %>%
      merge_at(i = join_vert, j = 7, part = "body") %>%
      merge_at(i = join_vert, j = 8, part = "body") %>%
      merge_at(i = join_vert, j = 9, part = "body") %>%
      merge_at(i = join_vert, j = 10, part = "body") %>%
      merge_at(i = join_vert, j = 11, part = "body") %>%
      merge_at(i = join_vert, j = 12, part = "body")

    if(h_lines == TRUE){

      out <- hline(out, i = join_vert[2], border = fp_border(width = 1))

    }
  }

  out <- out %>%
    width(j = 7, width = 0.80) %>%
    width(j = 8, width = 0.30) %>%
    align(j = 7, align = "right", part = "body") %>%
    align(j = 8, align = "left", part = "body") %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "all") %>%
    hline(border = fp_border(width = 2), i = length(vars)*2, part = "body") %>%
    fix_border_issues(part = "all") %>%
    padding(j = 7:8, padding = 1)



  return(out)

}

#' Calculate multiple linear regression and generate flextable.
#'
#' @name flex_reg
#'
#' @param data data.frame
#' @param ind_var name of column containing dependent variable
#' @param pred_vars names of columns containing indepented (predictory) variables
#' @param ind_name name of dependent variable for printing
#' @param pred_names names of independent (predictory) variables for printing
#' @param formula optionally formula to use - if there are something more complex to calculate
#' @param lang defaults to "english", may also be "polish"
#' @param width_adjust for adjusting width of columns
#' @param fontsize defaults to 9
#'
#' @import lm.beta
#' @import broom
#' @import flextable
#' @import tidyverse
#'

flex_reg <- function(data,
                     ind_var,
                     pred_vars,
                     ind_name,
                     pred_names,
                     formula = NULL,
                     lang = "english",
                     width_adjust = 1,
                     fontsize = 9){

  if(is.null(formula) == T) {

    formula <- paste(ind_var,
                     paste(pred_vars, collapse = " + "),
                     sep = " ~ ")

  }


  require(tidyverse)
  require(lm.beta)
  require(broom)
  require(flextable)

  sig_check <- function(x) {
    ifelse(x >= 0.1, " ",
           ifelse(x >= 0.05, ".",
                  ifelse(x >= 0.01, "*",
                         ifelse(x >= 0.001, "**", "***"))))
  }

  formula_grepl <- formula

  vars <- c(ind_var, pred_vars)

  n_name <- 0

  for(name in c(ind_name, pred_names)){

    n_name <- n_name + 1

    formula_grepl <- sub(pattern = vars[n_name],
                         replacement = name,
                         x = formula_grepl,
                         fixed = TRUE)

  }

  reg_test <- lm.beta(lm(formula = formula,
                         data = data))

  reg_test_1 <- tidy(reg_test) %>%
    mutate(signif = sig_check(p.value))

  n_name <- 0

  reg_test_1$term[1] <- if_else(lang == "polish", "(punkt przecięcia)", "(Intercept)")

  for(name in (pred_names)){

    n_name <- n_name + 1

    reg_test_1$term <- sub(pattern = pred_vars[n_name],
                           replacement = name,
                           x = reg_test_1$term,
                           fixed = TRUE)


  }

  reg_test_1$term <- sub(pattern = ":",
                         replacement = " x ",
                         x = reg_test_1$term,
                         fixed = TRUE)

  if(lang == "polish"){

    label_names <- data_frame(col_keys = c("term", "estimate", "std_estimate", "std.error",
                                           "statistic", "p.value", "signif"),
                              names = c("predyktor", "współczynnik", "β", "błąd standardowy",
                                        "statystyka T", "p", "p"))

  } else {


    label_names <- data_frame(col_keys = c("term", "estimate", "std_estimate", "std.error",
                                           "statistic", "p.value", "signif"),
                              names = c("predictor", "estimate", "β", "std.error",
                                        "T statistic", "p", "p"))

  }

  reg_test_1 <- reg_test_1 %>%
    mutate_at(.vars = c(2:5),
              .funs = round,
              digits = 2) %>%

    mutate_at(.vars = c(2:5),
              .funs = sprintf,
              fmt = "%0.2f",
              how = "replace")

  reg_test_1 <- reg_test_1 %>%
    mutate(p.value = round(p.value, 3)) %>%
    mutate(p.value = sprintf(p.value, fmt = "%0.3f", how = "replace"))

  reg_test_1$p.value[reg_test_1$p.value == "0.000"] <- "< 0.001"

  if(lang == "polish"){

    formula_label = paste("wzór", formula_grepl, sep = ": ")

  } else {

    formula_label = paste("formula", formula_grepl, sep = ": ")

  }

  reg_test_2 <- glance(reg_test) %>%
    mutate(signif = sig_check(p.value))

  if(lang == "polish"){

    subtop_label = data_frame(term = c("błąd standardowy reszt", as.character(formatC(reg_test_2$sigma, digits = 2, format = "f"))),
                              estimate = c("R2", as.character(formatC(reg_test_2$r.squared, digits = 2, format = "f"))),
                              std_estimate = c("R2 poprawione", as.character(formatC(reg_test_2$adj.r.squared, digits = 2, format = "f"))),
                              std.error = c("statystyka F", as.character(formatC(reg_test_2$statistic, digits = 2, format = "f"))),
                              statistic = c("df", as.character(paste(reg_test_2$df, reg_test_2$df.residual, sep = ", "))),
                              p.value = c("p", as.character(if_else(reg_test_2$p.value < 0.001, "< 0.001",
                                                                    as.character(formatC(as.numeric(reg_test_2$p.value)[1], digits = 3, format = "f"))))),
                              signif = c("p", as.character(reg_test_2$signif)))

  } else {

    subtop_label = data_frame(term = c("residual standard error", as.character(formatC(reg_test_2$sigma, digits = 2, format = "f"))),
                              estimate = c("R2", as.character(formatC(reg_test_2$r.squared, digits = 2, format = "f"))),
                              std_estimate = c("R2 adjusted", as.character(formatC(reg_test_2$adj.r.squared, digits = 2, format = "f"))),
                              std.error = c("F statistic", as.character(formatC(reg_test_2$statistic, digits = 2, format = "f"))),
                              statistic = c("df", as.character(paste(reg_test_2$df, reg_test_2$df.residual, sep = ", "))),
                              p.value = c("p", as.character(if_else(reg_test_2$p.value < 0.001, "< 0.001",
                                                                    as.character(formatC(as.numeric(reg_test_2$p.value)[1], digits = 3, format = "f"))))),
                              signif = c("p", as.character(reg_test_2$signif)))

  }

  out <- flextable(reg_test_1) %>%
    set_header_df(mapping = label_names, key = "col_keys") %>%
    add_header(values = subtop_label) %>%
    add_header_lines(values = formula_label) %>%
    theme_booktabs() %>%
    width(j = 1, width = 1.8 * width_adjust) %>%
    width(j = 2:5, width = 0.85 * width_adjust) %>%
    align(align = "center", part = "all") %>%
    merge_at(i = 2, j = 6:7, part = "header") %>%
    merge_at(i = 4, j = 6:7, part = "header") %>%
    width(j = 6, width = 0.6 * width_adjust) %>%
    width(j = 7, width = 0.25 * width_adjust) %>%
    align(j = 6, align = "right", part = "body") %>%
    align(j = 7, align = "left", part = "body") %>%
    align(i = 3, j = 6, align = "right", part = "header") %>%
    align(i = 3, j = 7, align = "left", part = "header") %>%
    padding(j = 6, padding.right = 1, part = "body") %>%
    padding(j = 7, padding.left = 1, part = "body") %>%
    padding(i = 3, j = 6, padding.right = 1, part = "header") %>%
    padding(i = 3, j = 7, padding.left = 1, part = "header") %>%
    align(j = 2:5, align = "right", part = "body") %>%
    padding(j = 2:5, padding.right = 22 * width_adjust) %>%
    font(j = 1:6, fontname = "Times New Roman", part = "all") %>%
    font(j = 7, fontname = "Arial", part = "all") %>%
    fontsize(size = fontsize, part = "all")

}

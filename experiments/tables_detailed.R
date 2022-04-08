library(tidyverse)
library(kableExtra)

idir <- sprintf("results_hpc/detailed/", dataset)
ifile.list <- list.files(idir)
results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    if("gap" %in% colnames(df)) {
        df <- df %>% select(-gap)
    }
    if("sketch" %in% colnames(df)) {
        df <- df %>% select(-sketch)
    }
    if("tracking" %in% colnames(df)) {
        df <- df %>% select(-tracking)
    }
    if(endsWith(ifile, "_95.txt")) {
        df <- df %>% mutate(confidence=0.95)
    } else if (endsWith(ifile, "_50.txt")) {
        df <- df %>% mutate(confidence=0.50)
    }
}))

# "conformal-proportional"
method.values = c("classical", "bayesian-dp-?","bootstrap", "conformal-constant", "conformal-adaptive") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap", "Fixed", "Adaptive") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
color.scale <- cbPalette[c(1,2,3,4,4,4)]
shape.scale <- c(1,2,0,9,8,7,3)

plot.w <- 5000

list.lower <- results %>%
    group_by(data) %>%
    filter(method=="classical", count < 25) %>%
    sample_n(10)

df.upper <- results %>%
    distinct(method, w, x, .keep_all=TRUE) %>%
    mutate(Method=factor(method, method.values, method.labels)) %>% select(-method) %>%
    filter(seed==1, n_bins %in% c(NA,5), w==plot.w) %>%
    group_by(data, Method, d, w, confidence, seed) %>%
    arrange(data, Method, d, w, confidence, seed, desc(count), x) %>%    
    slice_max(count, with_ties=FALSE, n=10) %>%
    ungroup() %>%
    select(-mean, -median, -mode, -seen, -`method-unique`, -posterior, -n_bins, -n, -seed, -n_track, -d, -w, -confidence) %>%
    group_by(data,x,count,upper) %>%
    mutate(lower.max=max(lower[lower<=count])) %>%
    mutate(lower.max.2=sort(lower[lower<=count],decreasing=T)[2]) %>%
    ungroup() %>%
    mutate(lower.str = ifelse(lower>count, sprintf("\\textcolor{red}{%s}",lower), sprintf("\\textcolor{darkgreen}{%s}",lower))) %>%
    mutate(lower = ifelse((lower==lower.max)*(lower<=count)*(lower>lower.max.2), sprintf("\\textbf{%s}",lower.str), lower.str)) %>%    
    select(-lower.str, -lower.max, -lower.max.2) %>%    
    pivot_wider(names_from=c("Method"), values_from=c("lower")) %>%
    select(data, everything())
df.lower <- results %>%
    distinct(method, w, x, .keep_all=TRUE) %>%
    mutate(Method=factor(method, method.values, method.labels)) %>% select(-method) %>%
    filter(seed==1, n_bins %in% c(NA,5), w==plot.w) %>%
    group_by(data, Method, d, w, confidence, seed) %>%
    arrange(data, Method, d, w, confidence, seed, desc(count), x) %>%
#    slice_min(count, with_ties=FALSE, n=10) %>%
    filter(x %in% list.lower$x) %>%
    ungroup() %>%
    select(-mean, -median, -mode, -seen, -`method-unique`, -posterior, -n_bins, -n, -seed, -n_track, -d, -w, -confidence) %>%
    group_by(data,x,count,upper) %>%
    mutate(lower.max=max(lower[lower<=count])) %>%
    mutate(lower.max.2=sort(lower[lower<=count],decreasing=T)[2]) %>%
    ungroup() %>%
    mutate(lower.str = ifelse(lower>count, sprintf("\\textcolor{red}{%s}",lower), sprintf("\\textcolor{darkgreen}{%s}",lower))) %>%
    mutate(lower = ifelse((lower==lower.max)*(lower<=count)*(lower>lower.max.2), sprintf("\\textbf{%s}",lower.str), lower.str)) %>%
    select(-lower.str, -lower.max, -lower.max.2) %>%
    pivot_wider(names_from=c("Method"), values_from=c("lower")) %>%
    select(data, everything())


df <- rbind(df.upper, df.lower) %>%
    arrange(data, desc(count), x)

tb1 <- df %>%
    select(-data) %>%
    kbl("latex", booktabs=TRUE, longtable = FALSE, escape = FALSE,
        col.names = c("Data", "Frequency", "Upper bound", method.labels)) %>%
    pack_rows(index = table(df$data)) %>%
    add_header_above(c(" " = 6, "Conformal" = 2)) %>%
    add_header_above(c(" " = 3, "Lower bound" = length(method.labels)))
tb1


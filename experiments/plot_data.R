library(tidyverse)

dataset = "words"

# "conformal-proportional"
method.values = c("classical", "bayesian-dp-?","bootstrap", "conformal-constant", "conformal-adaptive") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap", "Conformal (fixed)", "Conformal (adaptive)") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
color.scale <- cbPalette[c(1,2,3,4,4,4)]
shape.scale <- c(1,2,0,9,8,7,3)

idir <- sprintf("results_hpc/%s/marginal/", dataset)
ifile.list <- list.files(idir)
results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    if("gap" %in% colnames(df)) {
        df <- df %>% select(-gap)
    }
    if(endsWith(ifile, "_95.txt")) {
        df <- df %>% mutate(confidence=0.95)
    } else if (endsWith(ifile, "_50.txt")) {
        df <- df %>% mutate(confidence=0.50)
    }
}))

# Default arguments
plot_n = 1000000
model_name = dataset
plot_method_unique=TRUE
plot_unique=TRUE
plot_n_bins=5

make_plot_marginal = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    data.list <- dataset
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==0.95,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=w, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
                                        #    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
        facet_wrap(.~key, scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') +
        xlab("Hash width") +
        ylab("") +
        theme_bw() +
        theme(legend.position = "right",
              text = element_text(size = 10),
              axis.text=element_text(size=10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              strip.text.x = element_text(size = 10),
              strip.text.y = element_text(size = 10),
              axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust=1)
              )
    suffix_seen = ifelse(plot_seen==T, "", "_unseen")
    filename <- sprintf("figures/%s_%s_w_mu%s_u%s_bins%s%s.pdf", dataset, plot_sketch, plot_method_unique, plot_unique, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=2, width=5.5, units="in")
    #pp
}
make_plot_marginal(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5, plot_sketch="cms-cu")
make_plot_marginal(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5, plot_sketch="cms")
make_plot_marginal(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5, plot_sketch="cms-cu", plot_seen=FALSE)
make_plot_marginal(dataset, plot_method_unique=FALSE, plot_unique=TRUE, plot_n_bins=5, plot_sketch="cms-cu", plot_seen=FALSE)
make_plot_marginal(dataset, plot_method_unique=FALSE, plot_unique=TRUE, plot_n_bins=5, plot_sketch="cms-cu")

#######################
## Compares sketches ##
#######################

make_plot_marginal_cms = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    data.list <- dataset
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`, sketch %in% c("cms", "cms-cu"),
               data %in% data.list, confidence==0.95,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        gather(Coverage, Length, key="key", value="value") %>%
        mutate(sigma = as.numeric(sigma),
               Sketch = factor(sketch, c("cms-cu", "cms"), c("CMS-CU", "CMS"))) %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, Sketch, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=w, y=value, color=Method, shape=Method, linetype=Sketch)) +
        #geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
        facet_grid(key~Method, scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') +
        xlab("Hash width") +
        ylab("") +
        guides(shape=FALSE, color=FALSE) +
        theme_bw() +
        theme(legend.position = "bottom",
              text = element_text(size = 10),
              axis.text=element_text(size=10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              strip.text.x = element_text(size = 10),
              strip.text.y = element_text(size = 10),
              axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust=1)
              )
    filename <- sprintf("figures/%s_sketch_w_mu%s_u%s_bins%s.pdf", dataset, plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=3.5, width=8, units="in")    
}
make_plot_marginal_cms(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5)


######################
## Estimation plots ##
######################

make_plot_estim = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    data.list <- dataset
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`, sketch=="cms-cu",
               data %in% data.list, (confidence==0.5) | (!startsWith(results$method, "conformal-")),
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        mutate(Error = `error-median`) %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, `include-seen`) %>%
        summarise(Error=mean(Error))
    pp <- df %>%
        filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=w, y=Error, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') +
        xlab("Hash width") +
        ylab("MAD") +
        theme_bw() +
        theme(legend.position = "right",
              text = element_text(size = 10),
              axis.text=element_text(size=10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              strip.text.x = element_text(size = 10),
              strip.text.y = element_text(size = 10)
              )
    filename <- sprintf("figures/%s_error_w_mu%s_u%s_bins%s.pdf", dataset, plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=2, width=5.5, units="in")    
}
make_plot_estim(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5)



make_plot_cms_estim = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    data.list <- dataset
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`, sketch %in% c("cms", "cms-cu"),
               data %in% data.list, (confidence==0.5) | (!startsWith(results$method, "conformal-")),
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        mutate(Error = `error-median`, Sketch = factor(sketch, c("cms-cu", "cms"), c("CMS-CU", "CMS"))) %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, `include-seen`, Sketch) %>%
        summarise(Error=mean(Error))
    pp <- df %>%
        filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=w, y=Error, color=Method, shape=Method, linetype=Sketch)) +
        #geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        facet_grid(.~Method, scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') +
        guides(shape=FALSE, color=FALSE) +
        xlab("Hash width") +
        ylab("MAD") +
        theme_bw() +
        theme(legend.position = "right",
              text = element_text(size = 10),
              axis.text=element_text(size=10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              strip.text.x = element_text(size = 10),
              strip.text.y = element_text(size = 10),
              axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust=1)
              )
    filename <- sprintf("figures/%s_cms_error_w_mu%s_u%s_bins%s.pdf", dataset, plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=2, width=5.5, units="in")    
}
#make_plot_cms_estim(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=5)


#######################
## Conditional plots ##
#######################

idir <- sprintf("results_hpc/%s/conditional/", dataset)
ifile.list <- list.files(idir)
results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    if("gap" %in% colnames(df)) {
        df <- df %>% select(-gap)
    }
    if(endsWith(ifile, "_95.txt")) {
        df <- df %>% mutate(confidence=0.95)
    } else if (endsWith(ifile, "_50.txt")) {
        df <- df %>% mutate(confidence=0.50)
    }
}))

make_plot_cond = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    data.list <- dataset
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data %in% data.list,
               n_bins %in% c(NA,plot_n_bins), sketch=="cms-cu", confidence==0.95) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, `count-bin`, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(Freq=factor(`count-bin`, bin.values, bin.labels)) %>%
        ggplot(aes(x=w, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), linetype=2, alpha=0) +    
        geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1, alpha=0.5) +
        facet_grid(key~Freq, scales="free", labeller = labeller(Freq = label_both, key = label_value)) +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') +
        xlab("Hash width") +
        ylab("") +
        theme_bw() +
        theme(legend.position = "bottom",
              text = element_text(size = 10),
              axis.text=element_text(size=10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              strip.text.x = element_text(size = 10),
              strip.text.y = element_text(size = 10),
              axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust=1)
              )
    filename <- sprintf("figures/%s_conditional_w_mu%s_u%s_bins%s.pdf", dataset, plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=4.5, width=8, units="in")
}
#make_plot_cond(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1)
make_plot_cond(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5)


##########################
## Frequency histograms ##
##########################

df.covid <- read_csv("../../notebooks/freq_covid.txt", col_names="Frequency") %>% mutate(Data="SARS-CoV-2")
df.words <- read_csv("../../notebooks/freq_words.txt", col_names="Frequency") %>% mutate(Data="Novels")

df <- rbind(df.covid, df.words) 

pp <- df %>%
    mutate(Data=factor(Data, c("SARS-CoV-2", "Novels"), c("SARS-CoV-2", "Literature"))) %>%
    ggplot(aes(x=Frequency))+
    geom_histogram()+
    facet_grid(.~Data, scales="free") +
        ylab("") +
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10') +
    theme_bw() +
    theme(legend.position = "bottom",
          text = element_text(size = 10),
          axis.text=element_text(size=10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 10),
          axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust=1)
          )
pp %>% ggsave(file="figures/data_frequencies.pdf", height=2.5, width=6, units="in")
    

## ########################
## ## Conditional Tables ##
## ########################
## library(kableExtra)


## make_table_cond = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
##     data.list <- dataset
##     df <- results %>%
##         mutate(Coverage = coverage, Length = length) %>%
##         filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data %in% data.list,
##                n_bins %in% c(NA,plot_n_bins)) %>%
##         separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
##         mutate(sigma = as.numeric(sigma)) %>%
##         group_by(sigma, d, w, n, `count-bin`, method, posterior, n_bins, n_track, `include-seen`) %>%
##         summarise(Coverage=mean(Coverage), Length=mean(Length)) %>%
##         ungroup() %>%
##         filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
##         mutate(Method=factor(method, method.values, method.labels)) %>%
##         mutate(Frequency=factor(`count-bin`, bin.values, bin.labels)) %>%
##         select(w, method, Coverage, Length,`count-bin`)
        
##     print(df, n=1000)
##     ## pp <- df %>%
##     ##     filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
##     ##     mutate(Method=factor(method, method.values, method.labels)) %>%
##     ##     mutate(Frequency=factor(`count-bin`, bin.values, bin.labels)) %>%
##     ##     ggplot(aes(x=w, y=value, color=Method, shape=Method)) +
##     ##     geom_point(alpha=0.5) +
##     ##     geom_line(alpha=0.5) +
##     ##     geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
##     ##     geom_hline(data=df.dummy2, aes(yintercept=value), linetype=2, alpha=0) +    
##     ##                                     #    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
##     ##     facet_wrap(key~Frequency, scales="free", nrow=2) +
##     ##     scale_color_manual(values=color.scale) +
##     ##     scale_shape_manual(values=shape.scale) +
##     ##     xlab("Sigma") +
##     ##     ylab("") +
##     ##     theme_bw() +
##     ##     theme(legend.position = "bottom")
##     filename <- sprintf("tables/%s_conditional_w_mu%s_u%s_bins%s.tex", dataset, plot_method_unique, plot_unique, plot_n_bins)
##     #pp %>% ggsave(file=filename, height=4.5, width=8, units="in")
## }
## make_table_cond(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1)


## data.list <- dataset
## df <- results %>%
##     mutate(Coverage = coverage, Length = length) %>%
##     filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data %in% data.list,
##            n_bins %in% c(NA,plot_n_bins)) %>%
##     separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
##     mutate(sigma = as.numeric(sigma)) %>%
##     group_by(sigma, d, w, n, `count-bin`, method, posterior, n_bins, n_track, `include-seen`) %>%
##     summarise(Coverage=mean(Coverage), Length=mean(Length)) %>%
##     ungroup() %>%
##     filter(posterior=="mcmc", n==plot_n, method %in% method.values) %>%
##     mutate(Method=factor(method, method.values, method.labels)) %>%
##     mutate(Frequency=factor(`count-bin`, bin.values, bin.labels)) %>%
##     select(w, Method, Coverage, Length,`count-bin`) %>%
##     arrange(w, `count-bin`, Method)

## df %>% pivot_wider(names_from=Method, values_from=c("Length", "Coverage"))

library(tidyverse)
library(gridExtra)

dataset <- "zipf"
dataset <- "pyp"

if(dataset=="zipf") {
    data.list <- paste("zipf",c("1.05","1.1","1.2", "1.3","1.4", "1.5","1.6","1.7","1.8","1.9","2.0","2.5","3.0"),sep="-")
    plot_n = 1000000
    plot_w = 1000
} else {
    data.list <- paste("pyp-5000",c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","0.95"),sep="-")
    plot_n = 100000
    plot_w = 1000
}

results <- do.call("rbind", lapply(data.list, function(data.name) {
    idir <- sprintf("results_hpc/%s/conditional", data.name)
    ifile.list <- list.files(idir)
    res <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
        if(endsWith(ifile, "_95.txt")) {
            df <- df %>% mutate(confidence=0.95)
        } else if (endsWith(ifile, "_50.txt")) {
            df <- df %>% mutate(confidence=0.50)
        }
    }))
}))


# "conformal-proportional"
method.values = c("classical", "bayesian-dp-?","bootstrap", "conformal-constant", "conformal-adaptive") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap", "Conformal (fixed)", "Conformal (adaptive)") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
color.scale <- cbPalette[c(1,2,3,4,4,4)]
shape.scale <- c(1,2,0,9,8,7,3)

bin.values <- c(0,1,2,3,4)
bin.labels <- c("(0,20]%", "(20,40]%", "(40,60]%", "(60,80]%", "(80,100]%")

# Default arguments
model_name = dataset
plot_method_unique=TRUE
plot_unique=TRUE
plot_sketch="cms-cu"
plot_n_bins=5

#########################
## Plots vs data model ##
#########################

# Compare different lower bounds

plot_param_lb = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu") {
    if(dataset=="zipf") {
        xlabel <- TeX('Parameter $a$')
        xlims = c(1,1.75)
        xbreaks = c(1,1.25,1.5,1.75)
        param.split = c("tmp1", "sigma")
        scale_y = scale_y_log10()
        bin.values <- c(0,1,2,3)
        bin.labels <- c("Low frequencies", "Low frequencies", "High frequencies", "High frequencies")
    } else if(dataset=="pyp") {
        xlabel <- TeX("Parameter $\\sigma$")
        xlims = c(0,0.4)
        xbreaks = c(0,0.2,0.4)
        param.split = c("tmp1", "tmp2", "sigma")
        scale_y = scale_y_continuous()
        #bin.values <- c(0,1,2,3)
        #bin.labels <- c("Low frequencies", "Low frequencies", "High frequencies", "High frequencies")
        bin.values <- c(0,1,2,3,4)
        bin.labels <- c("(0,20]%", "(20,40]%", "(40,60]%", "(60,80]%", "(80,100]%")
    } else {
        xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`, sketch==plot_sketch,
               data %in% data.list, confidence==0.95,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        mutate(Freq=factor(`count-bin`, bin.values, bin.labels)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, Freq, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value)) %>%
        ungroup()
    pp <- df %>%
        filter(posterior=="mcmc", w==plot_w, n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        filter(!is.na(Freq)) %>%
        ggplot(aes(x=sigma, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
                                        #    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
        facet_grid(key~Freq, scales="free", labeller = labeller(Freq = label_value, key = label_value)) +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y +
        scale_x_continuous(lim=xlims, breaks=xbreaks) +
        xlab(xlabel) +
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
    filename <- sprintf("figures/%s_%s_conditional_param_mu%s_u%s_bins%s.pdf", dataset, plot_sketch, plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=3, width=8, units="in")
}
#plot_param_lb(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1)
plot_param_lb(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5, plot_sketch="cms-cu")
#plot_param_lb(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=5, plot_sketch="cms")




if(TRUE) {
    idir <- "results_hpc/results_conditional/"
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
        if(endsWith(ifile, "_95.txt")) {
            df <- df %>% mutate(confidence=0.95)
        } else if (endsWith(ifile, "_50.txt")) {
            df <- df %>% mutate(confidence=0.50)
        }
    }))
    if(FALSE) {
        results %>% write_delim("summary/summary.txt")
    }
} else {
    results <- read_delim("summary/summary.txt", delim=" ")
}

method.values = c("classical", "bootstrap", "conformal-constant", "conformal-adaptive", "conformal-bayesian-dp", "bayesian-dp-?")
method.labels = c("Classical", "Bootstrap", "Conformal (constant)", "Conformal (adaptive)", "Conformal (Bayesian)", "Bayesian (DP)")
bin.values <- c(0,1,2,3,4)
bin.labels <- c("(0,20]%", "(20,40]%", "(40,60]%", "(60,80]%", "(80,100]%")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
color.scale <- cbPalette[c(1,2,3,3,3,6,6)]
shape.scale <- c(1,2,4,8,7,6,3)

make_plot_param = function(model_name, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    if(model_name=="pyp") {
        data.list <- c("pyp-100-0.0", paste("pyp-100-", seq(0.1,0.9,0.1), sep=""))
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data %in% data.list,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, `count-bin`, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    df.dummy2 <- tibble(key="Coverage", value=0.5)
    pp <- df %>%
        filter(posterior=="mcmc", w==333, n==100000, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(Frequency=factor(`count-bin`, bin.values, bin.labels)) %>%
        ggplot(aes(x=sigma, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), linetype=2, alpha=0) +    
                                        #    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
        facet_wrap(key~Frequency, scales="free", nrow=2) +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        xlab("Sigma") +
        ylab("") +
        theme_bw() +
        theme(legend.position = "bottom")
    filename <- sprintf("figures/pyp_conditional_sigma_mu%s_u%s_bins%s.pdf", plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=4.5, width=8, units="in")
}
make_plot_param("pyp", plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1)




data.name <- "dp-100"

df <- results %>%
    mutate(Coverage = coverage, Length = length) %>%
    filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data==data.name) %>%
    separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
    mutate(sigma = as.numeric(sigma)) %>%
    gather(Coverage, Length, key="key", value="value") %>%
    group_by(sigma, d, w, n, `count-bin`, method, posterior, n_bins, n_track, key, `include-seen`) %>%
    summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
df.dummy2 <- tibble(key="Coverage", value=0)
pp <- df %>%
    filter(posterior=="mcmc", w==333, method %in% method.values) %>%
    mutate(Method=factor(method, method.values, method.labels)) %>%
    mutate(Frequency=factor(`count-bin`, bin.values, bin.labels)) %>%
    ggplot(aes(x=n, y=value, color=Method, shape=Method)) +
    geom_point(alpha=0.5) +
    geom_line(alpha=0.5) +
    geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
    geom_hline(data=df.dummy2, aes(yintercept=value), linetype=2, alpha=0) +    
#    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
    facet_wrap(key~Frequency, scales="free", nrow=2) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    xlab("Sigma") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "bottom")
filename <- sprintf("figures/pyp_conditional_%s_n_mu%s_u%s.pdf", data.name, plot_method_unique, plot_unique)
pp %>% ggsave(file=filename, height=4.5, width=8, units="in")


if(FALSE) {
    data.list <- c("dp-1", "dp-100", "pyp-1-0.5", "pyp-100-0.5")

    df <- results %>%
        filter(`include-seen`, data %in% data.list) %>%
        gather(coverage, length, key="key", value="value") %>%
        group_by(data, d, w, n, `count-bin`, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))

    df.dummy2 <- tibble(key="Coverage", value=0)
    pp <- df %>%
        filter(key=="coverage", posterior=="mcmc", w==333) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(Frequency=factor(`count-bin`, bin.values, bin.labels)) %>%
        ggplot(aes(x=n, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), linetype=2, alpha=0) +    
        geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
        facet_grid(Frequency~data, scales="free", ncol=length(data.list)) +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        xlab("Sigma") +
        ylab("") +
        theme_bw() +
        theme(legend.position = "bottom")
    pp

}

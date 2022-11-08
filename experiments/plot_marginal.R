library(tidyverse)
library(latex2exp)

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
    idir <- sprintf("results_hpc/%s/marginal", data.name)
    ifile.list <- list.files(idir)
    res <- do.call("rbind", lapply(ifile.list, function(ifile) {
        if(startsWith(ifile, "exp1")) {
            df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
            if(endsWith(ifile, "_95.txt")) {
                df <- df %>% mutate(confidence=0.95)
            } else if (endsWith(ifile, "_50.txt")) {
                df <- df %>% mutate(confidence=0.50)
            }
        } else {
            df <- tibble()
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

# Default arguments
model_name = dataset
plot_method_unique=FALSE
plot_unique=FALSE
plot_n_bins=5

#########################
## Plots vs data model ##
#########################

# Compare different lower bounds

plot_param_lb = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    if(dataset=="zipf") {
        xlabel <- TeX('Parameter $a$')
        param.split = c("tmp1", "sigma")
        xlims = c(1,1.75)
        xbreaks = c(1,1.25,1.5,1.75)
        title <- "(a) Zipf"
    } else if(dataset=="pyp") {
        xlabel <- TeX("Parameter $\\sigma$")
        param.split = c("tmp1", "tmp2", "sigma")
        xlims = c(0,1)
        xbreaks = c(0,0.5,1)
        title <- "(b) PYP"
    } else {
        xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==0.95,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(posterior=="mcmc", w==plot_w, n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=sigma, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
        ##    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
        facet_wrap(.~key, scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
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
              strip.text.y = element_text(size = 10)
              )
    suffix_seen = ifelse(plot_seen==T, "", "_unseen")
    filename <- sprintf("figures/%s_%s_param_mu%s_u%s_bins%s%s.pdf", dataset, plot_sketch, plot_method_unique, plot_unique, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=2, width=6, units="in")
}
#plot_param_lb(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1)
plot_param_lb(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5, plot_sketch = "cms-cu")
plot_param_lb(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5, plot_sketch = "cms")
plot_param_lb(dataset, plot_method_unique=FALSE, plot_unique=TRUE, plot_n_bins=5, plot_sketch = "cms-cu", plot_seen=FALSE)
plot_param_lb(dataset, plot_method_unique=FALSE, plot_unique=TRUE, plot_n_bins=5, plot_sketch = "cms-cu")


# Compares sketches
plot_param_cms = function(dataset, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    if(dataset=="zipf") {
        xlabel <- TeX('Parameter $a$')
        param.split = c("tmp1", "sigma")
        xlims = c(1,1.75)
        xbreaks = c(1,1.25,1.5,1.75)
        title <- "(a) Zipf"
    } else if(dataset=="pyp") {
        xlabel <- "Sigma"
        param.split = c("tmp1", "tmp2", "sigma")
        xlims = c(0,1)
        xbreaks = c(0,0.5,1)
        title <- "(b) PYP"
    } else {
        xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`==T,  sketch %in% c("cms", "cms-cu"),
               data %in% data.list, confidence==0.95,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma),
               Sketch = factor(sketch, c("cms-cu", "cms"), c("CMS-CU", "CMS"))) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, Sketch, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(posterior=="mcmc", w==plot_w, n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=sigma, y=value, color=Method, shape=Method, linetype=Sketch)) +
        #geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
        facet_grid(key~Method, scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
        scale_x_continuous(lim=xlims, breaks=xbreaks) +
        xlab(xlabel) +
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
    filename <- sprintf("figures/%s_sketch_mu%s_u%s_bins%s.pdf", dataset, plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=3.5, width=8, units="in")    
}
plot_param_cms(dataset, plot_method_unique=FALSE, plot_unique=FALSE, plot_n_bins=5)





make_plot_param = function(model_name, plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=1) {
    if(model_name=="pyp") {
        data.list <- c("pyp-100-0.0", paste("pyp-100-", seq(0.1,0.9,0.1), sep=""))
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique),
               `include-seen`,
               data %in% data.list, confidence==0.95,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, posterior, n_bins, n_track, key, `include-seen`) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(posterior=="mcmc", w==333, n==100000, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        ggplot(aes(x=sigma, y=value, color=Method, shape=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
                                        #    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
        facet_wrap(.~key, scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        xlab("Sigma") +
        ylab("") +
        theme_bw() +
        theme(legend.position = "right")
    filename <- sprintf("figures/pyp_sigma_mu%s_u%s_bins%s.pdf", plot_method_unique, plot_unique, plot_n_bins)
    pp %>% ggsave(file=filename, height=2.5, width=7, units="in")    
}
make_plot_param("pyp", plot_method_unique=TRUE, plot_unique=TRUE, plot_n_bins=5)





data.list <- c("dp-1", "dp-100", "pyp-1-0.5", "pyp-100-0.5")
df <- results %>%
    mutate(Coverage = coverage, Length = length, `Length (log)`=log(Length)) %>%
    filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data %in% data.list, confidence==0.95) %>%
    gather(Coverage, `Length`, key="key", value="value") %>%
    group_by(data, d, w, n, method, posterior, n_bins, n_track, key, `include-seen`) %>%
    summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
pp <- df %>%
    filter(posterior=="mcmc", w==333, method %in% method.values) %>%
    mutate(Method=factor(method, method.values, method.labels)) %>%
    ggplot(aes(x=n, y=value, color=Method, shape=Method)) +
    geom_point(alpha=0.75) +
    geom_line(alpha=0.75) +
#    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
    geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
    facet_grid(key~data, scales="free") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Sample size") +
    ylab("Performance") +
    theme_bw() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
filename <- sprintf("figures/pyp_n_mu%s_u%s.pdf", plot_method_unique, plot_unique)
pp %>% ggsave(file=filename, height=4, width=8, units="in")


data.list <- c("zipf-1.1", "zipf-1.25", "zipf-1.5", "zipf-1.75", "zipf-2.0", "zipf-2.5", "zipf-3.0")
df <- results %>%
    mutate(Coverage = coverage, Length = length) %>%
    filter(unique==plot_unique, `method-unique`%in% c(NA,plot_method_unique), `include-seen`, data %in% data.list, confidence==0.95) %>%
    separate(data, c("tmp1", "sigma"), sep="-") %>%
    mutate(sigma = as.numeric(sigma)) %>%
    gather(Coverage, Length, key="key", value="value") %>%
    group_by(sigma, d, w, n, method, posterior, n_bins, n_track, key, `include-seen`) %>%
    summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
pp <- df %>%
    filter(posterior=="mcmc", w==333, n==100000, method %in% method.values) %>%
    mutate(Method=factor(method, method.values, method.labels)) %>%
    ggplot(aes(x=sigma, y=value, color=Method, shape=Method)) +
    geom_point(alpha=0.5) +
    geom_line(alpha=0.5) +
    geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
    geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
#    geom_errorbar(aes(ymin=value-se, ymax=value+se)) +
    facet_wrap(.~key, scales="free") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_y_log10() +
    xlab("Sigma") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "right")
filename <- sprintf("figures/zipf_alpha_mu%s_u%s.pdf", plot_method_unique, plot_unique)
pp %>% ggsave(file=filename, height=2.5, width=7, units="in")


##############################
## Plot estimation accuracy ##
##############################

# "conformal-proportional"
method.values = c("bayesian-dp-?", "conformal-bayesian-dp", "conformal-constant", "conformal-adaptive")
method.labels = c("Bayesian (DP)", "Conformal (Bayesian)", "Conformal (constant)", "Conformal (adaptive)")

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df.dummy <- tibble(key="Coverage", value=0.95)
color.scale <- cbPalette[c(2,3,3,3,3)]
shape.scale <- c(2,4,8,7,6)

data.list <- c("pyp-100-0.0", paste("pyp-100-", seq(0.1,0.9,0.1), sep=""))
df <- results %>%
    filter(!(startsWith(method, "conformal")*(confidence==0.95))) %>%
    mutate(Error = `error-median`) %>%
    filter(`include-seen`, data %in% data.list) %>%
    separate(data, c("tmp1", "tmp2", "sigma"), sep="-") %>%
    mutate(sigma = as.numeric(sigma)) %>%
    group_by(sigma, d, w, n, method, posterior, n_bins, n_track, `include-seen`) %>%
    summarise(se=2*sd(Error)/sqrt(n()), Error=mean(Error))
pp <- df %>%
    filter(posterior=="mcmc", w==333, n==100000, method %in% method.values) %>%
    mutate(Method=factor(method, method.values, method.labels)) %>%
    ggplot(aes(x=sigma, y=Error, color=Method, shape=Method)) +
    geom_point(alpha=0.5) +
    geom_line(alpha=0.5) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    ylab("Absolute error") +
    xlab("Sigma") +
    theme_bw() +
    theme(legend.position = "right")
pp %>% ggsave(file="figures/pyp_sigma_estimation.pdf", height=2, width=4.5, units="in")

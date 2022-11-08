options(width=175)

library(tidyverse)
library(latex2exp)

dataset <- "zipf"

if(dataset=="zipf") {
    data.list <- paste("zipf",c("1.05","1.1","1.2", "1.3","1.4", "1.5","1.6","1.7","1.8","1.9","2.0","2.5","3.0"),sep="-")
    plot_n = 1000000
    plot_w = 1000
    plot_confidence = 0.95
}

results <- do.call("rbind", lapply(data.list, function(data.name) {
    idir <- sprintf("results_hpc/%s/marginal", data.name)
    ifile.list <- list.files(idir)
    res <- do.call("rbind", lapply(ifile.list, function(ifile) {
        if(startsWith(ifile, "exp4")) {
            df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
        } else {
            df <- tibble()
        }
    }))
}))




# "conformal-proportional"
method.values = c("classical", "bayesian-dp-?","bootstrap","conformal-constant", "theory") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap","Conformal (fixed)", "Lower bound (asymptotic)") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#008a65", "#007756", "#006348", "#004f39", "#003b2b", "#00281d", "#00140e")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.8)
color.scale <- cbPalette[c(1,4,4,4)]
shape.scale <- c(1,9,NA,7,3)
linetype.scale <- c(1,1,3)

# Default arguments
model_name = dataset
plot_n_bins=1
plot_sketch="cms-cu"
plot_seen=TRUE

plot_param_lb_M = function(dataset, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    if(dataset=="zipf") {
        xlabel <- TeX('Distribution shift')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(breaks=c(0,0.5,1))
        title <- "(a) Zipf"        
    } else if(dataset=="pyp") {
        xlabel <- TeX('Distribution shift')
        param.split = c("tmp1", "tmp2", "sigma")
        x_scale = scale_x_continuous(breaks=c(0,0.5,1))
        title <- "(b) PYP"
    } else if(dataset=="dp") {
        xlabel <- TeX('Distribution shift')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(breaks=c(0,0.5,1))
        title <- "(b) DP"
    } else {
        xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        mutate(Intervals = factor(two_sided, c(FALSE,TRUE), c("One-sided", "Two-sided"))) %>%
        filter(((`method-unique`==100)*(unique==TRUE)==1) + ((`method-unique`==1)*(unique==FALSE)==1)>0,
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==plot_confidence,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, `method-unique`, posterior, n_bins, n_track, key, `include-seen`, Intervals, shift) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    M.vals <- unique(na.omit(df$`method-unique`))
    pp <- df %>%
#        rbind(df.classical) %>%
        filter(key=="Coverage") %>%
        filter(sigma %in% c(1.1,1.2,1.3)) %>%
        filter(sigma<1.8) %>%
        mutate(`Tail parameter` = sigma) %>% 
        filter(posterior=="mcmc", method=="conformal-constant", w==plot_w, n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(Calibration=factor(`method-unique`, c(1,100), c("Marginal", "Unique"))) %>%
        ggplot(aes(x=shift, y=value, color=Calibration, shape=Calibration)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
        geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
        facet_grid(key~`Tail parameter`, scales="free_y", labeller = labeller(.rows=label_value, .cols=label_both)) +
        x_scale +
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
    ## suffix_seen = ifelse(plot_seen==T, "", "_unseen")
    filename <- sprintf("figures/%s_%s_bins%s%s_unique_M_shift.pdf", dataset, plot_sketch, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=1.85, width=6, units="in")
    pp
}

if(dataset=="zipf") {
    plot_param_lb_M(dataset, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE)
}




dataset = "words"

if(dataset=="words") {
    data.list <- c("words")
    plot_n = 1000000
    plot_w = c(5000)
    plot_confidence = 0.95
} else if(dataset=="covid") {
    data.list <- c("covid")
    plot_n = 1000000
    plot_w = c(5000)
    plot_confidence = 0.95
}

idir <- sprintf("results_hpc/%s/marginal/", dataset)
ifile.list <- list.files(idir)
results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    if("gap" %in% colnames(df)) {
        df <- df %>% select(-gap)
    }
    if(startsWith(ifile, "exp4")) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    } else {
        df <- tibble()
    }
}))


# "conformal-proportional"
method.values = c("classical", "bayesian-dp-?","bootstrap","conformal-constant", "theory") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap","Conformal (fixed)", "Lower bound (asymptotic)") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#008a65", "#007756", "#006348", "#004f39", "#003b2b", "#00281d", "#00140e")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.8)
color.scale <- cbPalette[c(1,2,3,4,4,4)]
shape.scale <- c(1,2,0,9,NA,7,3)
linetype.scale <- c(1,1,1,1,3)

# Default arguments
model_name = dataset
plot_n_bins=1
plot_sketch="cms-cu"
plot_seen=TRUE

plot_param_lb_M = function(dataset, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    if(dataset=="zipf") {
        xlabel <- TeX('Distribution shift')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(breaks=c(0,0.5,1))
        title <- "(a) Zipf"        
    } else if(dataset=="words") {
        xlabel <- TeX('Distribution shift')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(breaks=c(0,0.5,1))
        title <- "(a) Zipf"        
    } else if(dataset=="covid") {
        xlabel <- TeX('Distribution shift')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(breaks=c(0,0.5,1))
        title <- "(a) Zipf"        
    } else { 
       xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        mutate(Intervals = factor(two_sided, c(FALSE,TRUE), c("One-sided", "Two-sided"))) %>%
        filter(((`method-unique`==100)*(unique==TRUE)==1) + ((`method-unique`==1)*(unique==FALSE)==1)>0,
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==plot_confidence,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(d, w, n, method, `method-unique`, posterior, n_bins, n_track, key, `include-seen`, Intervals, shift) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    M.vals <- unique(na.omit(df$`method-unique`))
    pp <- df %>%
        filter(posterior=="mcmc", method=="conformal-constant", w==plot_w, n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(Calibration=factor(`method-unique`, c(1,100), c("Marginal", "Unique"))) %>%
        ggplot(aes(x=shift, y=value, color=Calibration, shape=Calibration)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
        geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
        facet_wrap(.~key, scales="free") +
        ## scale_color_manual(values=color.scale) +
        ## scale_shape_manual(values=shape.scale) +
        ## scale_linetype_manual(values=linetype.scale) +
#        scale_y_continuous(trans='log10') +
        x_scale +
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
    filename <- sprintf("figures/%s_%s_bins%s%s_unique_M_shift.pdf", dataset, plot_sketch, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=2.25, width=6, units="in")
    pp
}

plot_param_lb_M(dataset, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE)

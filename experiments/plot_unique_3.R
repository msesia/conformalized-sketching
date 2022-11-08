options(width=175)

library(tidyverse)
library(latex2exp)

dataset <- "zipf"

if(dataset=="zipf") {
    data.list <- paste("zipf",c("1.05","1.1","1.2", "1.3","1.4", "1.5","1.6","1.7","1.8","1.9","2.0","2.5","3.0"),sep="-")
    plot_n = 1000000
    plot_w = 1000
    plot_confidence = 0.95
} else if (dataset=="dp") {
    data.list <- paste("dp",c("1","10","100","1000"),sep="-")
    plot_n = 10000
    plot_sigma = 1000
    plot_confidence = 0.95
} else {
    data.list <- paste("pyp-5000",c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","0.95"),sep="-")
    plot_n = 100000
    plot_w = 1000
    plot_confidence = 0.95
}

results <- do.call("rbind", lapply(data.list, function(data.name) {
    idir <- sprintf("results_hpc/%s/marginal", data.name)
    ifile.list <- list.files(idir)
    res <- do.call("rbind", lapply(ifile.list, function(ifile) {
        if(startsWith(ifile, "exp3")) {
            df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
        } else {
            df <- tibble()
        }
    }))
}))



# "conformal-proportional"
method.values = c("conformal-constant") #, "conformal-bayesian-dp")
method.labels = c("Conformal (fixed)") #, "Conformal (Bayesian)")
cbPalette <- c("#009E73", "#008a65", "#007756", "#006348", "#004f39", "#003b2b", "#00281d", "#00140e")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
color.scale <- cbPalette[c(1,4,7)]
shape.scale <- c(65,66,67,68,69)

# Default arguments
model_name = dataset
plot_unique=FALSE
plot_sketch="cms-cu"
plot_n_bins=1
plot_seen=TRUE

#########################
## Plots vs data model ##
#########################

## Compare different lower bounds
plot_param_lb = function(dataset, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    if(dataset=="zipf") {
        xlabel <- TeX('Parameter $a$')
        param.split = c("tmp1", "sigma")
        xlims = c(1,1.75)
        xbreaks = c(1,1.25,1.5,1.75)
        x_scale = scale_x_continuous(lim=xlims, breaks=xbreaks)
        title <- "(a) Zipf"        
    } else if(dataset=="pyp") {
        xlabel <- TeX("Parameter $\\sigma$")
        param.split = c("tmp1", "tmp2", "sigma")
        xlims = c(0,1)
        xbreaks = c(0,0.5,1)
        x_scale = scale_x_continuous(lim=xlims, breaks=xbreaks)
        title <- "(b) PYP"
    } else if(dataset=="dp") {
        xlabel <- TeX("Parameter $\\theta$")
        param.split = c("tmp1", "sigma")
        xlims = c(0,1000)
        xbreaks = c(1,10,100,1000)
        x_scale = scale_x_continuous(trans="log10")
        title <- "(b) DP"
    } else {
        xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        mutate(Intervals = factor(two_sided, c(FALSE,TRUE), c("One-sided", "Two-sided"))) %>%
        filter(unique==plot_unique,
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==plot_confidence,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, `method-unique`, posterior, n_bins, n_track, key, `include-seen`, Intervals) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    pp <- df %>%
        filter(`method-unique` %in% c(1,10,100)) %>%
        filter(sigma<1.8) %>%
        filter(posterior=="mcmc", w==plot_w, n==plot_n, method %in% method.values) %>%
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(M=factor(`method-unique`)) %>%
        ggplot(aes(x=sigma, y=value, color=M, shape=M)) +
        geom_point(alpha=1, size=2) +
        geom_line(alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
#        geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
        facet_wrap(key~., scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_y_continuous(trans='log10') +
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
    suffix_seen = ifelse(plot_seen==T, "", "_unseen")
    filename <- sprintf("figures/%s_%s_bins%s%s_unique.pdf", dataset, plot_sketch, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=2, width=6, units="in")
    pp
}

if(dataset=="zipf") {
    plot_param_lb(dataset, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE)
}


# "conformal-proportional"
method.values = c("classical", "bayesian-dp-?","bootstrap","conformal-constant", "theory") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap","Conformal (fixed)", "Lower bound (asymptotic)") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#008a65", "#007756", "#006348", "#004f39", "#003b2b", "#00281d", "#00140e")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.8)
color.scale <- cbPalette[c(4,4,4)]
shape.scale <- c(9,NA,7,3)
linetype.scale <- c(1,3)

# Default arguments
model_name = dataset
plot_unique=FALSE
plot_sketch="cms-cu"
plot_n_bins=1
plot_seen=TRUE

nu <- function(a) {
    return(a^(-1/(a-1)) * (1-1/a))
}
compute_coverage_theoretical <- function(M, M1) {
    alpha = 1.0 - plot_confidence
    a = M / M1
    alpha_eff = alpha + 2 * nu(a)
    alpha_eff = pmin(alpha_eff, 1)
    return(1.0 - alpha_eff)
}

plot_param_lb_M = function(dataset, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    if(dataset=="zipf") {
        xlabel <- TeX('Size $M$\' of calibration subsets')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(trans="log10")
        title <- "(a) Zipf"        
    } else if(dataset=="pyp") {
        xlabel <- TeX("Parameter $\\sigma$")
        param.split = c("tmp1", "tmp2", "sigma")
        x_scale = scale_x_continuous(trans="log10")
        title <- "(b) PYP"
    } else if(dataset=="dp") {
        xlabel <- TeX("Parameter $\\theta$")
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(trans="log10")
        title <- "(b) DP"
    } else {
        xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        mutate(Intervals = factor(two_sided, c(FALSE,TRUE), c("One-sided", "Two-sided"))) %>%
        filter(unique==plot_unique,
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==plot_confidence,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, `method-unique`, posterior, n_bins, n_track, key, `include-seen`, Intervals) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    M.vals <- unique(na.omit(df$`method-unique`))
    df.classical <- df %>%
        filter(method %in% c("classical", "bayesian-dp-?", "bootstrap"))
    df.classical <- do.call("rbind", lapply(1:nrow(df.classical), function(r) {
        tmp = df.classical[r,] %>% slice(rep(1:n(), each = length(M.vals)))
        tmp$`method-unique` = M.vals
        return(tmp)
    })) %>%
        mutate(Method=factor(method, method.values, method.labels))
    M <- max(df$`method-unique`, na.rm=T)
    df.theory <- tibble(`method-unique`=seq(1,M), key="Coverage") %>%
        mutate(value = compute_coverage_theoretical(M, `method-unique`),
               M = `method-unique`,
               method = "theory") %>%
        mutate(Method=factor(method, method.values, method.labels))
    pp <- df %>%
        filter(method %in% c("conformal-constant")) %>%
#        rbind(df.classical) %>%
        filter(sigma %in% c(1.5)) %>%
        filter(sigma<1.8) %>%
        filter(posterior=="mcmc", w==plot_w, n==plot_n, method %in% method.values) %>%
        rbind(df.theory) %>%       
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(M=`method-unique`) %>%
        ggplot(aes(x=M, y=value, color=Method, shape=Method, linetype=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_line(data=df.theory, aes(x=M, y=value, color=Method, shape=Method, linetype=Method), alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
#        geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
        facet_wrap(key~., scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
        scale_linetype_manual(values=linetype.scale) +
#        scale_y_continuous(trans='log10', limits=c(0.5,NA)) +
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
    suffix_seen = ifelse(plot_seen==T, "", "_unseen")
    filename <- sprintf("figures/%s_%s_bins%s%s_unique_M.pdf", dataset, plot_sketch, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=2, width=6, units="in")
    pp
}

if(dataset=="zipf") {
    plot_param_lb_M(dataset, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE)
}




dataset = "covid"

idir <- sprintf("results_hpc/%s/marginal/", dataset)
ifile.list <- list.files(idir)
results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    if("gap" %in% colnames(df)) {
        df <- df %>% select(-gap)
    }
    if(startsWith(ifile, "exp3")) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols())
    } else {
        df <- tibble()
    }
}))

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


method.values = c("classical", "bayesian-dp-?","bootstrap","conformal-adaptive", "theory") #, "conformal-bayesian-dp")
method.labels = c("Classical", "Bayesian", "Bootstrap","Conformal (adaptive)", "Lower bound (asymptotic)") #, "Conformal (Bayesian)")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#008a65", "#007756", "#006348", "#004f39", "#003b2b", "#00281d", "#00140e")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.8)
color.scale <- cbPalette[c(4,4,4)]
shape.scale <- c(8,NA)
linetype.scale <- c(1,3)


plot_param_lb_M_data = function(dataset, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE) {
    if(dataset=="zipf") {
        xlabel <- TeX('Size $M$\' of calibration subsets')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(trans="log10")
        title <- "(a) Zipf"        
    } else if(dataset=="words") {
        xlabel <- TeX('Size $M$\' of calibration subsets')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(trans="log10")
        title <- "(a) Zipf"        
    } else if(dataset=="covid") {
        xlabel <- TeX('Size $M$\' of calibration subsets')
        param.split = c("tmp1", "sigma")
        x_scale = scale_x_continuous(trans="log10")
        title <- "(a) Zipf"        
    } else { 
       xlabel <- "NA"
    }
    df <- results %>%
        mutate(Coverage = coverage, Length = length) %>%
        mutate(Intervals = factor(two_sided, c(FALSE,TRUE), c("One-sided", "Two-sided"))) %>%
        filter(unique==plot_unique,
               `include-seen`==plot_seen, sketch==plot_sketch,
               data %in% data.list, confidence==plot_confidence,
               n_bins %in% c(NA,plot_n_bins)) %>%
        separate(data, param.split, sep="-") %>%
        mutate(sigma = as.numeric(sigma)) %>%
        gather(Coverage, Length, key="key", value="value") %>%
        group_by(sigma, d, w, n, method, `method-unique`, posterior, n_bins, n_track, key, `include-seen`, Intervals) %>%
        summarise(se=2*sd(value)/sqrt(n()), value=mean(value))
    M.vals <- unique(na.omit(df$`method-unique`))
    df.classical <- df %>%
        filter(method %in% c("classical", "bayesian-dp-?", "bootstrap"))
    df.classical <- do.call("rbind", lapply(1:nrow(df.classical), function(r) {
        tmp = df.classical[r,] %>% slice(rep(1:n(), each = length(M.vals)))
        tmp$`method-unique` = M.vals
        return(tmp)
    })) %>%
        mutate(Method=factor(method, method.values, method.labels))
    M <- max(df$`method-unique`, na.rm=T)
    df.theory <- tibble(`method-unique`=seq(1,M), key="Coverage") %>%
        mutate(value = compute_coverage_theoretical(M, `method-unique`),
               M = `method-unique`,
               method = "theory",
               Method = factor(method, method.values, method.labels))
    pp <- df %>%
        filter(method %in% c("conformal-adaptive")) %>%
#        rbind(df.classical) %>%       
        filter(posterior=="mcmc", w %in% plot_w, n==plot_n, method %in% method.values) %>%
        rbind(df.theory) %>%       
        mutate(Method=factor(method, method.values, method.labels)) %>%
        mutate(M=`method-unique`) %>%
        ggplot(aes(x=M, y=value, color=Method, shape=Method, linetype=Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_line(data=df.theory, aes(x=M, y=value, color=Method, linetype=Method), alpha=0.5) +
        geom_hline(data=df.dummy, aes(yintercept=value), linetype=2, alpha=0.5) +    
        geom_hline(data=df.dummy2, aes(yintercept=value), alpha=0) +    
        geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
        facet_wrap(key~., scales="free") +
        scale_color_manual(values=color.scale) +
        scale_shape_manual(values=shape.scale) +
#        scale_y_continuous(trans='log10', limits=c(0.8,NA)) +
        scale_linetype_manual(values=linetype.scale) +
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
    suffix_seen = ifelse(plot_seen==T, "", "_unseen")
    filename <- sprintf("figures/%s_%s_bins%s%s_unique_M.pdf", dataset, plot_sketch, plot_n_bins, suffix_seen)
    pp %>% ggsave(file=filename, height=2, width=6, units="in")
    pp
}

plot_param_lb_M_data(dataset, plot_unique=TRUE, plot_n_bins=1, plot_sketch="cms-cu", plot_seen=TRUE)

library(knitr)
library(kableExtra)
library(broom)
library(stringr)
library(scales)
options(knitr.table.NA = '-')

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


compute_summary_IC = function(dat){
    # Aggregate IC results and put into a summary table
    out = aggregate(dat$ARL, list("Type" = dat$um), mean)
    out = cbind(out, aggregate(dat$ARL, list("Type" = dat$um), sd)[, 2])
    out = cbind(out, aggregate(dat$ARL, list("Type" = dat$um), function(x) mean(x <= dat$Arl0[1]))[, 2])
    colnames(out) = c("Type", "AARL", "SDARL", "$\\text{Pr}(\\text{ARL}_0 \\leq \\text{a})$")

    tex_IC <- kable(out, format="latex", booktabs=TRUE, digits = 3, row.names=FALSE, escape=FALSE, align='c', linesep = "",
        caption = "Summary of the in-control performance of the control chart using the fixed-parameter, adaptive estimator, and the proposed cautious learning update rules."
    ) %>%
        kable_styling(latex_options = "hold_position")
    return(list("df" = out, "tab" = tex_IC))
}

compute_summary_OC = function(dat){
    out = aggregate(dat$ARL, list("Type" = dat$um, "delta"=dat$delta, "tau"=dat$tau), summary)
    colnames(out) = c("Type", "delta", "tau", "")
    # rearrange
    out = cbind(out[, c(3, 2, 1)], out[, 4])
    # vecBest = c()
    # uniquetypes = length(unique(out$Type))
    # for(i in 1:(NROW(out)/uniquetypes)){
    #     idx = (i-1)*uniquetypes + 1:3
    #     best = (i-1)*uniquetypes + which.min(out$AARL[idx])
    #     vecBest = c(vecBest, idx == best)
    # }
    # out$AARL = round(out$AARL, 2)
    # out$AARL = cell_spec(out$AARL, bold = vecBest, format="latex")

    tex_OC <- kable(out, format="latex", booktabs=TRUE, digits = 3, row.names=FALSE, escape=FALSE, align='c', linesep = "",
        caption = "Summary of the out-of-control performance of the control chart using the fixed-parameter, adaptive estimator, and the proposed cautious learning update rules."
    ) %>%
        kable_styling(latex_options = "hold_position")

    tex_OC = collapse_rows(tex_OC, columns = 1:2, valign = "top", latex_hline = "custom", custom_latex_hline = 2)
    return(list("df" = out, "tab" = tex_OC))
}

plotProfiles = function(dat, only_oc = TRUE){
    # Plot median CARL profiles for each control chart
	#
    # @param dat: output dataset
	#
    # @return TODO

    require(dplyr)
    require(latex2exp)
    require(ggpubr)
    require(gridExtra)
    proc = dat %>% 
        group_by(um, delta, tau) %>% 
        summarise_at(vars(ARL), list(stat = median))

    ut = unique(proc$tau)
    ut = setdiff(ut, 0)
    nt = length(ut)
    plotList = vector(mode="list", length=nt)
    for(i in 1:nt){
        tau = ut[i]
        df_plot = proc[proc$tau == tau, ]
        if(!only_oc){
            df_plot = rbind(df_plot, proc[proc$delta == 0, ])
        }
        p <- ggplot(df_plot, aes(x = delta, y = stat, colour = um, shape = um)) +
            xlab(TeX("$\\delta$")) + 
            ylab("median") + 
            geom_line() +
            geom_point() +
            guides(shape = "none") +
            theme_bw() + 
            scale_color_discrete(labels=c('AE', 'MCL', "C&M", 'FE'))+
            scale_x_continuous(breaks = seq(min(df_plot$delta), max(df_plot$delta), by = 0.25)) +
            ggtitle(TeX(paste0("$\\tau = ", as.character(tau), "$"))) +
            theme(legend.title=element_blank(), legend.direction="horizontal") 
        plotList[[i]] = p
    }
    shared_legend <- get_legend(plotList[[1]])
    plt = grid.arrange(arrangeGrob(
        plotList[[1]] + theme(legend.position="none"), 
        plotList[[2]] + theme(legend.position="none"),
        plotList[[3]] + theme(legend.position="none"),
        nrow=3),
        shared_legend,
        nrow = 3,
        heights = c(10, 1, 1))

    # plt = ggarrange(plotList[[1]], plotList[[2]], plotList[[3]], ncol = 1, nrow = length(plotList), common.legend = TRUE, legend="top")+
            # labs(color = "Dose (mg)") +

            # scale_fill_discrete(labels=c('AE', 'MCL', "C&M", 'FE'))+
    return(plt)
}

library(ggplot2)
transp = 0.3
# Load data and create the same boxplot as in Capizzi and Masarotto (2020)
setwd("/home/dede/Documents/git/SPC/CautiousLearning/")

for(sims_folder in list.files(path = "data/sims", pattern = "theta*", full.names=TRUE)){
    limitsFile = paste0(sims_folder, "/output/limits_R.csv")
    outputFile = paste0(sims_folder, "/output/output_R.csv")
    outputFolder = paste0("plots/sims/", basename(sims_folder))
    dir.create(outputFolder, recursive = FALSE, showWarnings = FALSE)
    if(file.exists(limitsFile)){
        limits = read.csv(limitsFile)
        p = ggplot(limits, aes(x = L, group = um)) +
            geom_density(aes(y = ..scaled.., fill=um), alpha=transp) + 
            scale_x_continuous(expand = expansion(mult=0.05)) +
            theme_bw() + 
            theme(legend.position="top")+
            scale_fill_discrete(labels=c('AE', 'MCL', "C&M", 'FE'))+
            guides(fill=guide_legend(title=NULL))

        ggsave(paste0(outputFolder, "/limits.png"), p)
    }

    if(file.exists(outputFile)){
        allMethods = c("AE", "MCL", "C&M", "FE")
        dat = read.csv(outputFile)
        # Boxplot ARL results
        png(paste0(outputFolder, "/IC.png"), height=1600, width = 2440, res=300)
        IC = dat[dat$tau == 0, ]
        tabIC = compute_summary_IC(IC)
        uniqueMethods = unique(IC$um)
        uniqueNames = allMethods[1:length(uniqueMethods)]

        #! Check order of xlab in boxplots
        colour = alpha(hue_pal()(length(uniqueMethods)), transp)
        boxplot(ARL ~ um, data=IC, xlab="", names=rep(uniqueNames, 1), col = colour, outline=FALSE)
        abline(h = IC$Arl0[1], lty="dashed")
        dev.off()
        write(tabIC$tab, file=paste0(outputFolder, "/IC-summary.tex"))

        OC = dat[dat$tau != 0, ]
        for(d in unique(OC$delta)){
            subdf = OC[OC$delta == d, ]
            # Get limits for the boxplot
            tmp = boxplot(ARL ~ um + tau, data=subdf, xlab="", names=rep(uniqueNames, 3), outline=FALSE, plot=FALSE)
            max_y = max(tmp$stats)
            vertical_offset = 0.10*max_y
            top_position = max_y + vertical_offset
            num_up = length(unique(subdf$um))
            num_tau = length(unique(subdf$tau))
            png(paste0(outputFolder, "/delta=", sprintf("%.2f", d), ".png"), height=1600, width = 2440, res=300)
            boxplot(ARL ~ um + tau, data=subdf, xlab="", names=rep(uniqueNames, 3), col = rep(colour, 3), outline=FALSE, ylim = c(1, top_position))
            abline(v = 0.5 + num_up*(2:num_tau-1), lty="dashed")
            text(x = 2 + num_up*(1:num_tau-1), y = max_y + 0.75*vertical_offset, parse(text = paste0("tau ==", unique(OC$tau))))
            dev.off()
        }
        
        # Plot OC median CARL profiles
        only_oc = TRUE
        plt = plotProfiles(dat, only_oc)
        ggsave(paste0(outputFolder, "/OC-profiles.png"), plt)
        tabOC = compute_summary_OC(OC)
        write(tabOC$tab, file=paste0(outputFolder, "/OC-summary.tex"))
    }
}


comparison_table = function(df1, df2, name){
    # Comparison between two control charts, depending on the update mechanisms
	#
    # @param df1: TODO
    # @param df2: TODO
    # @param name: TODO
	#
    # @return TODO

    df = data.frame()
    for(um in unique(df1$Type)){
        sub1 = df1[df1$Type == um, ]
        sub2 = df2[df2$Type == um, ]
        df = rbind(df, data.frame(rep(um, nrow(sub1)), sub1$tau, sub1$delta, sub1$Median, sub2$Median))
    }
    names(df) = c("um", "delta", "tau", name)
    return(df)
}

# table with comparison between EWMA and AEWMA depending on the update mechanism
for(sims_folder_AEWMA in list.files(path = "data/sims", pattern = "*AEWMA", full.names=TRUE)){
    sims_folder_EWMA = str_replace_all(sims_folder_AEWMA, "AEWMA", "EWMA")
    reg = regex("k = \\d.\\d+, ")
    sims_folder_EWMA = str_replace_all(sims_folder_EWMA, reg, "")
    dat_AEWMA = read.csv(paste0(sims_folder_AEWMA, "/output/output_R.csv"))
    dat_EWMA = read.csv(paste0(sims_folder_EWMA, "/output/output_R.csv"))


    oc_AEWMA = compute_summary_OC(dat_AEWMA)$df
    oc_EWMA = compute_summary_OC(dat_EWMA)$df

    tab = comparison_table(oc_EWMA, oc_AEWMA, c("EWMA", "AEWMA"))
    print(tab)
}



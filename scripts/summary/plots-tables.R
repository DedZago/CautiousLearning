library(knitr)
library(kableExtra)
library(broom)
library(stringr)
options(knitr.table.NA = '-')

compute_summary_IC = function(dat){
    # Aggregate IC results and put into a summary table
    out = aggregate(dat$ARL, list("Type" = dat$um), mean)
    out = cbind(out, aggregate(dat$ARL, list("Type" = dat$um), sd)[, 2])
    out = cbind(out, aggregate(dat$ARL, list("Type" = dat$um), function(x) mean(x <= dat$Arl0[1]))[, 2])
    colnames(out) = c("Type", "AARL", "SDARL", "$\\text{Pr}(\\text{ARL}_0 \\leq \\text{a})$")

    tex_IC <- kable(out, format="latex", booktabs=TRUE, digits = 2, row.names=FALSE, escape=FALSE, align='c', linesep = "",
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

    tex_OC <- kable(out, format="latex", booktabs=TRUE, digits = 2, row.names=FALSE, escape=FALSE, align='c', linesep = "",
        caption = "Summary of the out-of-control performance of the control chart using the fixed-parameter, adaptive estimator, and the proposed cautious learning update rules."
    ) %>%
        kable_styling(latex_options = "hold_position")

    tex_OC = collapse_rows(tex_OC, columns = 1:2, valign = "top", latex_hline = "custom", custom_latex_hline = 2)
    return(list("df" = out, "tab" = tex_OC))
}


library(ggplot2)
library(ggformula)
# Load data and create the same boxplot as in Capizzi and Masarotto (2020)
setwd("/home/dede/Documents/git/SPC/CautiousLearning/")
for(sims_folder in list.files(path = "data/sims", pattern = "theta*", full.names=TRUE)){
    outputFile = paste0(sims_folder, "/output/output_R.csv")
    outputFolder = paste0("plots/sims/", basename(sims_folder))
    dir.create(outputFolder, showWarnings = FALSE)
    if(file.exists(outputFile)){
        dat = read.csv(outputFile)
        IC = dat[dat$tau == 0, ]
        tabIC = compute_summary_IC(IC)

        p = ggplot(IC, aes(x = L, group = um)) +
            geom_density(aes(y = ..scaled.., fill=um), alpha=0.3) + 
            scale_x_continuous(expand = expansion(mult=0.05)) +
            theme_bw() + 
            theme(legend.position="top")+
            guides(fill=guide_legend(title=NULL))

        ggsave(paste0(outputFolder, "/limits.png"), p)

        # Boxplot ARL results
        png(paste0(outputFolder, "/IC.png"))
        boxplot(ARL ~ um, data=IC, xlab="", names=rep(c("AE", "CL", "FP"), 1), outline=FALSE)
        abline(h = IC$Arl0[1], lty="dashed")
        dev.off()
        write(tabIC$tab, file=paste0(outputFolder, "/IC-summary.tex"))

        OC = dat[dat$tau != 0, ]
        for(d in unique(OC$delta)){
            subdf = OC[OC$delta == d, ]
            # Get limits for the boxplot
            tmp = boxplot(ARL ~ um + tau, data=subdf, xlab="", names=rep(c("AE", "CL", "FP"), 3), outline=FALSE, plot=FALSE)
            max_y = max(tmp$stats)
            vertical_offset = 0.10*max_y
            top_position = max_y + vertical_offset
            num_up = length(unique(subdf$um))
            num_tau = length(unique(subdf$tau))
            png(paste0(outputFolder, "/delta=", sprintf("%.2f", d), ".png"))
            boxplot(ARL ~ um + tau, data=subdf, xlab="", names=rep(c("AE", "CL", "FP"), 3), outline=FALSE, ylim = c(1, top_position))
            abline(v = 0.5 + num_up*(2:num_tau-1), lty="dashed")
            text(x = 2 + num_up*(1:num_tau-1), y = max_y + 0.75*vertical_offset, parse(text = paste0("tau ==", unique(OC$tau))))
            dev.off()
        }
        tabOC = compute_summary_OC(OC)
        write(tabOC$tab, file=paste0(outputFolder, "/OC-summary.tex"))
    }
}

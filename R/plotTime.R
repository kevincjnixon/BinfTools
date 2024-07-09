.time_sum<-function (data, eb)
{
  summary_func <- function(x, col, eb) {
    sum <- NULL
    if (eb == "sd") {
      sum <- c(mean = mean(as.numeric(x[[col]]), na.rm = TRUE),
               eb = sd(as.numeric(x[[col]]), na.rm = TRUE))
    }
    if (eb == "se") {
      sum <- c(mean = mean(as.numeric(x[[col]]), na.rm = TRUE),
               eb = (sd(as.numeric(x[[col]]), na.rm = TRUE)/sqrt(length(x[[col]]))))
    }
    if (eb == 0) {
      sum <- c(mean = mean(as.numeric(x[[col]]), na.rm = TRUE),
               eb = 0)
    }
    return(sum)
  }
  if (eb == 0) {
    message("eb=0. Error bars will not be shown...")
  }
  data_sum <- plyr::ddply(data, c("group","condition","time"), .fun = summary_func,
                          "exp", eb)
  return(data_sum)
}
.testNorm<-function(sum){
  refTime<-NULL
  for(i in unique(sum$time)){
    tmp<-subset(sum, time==i)
    if(var(tmp$mean)==0 && var(tmp$eb)==0){
      refTime<-i
    }
  }
  return(refTime)
}
#' Plot time course data
#'
#' This function accepts a single gene or vector of multiple genes along with counts, conditions, and timepoints
#' to plot the difference of gene expression between conditions over time.
#'
#' @param genes character vector of at least one gene found in rownames(counts)
#' @param counts data frame of gene expression where rows are genes and columns are samples
#' @param condition character vector describing the condition of each sample in the same order of colnames(counts)
#' @param timePoints character vector describing the time points of each sample in the same order of colnames(counts)
#' @param tOrd character vector describing the order in which the timePoints should be plotted. Must be length of length(unique(timePoints))
#' @param title character describing the title of the plot
#' @param cols Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"
#' @param yax character indicating how the y-axis should be labeled - note that error bar notation will be added automatically if applicable
#' @param xax character indicating x-axis label
#' @param showStat Boolean indicating if significance should be shown on the plot.
#' @return plot showing expression changes over time
#' @export

plotTime<-function(genes, counts, condition, timePoints, tOrd, title="", cols="Dark2",
                   yax="Expression", xax="Time Point", eb="se", showStat=T){
  require(ggplot2)
  require(rstatix)
  require(ggpubr)
  require(dplyr)
  counts<-counts[which(rownames(counts) %in% genes),]
  singGene<-nrow(counts)==1
  singRep<-length(unique(paste0(condition, timePoints)))==ncol(counts)
  p<-NULL
  if(isTRUE(singGene)){
    message("Plotting a single gene...")
    if(isTRUE(singRep)){
      message("No replicates, plotting single replicate...")
      dat<-data.frame(exp=t(counts[1,]), condition=condition, time=factor(timePoints, levels=tOrd))
      colnames(dat)[1]<-"exp"
      p<-ggplot(data=dat, aes(x=time, y=exp, colour=condition, group=condition)) +
        geom_point() + geom_line() + labs(title=title, y=yax, x=xax) +
        scale_colour_manual(values=BinfTools::colPal(cols)) + theme_minimal()
    } else {
      message("Replicates, plotting mean +/-",eb)
      dat<-data.frame(exp=t(counts[1,]), condition=condition, time=factor(timePoints, levels=tOrd))
      dat$group<-paste(dat$condition, dat$time, sep="_")
      colnames(dat)[1]<-"exp"
      t<-.time_sum(dat, eb)
      p<-ggplot(data=t, aes(x=time, y=mean, colour=condition, group=condition)) +
        geom_point() + geom_line() + labs(title=title, y=paste0(yax, " (+/-",eb,")"), x=xax) +
        geom_errorbar(ggplot2::aes(ymin = mean - eb, ymax = mean + eb), width = 0.2) +
        scale_colour_manual(values=BinfTools::colPal(cols)) + theme_minimal()
      normRef<-.testNorm(t)
      if(!is.null(normRef)){
        message(normRef," time point is equal in all conditions, assuming data is normalized to this time point...")
        dat<-subset(dat, time != normRef)
      }
      pwc<- dat %>% group_by(time) %>% pairwise_wilcox_test(exp~condition) %>% adjust_pvalue(method="BH") %>%
        add_significance("p.adj")
      pwc <- pwc %>% add_xy_position(x="time")
      if(isTRUE(showStat)){
        p<-p+stat_pvalue_manual(pwc, label="p.adj.signif",tip.length = 0, inherit.aes = FALSE)
      }
    }
  } else {
    message("Plotting group of ", nrow(counts)," genes...")
    if(isTRUE(singRep)){
      message("No replicates, plotting single replicate...")
      dat<-tidyr::gather(counts, key="group", value="exp") %>% mutate(condition=as.factor(c(rep(condition, each=nrow(counts))))) %>%
        mutate(time=factor(c(rep(timePoints, each=nrow(counts))), levels=tOrd))
      t<-.time_sum(dat, eb)
      p<-ggplot(data=t, aes(x=time, y=mean, colour=condition, group=condition)) +
        geom_point() + geom_line() + labs(title=title, y=paste0(yax, " (+/-",eb,")"), x=xax) +
        geom_errorbar(ggplot2::aes(ymin = mean - eb, ymax = mean + eb), width = 0.2) +
        scale_colour_manual(values=BinfTools::colPal(cols)) + theme_minimal()
      normRef<-.testNorm(t)
      if(!is.null(normRef)){
        message(normRef," time point is equal in all conditions, assuming data is normalized to this time point...")
        dat<-subset(dat, time != normRef)
      }
      pwc<- dat %>% group_by(time) %>% pairwise_wilcox_test(exp~condition) %>% adjust_pvalue(method="BH") %>%
        add_significance("p.adj")
      pwc <- pwc %>% add_xy_position(x="time", fun=paste0("mean_",eb))
      if(isTRUE(showStat)){
        p<-p+stat_pvalue_manual(pwc, label="p.adj.signif",tip.length = 0, inherit.aes = FALSE)
      }
    } else {
      message("Replicates...")
      dat<-tidyr::gather(counts, key="group", value="exp") %>% mutate(condition=as.factor(c(rep(condition, each=nrow(counts))))) %>%
        mutate(time=factor(c(rep(timePoints, each=nrow(counts))), levels=tOrd))
      dat$group<-paste(dat$condition, dat$time, sep="_")
      t<-.time_sum(dat, eb)
      p<-ggplot(data=t, aes(x=time, y=mean, colour=condition, group=condition)) +
        geom_point() + geom_line() + labs(title=title, y=paste0(yax, " (+/-",eb,")"), x=xax) +
        geom_errorbar(ggplot2::aes(ymin = mean - eb, ymax = mean + eb), width = 0.2) +
        scale_colour_manual(values=BinfTools::colPal(cols)) + theme_minimal()
      normRef<-.testNorm(t)
      if(!is.null(normRef)){
        message(normRef," time point is equal in all conditions, assuming data is normalized to this time point...")
        dat<-subset(dat, time != normRef)
      }
      pwc<- dat %>% group_by(time) %>% pairwise_wilcox_test(exp~condition) %>% adjust_pvalue(method="BH") %>%
        add_significance("p.adj")
      pwc <- pwc %>% add_xy_position(x="time", fun=paste0("mean_",eb))
      if(isTRUE(showStat)){
        p<-p+stat_pvalue_manual(pwc, label="p.adj.signif",tip.length = 0, inherit.aes = FALSE)
      }
    }
  }
  print(p)
}

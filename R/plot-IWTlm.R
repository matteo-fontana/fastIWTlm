###################################################################################################################################################
# ggplot 2 for outputs of a IWTlm object
# 
# By: Matteo Fontana
# last modified: 30/07/2018
###################################################################################################################################################

# Inputs: 
# x:            Object of class 'IWTlm' obtained from function IWTlm
# xrange:       Range of the plotted x axis. Default is c(0,1) #should probably override this
# alpha1:       First level of significance used to select and display significant differences. 
#               Default is alpha1 = 0.05.
# alpha2:       Second level of significance used to select and display significant differences. 
#               Default is alpha2 = 0.01.
# plot.adjpval: A logical indicating wether the plots of adjusted p-values have to be done. 
#               Default is plot.adjpval = FALSE.
# col:          Vector of colors for the plot of functional data (first element), and functional coefficients (following elements).
#               Default is c(1,rainbow(dim(x$adjusted_pval_part)[1]))
# ylim:         Range of the y axis. If set to NULL, the range of the data evaluations is used.
#               Default is ylim = NULL.
# ylab:         Label of y axis of the plot of functional data. Default is "Functional Data".
# main:         An overall title for the plots (it will be pasted to "Functional Data and F-test" for the first plot 
#               and "t-test" for the other plots).
# lwd:          Line width for the plot of functional data. Default is lwd=1.
# type:         Type of lines to be used for the plot of functional data. Default is type='l'.
# ...:          Additional plotting arguments that can be used with function plot, such as graphical parameters (see par).


# Value:
# No value returned. The function produces a graphical output of the IWT results: 
# the plot of the functional data, functional regression coefficients, and IWT-adjusted p-values 
# for the F-test and t-tests. The basis components selected as significant by the tests at level 
# alpha1 and alpha2 are highlighted in the plot of the corrected p-values and in the one of 
# functional data by gray areas (light and dark gray, respectively). The plot of functional data 
# reports the gray areas corresponding to a significant F-test. The plots of functional regression 
# coefficients report the gray areas corresponding to significant t-tests for the corresponding 
# covariate.

autoplot.IWTlm <- function(object,

                      prompt=TRUE,
                       save_output=FALSE,
                       layout_vertical=FALSE,
                       xrange = c(0,1),
                       alpha1 = 0.05, 
                       alpha2 = 0.01, 
                       plot_adjpval = FALSE,
                       plot_unadjpval = FALSE,
                       ylim = NULL,
                       ylimcoeff= NULL,
                       ylab='Functional Data', 
                       xlab='abscissa_smooth',
                       main=NULL, 
                       lwd = 1, 
                       type='l',
                       cextotal=NULL,
                       ...) {
  
  if (!inherits(object, "IWTlm")) stop("object should be of the class IWTlm")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  if (plot_adjpval==FALSE & plot_unadjpval){
    warning('Plotting unadjusted p-values only is not recommended since domain selection is performed with adjusted p-values. The plots will include also adjusted p-values.')
    plot_adjpval = TRUE
  }
  p <- length(object$unadjusted_pval_F)
  J <- p
  n <- dim(object$data.eval)[1]
  minx <- xrange[1]
  maxx <- xrange[2]
  abscissa_pval <- seq(minx, maxx, len = p)
  abscissa_smooth <- seq(minx, maxx, len = J)
  
  devAskNewPage(ask = F)  
  
  main_F <- paste(main,': Functional Data and F-test')
  main_F <- sub("^ : +", "", main_F)
  
  

  
###set theme used for plotting:
  IWTtheme=theme_minimal()+
    theme(strip.background = element_blank(), strip.text.y = element_blank(), 
          panel.grid.major = element_line(color='grey80'),
          panel.grid.minor = element_line(color='grey80'),)
  

  
#####DATA PLOT#####
  
  diff_df=data.table(object$adjusted_pval_F,min_rect=abscissa_smooth)
  diff_df[,max_rect:=min_rect+(abscissa_pval[2] - abscissa_pval[1])]
  diff_df=melt(diff_df,measure.vars = 1,id.vars=c('min_rect','max_rect'))
  
  diff_df[value<alpha1,significance:=alpha1]
  diff_df[value<alpha2,significance:=alpha2]
  diff_df[min_rect>maxx,min_rect:=NA]
  diff_df[max_rect>maxx,max_rect:=NA]
  
  diff_df=na.omit(diff_df)
  
  diff_df=diff_df[,significance:=factor(significance,c(alpha2,alpha1))]
  
  
  
  data_dt=melt(data.table(t(object$data.eval), abscissa_smooth),id.vars='abscissa_smooth',measure.vars = 1:n)
  
  data_plot=ggplot(data=data_dt,aes(x=abscissa_smooth,y=value,group=variable)) +
          geom_rect(inherit.aes = F, data=diff_df, mapping=aes(xmin=min_rect,xmax=max_rect, ymin=-Inf,ymax=Inf,alpha=significance),size=0,col=NA) + 
          #geom_rect(inherit.aes = F, data=dt_rectF2, mapping=aes(xmin=min_rect,xmax=max_rect, ymin=-Inf,ymax=Inf),fill='blue',alpha=.3) + 
          scale_alpha_discrete(range=c(.6,.3),drop=FALSE) + 
          geom_line(size=.8) + 
          labs(title=main_F,y=ylab,x=xlab) +
          IWTtheme
            
  




data_plot  

#####COEFFICIENTS PLOT#####

npar=length(rownames(object$coeff.regr.eval))
coeff_name=rownames(object$coeff.regr.eval)
coeff_dt=data.table(t(object$coeff.regr.eval),abscissa_smooth)

coeff_dt=melt(coeff_dt,id.vars='abscissa_smooth',measure.vars = 1:npar)

    
main_t <- paste(main, ': Functional Coefficients & t-Tests')
main_t <- sub("^ : +", "", main_t)


diff_df=data.table(t(object$adjusted_pval_part),min_rect=abscissa_smooth)
diff_df[,max_rect:=min_rect+(abscissa_pval[2] - abscissa_pval[1])]
diff_df=melt(diff_df,measure.vars = 1:npar,id.vars=c('min_rect','max_rect'))

diff_df[value<alpha1,significance:=alpha1]
diff_df[value<alpha2,significance:=alpha2]

diff_df[min_rect>maxx,min_rect:=NA]
diff_df[max_rect>maxx,max_rect:=NA]

diff_df=na.omit(diff_df)
diff_df=diff_df[,significance:=factor(significance,c(alpha2,alpha1))]



coeff_plot=ggplot(data=coeff_dt,aes(x=abscissa_smooth,y=value,color=variable)) +
  geom_rect(inherit.aes = F, data=diff_df, mapping=aes(xmin=min_rect,xmax=max_rect, ymin=-Inf,ymax=Inf,alpha=significance,fill=variable),na.rm = T) + 
  scale_fill_discrete(drop=F)+
  scale_alpha_discrete(range=c(.6,.3),drop=FALSE) + 
  geom_line(size=1) + 
  facet_grid(variable~.) +
  labs(title=main_t,y='Functional Coefficient',x=xlab) +
  IWTtheme


coeff_plot

main_pval_t <- paste(main, ': t-Test P-Value Functions')
main_pval_t <- sub("^ : +", "", main_pval_t)



pval_dt_adj=data.table(t(object$adjusted_pval_part),abscissa_smooth,type='adjusted')
pval_dt_nonadj=data.table(t(object$unadjusted_pval_part),abscissa_smooth,type='unadjusted')
pval_dt_tot=rbind(pval_dt_adj,pval_dt_nonadj)

pval_dt_tot=melt(pval_dt_tot,measure.vars=1:npar,id.vars = c('abscissa_smooth','type'))
pval_dt_tot=pval_dt_tot[,type:=as.factor(type)]


pval_plot_t=ggplot(data=pval_dt_tot ,aes(x=abscissa_smooth,y=value,color=variable))+
            scale_linetype_manual(name='type',values=c('unadjusted'='dashed','adjusted'='solid')) +
            geom_line(aes(linetype=type),size=1) +
            facet_grid(variable~.) +
            geom_hline(yintercept = alpha1,alpha=.4) +
            geom_hline(yintercept = alpha2,alpha=.6) +
            labs(title=main_pval_t,y='P-Value Function',x=xlab) +
            ylim(0,1) +
            IWTtheme





#plot f-test pvalues

main_pval_f <- paste(main, ': F-Test P-Value Function')
main_pval_f <- sub("^ : +", "", main_pval_f)

pval_df_adj=data.table(value=object$adjusted_pval_F,abscissa_smooth,type='adjusted')
pval_df_nonadj=data.table(value=object$unadjusted_pval_F,abscissa_smooth,type='unadjusted')
pval_df_tot=rbind(pval_df_adj,pval_df_nonadj)
pval_df_tot=pval_df_tot[,type:=as.factor(type)]

pval_plot_f=ggplot(data=pval_df_tot ,aes(x=abscissa_smooth,y=value))+
  scale_linetype_manual(name='type',values=c('unadjusted'='dashed','adjusted'='solid')) +
  geom_line(aes(linetype=type),size=1) +
  geom_hline(yintercept = alpha1,alpha=.4) +
  geom_hline(yintercept = alpha2,alpha=.6) +
  labs(title=main_pval_f,y='P-Value Function',x=xlab) +
  ylim(0,1)+
  IWTtheme



plot_list=list(data_plot=data_plot,coeff_plot=coeff_plot,pval_t_plot=pval_plot_t,pval_plot_F=pval_plot_f)

if(prompt){par(ask=T)}
plot(data_plot)
plot(coeff_plot)
plot(pval_plot_f)
plot(pval_plot_t)

if(save_output) return(plot_list)
}

# S3 plot() method: convenience alias so plot(model) dispatches to autoplot.IWTlm().
plot.IWTlm <- function(x, ...) {
  autoplot.IWTlm(x, ...)
}







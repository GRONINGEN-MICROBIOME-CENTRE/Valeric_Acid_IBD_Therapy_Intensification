suppressPackageStartupMessages({library(survival);library(ggplot2);library(dplyr)})
args <- commandArgs(FALSE)
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', args, value=TRUE)), fixed=TRUE))
root <- dirname(script)
input <- Sys.getenv('VALERIC_DATA_DIR', file.path(root, 'data', 'private'))

m <- read.csv(file.path(input,'metadata.csv'),check.names=FALSE)
s <- read.csv(file.path(input,'SCFA_log2.csv'),check.names=FALSE)
d <- merge(m,s,by='ResearchID',sort=FALSE)
# Categorical diagnosis reproduces all 24 displayed P values in the final May PDF.
d$diagnosis_at_sampling <- factor(d$diagnosis_at_sampling)
cov5 <- c('ActiveDisease_at_sampling','diagnosis_at_sampling','age_at_sampling','resection_prior_sampling','clinical_BowelMovementADayDef')
outcomes <- c('therapy_intensification','medication_intensification','intensive_escalation')
times <- c('time_to_therapy_intensification','time_to_medication_intensification','time_to_intensive_escalation_m')
rows <- list()
for(i in seq_along(outcomes)) for(feat in setdiff(names(s),'ResearchID')) {
 z <- d[,c(times[i],paste0('event_',outcomes[i]),feat,cov5)]
 names(z)[1:3] <- c('time','event','x')
 z <- z[complete.cases(z),]
 fit <- coxph(reformulate(c('x',cov5),response='Surv(time,event)'),data=z,ties='efron')
 su <- summary(fit)
 rows[[length(rows)+1]] <- data.frame(outcome=outcomes[i],feature=feat,n=nrow(z),events=sum(z$event),HR=su$conf.int['x',1],lower95=su$conf.int['x',3],upper95=su$conf.int['x',4],p_value=su$coefficients['x',5])
}
r <- do.call(rbind,rows);r <- r[order(r$outcome,r$p_value),]
out <- file.path(root,'outputs','figureS2');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(r,file.path(out,'continuous_cox_five_covariates_categorical_diagnosis.csv'),row.names=FALSE)
feature_order <- c('Butyric_acid','Isovaleric_acid','Isobutyric_acid','Acetic_acid','X2.Methylbutyric_acid','Propionic_acid','Valeric_acid','Hexanoic_acid')
r$feature <- factor(r$feature,levels=rev(feature_order))
r$outcome <- factor(r$outcome,levels=outcomes)
r$sig <- ifelse(r$p_value<.05,'p < 0.05','ns')
pd <- position_dodge(width=.65)
p <- ggplot(r,aes(feature,HR,group=outcome))+geom_hline(yintercept=1,linetype='dashed',color='lightgray')+
 geom_errorbar(aes(ymin=lower95,ymax=upper95,color=sig),width=.2,linewidth=.8,position=pd)+
 geom_point(aes(shape=outcome,color=sig),size=3,position=pd)+
 geom_text(aes(y=upper95*1.03,label=sprintf('p=%.2f',p_value)),position=pd,hjust=0,size=2.65,color='gray45',show.legend=FALSE)+
 coord_flip()+scale_y_log10(limits=c(.5,3))+
 scale_x_discrete(labels=function(x) gsub('X2.Methylbutyric','2-Methylbutyric',gsub('_',' ',x),fixed=TRUE))+
 scale_shape_manual(values=c(15,17,16),labels=c('Therapy intensification','Medication intensification','Intensive escalation'))+
 scale_color_manual(values=c('ns'='gray45','p < 0.05'='darkblue'))+
 labs(x=NULL,y='Hazard Ratio (log scale)',shape='Outcome',color=NULL)+theme_classic(base_size=12)+theme(axis.text.y=element_text(size=11),legend.position='right')
ggsave(file.path(out,'FigureS2_continuous_cox_five_covariates.pdf'),p,width=10,height=7)
ggsave(file.path(out,'FigureS2_continuous_cox_five_covariates.png'),p,width=10,height=7,dpi=160)

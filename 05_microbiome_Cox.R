suppressPackageStartupMessages({library(survival);library(dplyr);library(readr)})
args <- commandArgs(FALSE)
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',args,value=TRUE)),fixed=TRUE))
root <- dirname(script)
input <- Sys.getenv('VALERIC_DATA_DIR',file.path(root,'data','private'))
m <- read_csv(file.path(input,'metadata.csv'),show_col_types=FALSE)
m$diagnosis_at_sampling <- factor(m$diagnosis_at_sampling)
covariates <- c('diagnosis_at_sampling','resection_prior_sampling','age_at_sampling','final_reads')
outcomes <- c('therapy_intensification','medication_intensification','intensive_escalation')
times <- c('time_to_therapy_intensification','time_to_medication_intensification','time_to_intensive_escalation_m')
results <- list(); warnings <- list()
for(layer in c('species','pathways')) {
 mat <- read_csv(file.path(input,paste0(layer,'_CLR.csv')),show_col_types=FALSE)
 d <- inner_join(m,mat,by='ResearchID')
 for(i in seq_along(outcomes)) for(feature in setdiff(names(mat),'ResearchID')) {
  z <- as.data.frame(d[,c(times[i],paste0('event_',outcomes[i]),feature,covariates)])
  names(z)[1:3] <- c('time','event','x');z <- z[complete.cases(z),]
  fit <- withCallingHandlers(coxph(reformulate(c('x',covariates),response='Surv(time,event)'),data=z,ties='efron'),warning=function(w){
   warnings[[length(warnings)+1]] <<- data.frame(dataset=layer,outcome=outcomes[i],feature=feature,message=conditionMessage(w))
   invokeRestart('muffleWarning')
  })
  su <- summary(fit)
  results[[length(results)+1]] <- data.frame(dataset=layer,outcome=outcomes[i],feature=feature,n=nrow(z),events=sum(z$event),HR=su$conf.int['x',1],CI_low=su$conf.int['x',3],CI_high=su$conf.int['x',4],p_value=su$coefficients['x',5])
 }
}
r <- bind_rows(results) %>% group_by(dataset,outcome) %>% mutate(p_adj=p.adjust(p_value,'BH')) %>% ungroup()
out <- file.path(root,'outputs','tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_csv(r,file.path(out,'Cox_results_metagenomic_features.csv'))
write_csv(if(length(warnings))bind_rows(warnings) else data.frame(message=character()),file.path(out,'microbial_cox_model_warnings.csv'))

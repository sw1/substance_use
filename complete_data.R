library(tidyverse)
library(Hmisc)
library(nhanesA)
library(broom)
library(rstanarm)
library(brms)
library(nnet)
library(effects)
options(mc.cores = 4)

get_labs <- function(x,y){
  labsx <- label(x)
  labsy <- label(y)
  overlap <- c(labsx[labsx != ''],labsy[labsy != ''])
  overlap <- overlap[!duplicated(overlap)]
  if (length(overlap) == 0) return(labsx) else return(overlap)
}

get_data <- function(query,type=NULL,translate=TRUE){
  out <- NULL
  suffix <- c('',paste0('_',LETTERS[2:10]))
  yrs <- seq(1999,1999+2*length(suffix),2)
  for (s in seq_along(suffix)){
    qname <- paste0(query,suffix[s])
    d  <- suppressMessages(nhanes(qname))
    if (is.null(d)) next
    
    d$year <- yrs[s]
    
    if (is.null(out)){
      out <- as_tibble(d) %>% mutate_if(is.factor,as.character.factor) %>% mutate_if(is.integer,as.character)
      labs <- get_labs(out,d)
    }else{
      labs <- get_labs(out,d)
      out <- out %>% bind_rows(d %>% mutate_if(is.factor,as.character.factor) %>% mutate_if(is.integer,as.character)) 
    }
    update <- rep('',ncol(out))
    names(update) <- colnames(out)
    update[names(labs)] <- labs
    label(out,self=FALSE) <- update
  }
  
  if (translate){
    if (!is.null(type)){
      d_vars <- NULL
      for (s in seq_along(suffix)){
        qname <- paste0(query,suffix[s])
        d_vars_tmp  <- try(nhanesTableVars(type,qname,namesonly=TRUE),silent=TRUE)
        if (class(d_vars_tmp) == 'try-error') next
        update <- setdiff(d_vars_tmp,d_vars)
        if (length(update) > 0){
          d_vars <- c(update,d_vars)
          out <- suppressWarnings(nhanesTranslate(qname,d_vars,data=out))
        }
      }
    }
  }
  
  return(out)
}





healthcond_surv <- get_data('HSQ','Q') %>% 
  select(SEQN,year,HSD010,HSQ470,HSQ480,HSQ490) 
bp_surv <- get_data('BPQ','Q') %>%
  select(SEQN,year,BPQ020,BPQ040A,BPQ080,BPQ090D) 
cv_surv <- get_data('CDQ','Q') %>%
  select(SEQN,year,CDQ001,CDQ002,CDQ003,CDQ008,CDQ010) 
dm_surv <- get_data('DIQ','Q') %>%
  select(SEQN,year,DIQ010,DIQ050,DIQ070,DIQ080,DIQ090,DIQ100,DIQ110,DIQ120,DIQ130,DIQ140,DIQ150) 
  mutate_at(vars(DIQ080:DIQ150),funs(ifelse(.=='Yes',1,0))) %>%
  mutate(periphsx=rowSums(select(.,DIQ080:DIQ150),na.rm=TRUE)) %>%
  mutate(periphsx=ifelse(periphsx > 0,1,0)) %>%
  select(SEQN:DIQ070,periphsx)
food_surv <- get_data('FSQ','Q') %>%
  select(SEQN,year,FSD032A)
vax_surv <- get_data('IMQ','Q') %>%
  select(SEQN,year,IMD010,IMQ020,IMQ030) 
conds_surv <- get_data('MCQ','Q') %>%
  select(SEQN,year,MCQ010,MCQ040,MCQ050,MCQ053,MCQ080,MCQ092,
         starts_with('MCQ120'),MCQ140,starts_with('MCQ160'),starts_with('MCQ170'),
         starts_with('MCQ180'),MCQ220,starts_with('MCQ230')) 
oral_surv <- get_data('OHQ','Q') %>%
  select(SEQN,year,OHQ030,OHQ630,OHQ680,OHQ770) 
bone_surv <- get_data('OSQ','Q') %>%
  select(SEQN,year,OSQ060,OSQ070,starts_with('OSQ010'))
  mutate_at(vars(starts_with('OSQ010')),funs(case_when(.=='1' ~ 1,
                                                       .=='2' ~ 0,
                                                       TRUE ~ NA_real_))) %>%
  mutate(brake=rowSums(select(.,starts_with('OSQ010')),na.rm=TRUE)) %>%
  mutate(brake=ifelse(brake > 0,1,0)) %>%
  select(SEQN:OSQ070,brake)
pe_surv <- get_data('PAQ','Q') %>%
  select(SEQN,year,PAD020,PAQ180) 
repro_surv <- get_data('RHQ','Q') %>%
  select(SEQN,year,RHD280)
resp_surv <- get_data('RDQ','Q') %>%
  select(SEQN,year,RDQ031,RDQ050,RDQ070,RDQ134,RDQ140) 
sex_surv <- get_data('SXQ','Q') %>%
  select(SEQN,year,SXQ260,SXQ265,SXQ270,SXQ272) 
smoke_surv <- get_data('SMQ','Q') %>%
  select(SEQN,year,
         SMQ020,SMD030,SMD055,SMD057,SMD070,SMD075,SMD650) 
support_surv <- get_data('SSQ','Q') %>%
  select(SEQN,year,SSQ011)
hear_surv <- get_data('AUQ','Q') %>%
  select(SEQN,year,AUQ130,AUQ140,AUQ150,AUQ160,AUQ190) 
hear_exam <- get_data('AUX','EXAM') %>% 
  select(SEQN,year,AUXOTSPL,AUXROTSP,AUXTMEPR,AUXTMEPL,AUAREQC,AUALEQC) 
bp_exam <- get_data('BPX','EXAM') %>%
  select(SEQN,year,BPXPLS,starts_with('BPXSY'),starts_with('BPXDI')) 
body_exam <- get_data('BMX','EXAM') %>%
  select(SEQN,year,BMXWT,BMXHT,BMXBMI)
feno_exam <- get_data('ENX','EXAM') %>%
  select(SEQN,year,ENXMEAN) 
spiro_exam <- get_data('SPX','EXAM') %>%
  select(SEQN,year,SPXNFEV1,SPXNFET,SPXNFVC,SPXBFEV1,SPXBFVC,SPXBFET) 
liver_exam <- get_data('LUX','EXAM') %>%
  select(SEQN,year,LUXCAPM,LUXCPIQR) 
vis_exam <- get_data('VIX','EXAM') %>%
  select(SEQN,year,VIDRVA,VIDLVA,VIDROVA,VIDLOVA) 
qual_surv <- get_data('HUQ','Q') %>%
  select(SEQN,year,HUQ010,HUQ020,HUQ030,HUQ040,HUQ050,HUQ060,HUQ071,HUD080,HUQ090)

dat <- hear_exam %>%
  full_join(bp_exam,by=c('SEQN','year')) %>%
  full_join(body_exam,by=c('SEQN','year')) %>%
  full_join(feno_exam,by=c('SEQN','year')) %>%
  full_join(spiro_exam,by=c('SEQN','year')) %>%
  full_join(liver_exam,by=c('SEQN','year')) %>%
  full_join(vis_exam,by=c('SEQN','year')) %>%
  full_join(healthcond_surv,by=c('SEQN','year')) %>%
  full_join(bp_surv,by=c('SEQN','year')) %>%
  full_join(cv_surv,by=c('SEQN','year')) %>%
  full_join(dm_surv,by=c('SEQN','year')) %>%
  full_join(food_surv,by=c('SEQN','year')) %>%
  full_join(vax_surv,by=c('SEQN','year')) %>%
  full_join(conds_surv,by=c('SEQN','year')) %>%
  full_join(oral_surv,by=c('SEQN','year')) %>%
  full_join(bone_surv,by=c('SEQN','year')) %>%
  full_join(pe_surv,by=c('SEQN','year')) %>%
  full_join(repro_surv,by=c('SEQN','year')) %>%
  full_join(resp_surv,by=c('SEQN','year')) %>%
  full_join(sex_surv,by=c('SEQN','year')) %>%
  full_join(smoke_surv,by=c('SEQN','year')) %>%
  full_join(support_surv,by=c('SEQN','year')) %>%
  full_join(hear_surv,by=c('SEQN','year')) %>%
  full_join(qual_surv,by=c('SEQN','year'))


saveRDS(dat,'/Research/Tiffany/dat_everything.rds')

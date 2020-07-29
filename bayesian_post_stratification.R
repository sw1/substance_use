library(tidyverse)
library(rstan)

invlogit <- function(x) 1/(1 + exp(-x))

rstan_options(auto_write=TRUE)
options(mc.cores=4)

set.seed(19581103)

library(tidyverse)
library(Hmisc)
library(nhanesA)
library(broom)
library(rstanarm)
library(brms)
library(nnet)
library(effects)
library(survey)
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


opioids <- c('CODEINE','OXYCODONE','TRAMADOL','MORPHINE',
             'HYDROMORPHONE','FENTANYL','CARFENTANIL','PROPOXYPHENE',
             'HYDROCODONE','HEROIN','PENTAZOCINE','DIHYDROCODEINE',
             'MEPERIDINE','TAPENTADOL','OPIUM') # 'BUPRENORPHINE','NALOXONE'
find_drugs <- function(rxs,lookup){
  lookup <- paste0('\\<',opioids,'\\>')
  hits <- vector(length=length(rxs))
  for (i in seq_along(rxs)){
    r <- rxs[i]
    if (is.na(r) | r=='') next
    for (l in lookup){
      if (grepl(l,r)){
        hits[i] <- TRUE
        next
      }
    }
  }
  return(as.integer(hits))
}

dat <- readRDS('D:/Downloads/tifdata_paper.rds')

subtract <- function(x,y) x-y
cancerdx <- dat$cancer %>% 
  select(SEQN,year,starts_with('MCQ240')) %>% 
  left_join(dat$demo %>% select(SEQN,year,age=RIDAGEYR),by=c('SEQN','year')) %>%
  mutate_if(is.factor,as.character.factor) %>%
  mutate_at(vars(starts_with('MCQ240')),~ifelse(.=="Don't know"|.=='99999'|.=='77777',NA,.)) %>%
  mutate_at(vars(starts_with('MCQ240')),~ifelse(.=='80 years or older','9999',.)) %>%
  mutate_at(vars(starts_with('MCQ240')),as.integer) %>%
  mutate(age=as.integer(age)) %>%
  mutate_at(vars(starts_with('MCQ240')),~subtract(age,.)) %>%
  mutate_at(vars(starts_with('MCQ240')),~ifelse(. < 1,1,0)) %>%
  mutate(cancer1yr=rowSums(select(.,starts_with('MCQ240')),na.rm=TRUE)) %>%
  filter(cancer1yr>0) %>%
  select(SEQN,year,cancer1yr) 

opioids <- dat$druginfo %>%
  filter(RXDDCI1A == '57' & RXDDCI1B == '58' & RXDDCI1C %in% c('60','191')|
           RXDDCI2A == '57' & RXDDCI2B == '58' & RXDDCI2C %in% c('60','191')|
           RXDDCI3A == '57' & RXDDCI3B == '58' & RXDDCI3C %in% c('60','191')) %>%
  select(id=RXDDRGID,drug=RXDDRUG) 
# filter(!grepl('NALOXONE|BUPRENORPHINE',drug))


drugusers <- dat$otherdruguse %>%
  select(SEQN,year,starts_with('DUQ270'),starts_with('DUQ350'),starts_with('DUQ310'),starts_with('DUQ400'),
         DUQ240,DUQ250,DUQ290,DUQ330,DUQ370) %>%
  filter(DUQ240 %in% c('Yes','No')) %>%
  # filter_at(vars(DUQ250,DUQ290,DUQ330,DUQ370),all_vars(. %in% c('Yes','No'))) %>%
  mutate_at(vars(ends_with('Q')),as.numeric) %>%
  mutate_at(vars(ends_with('Q')),~ifelse(.==7777 | .==9999,NA,.)) %>%
  # filter_at(vars(ends_with('Q')),any_vars(!is.na(.))) %>%
  mutate(days_cocaine=case_when(DUQ270U=='Days' ~ DUQ270Q,
                                DUQ270U=='Weeks' ~ DUQ270Q*7,
                                DUQ270U=='Months' ~ DUQ270Q*30,
                                DUQ270U=='Year' ~ DUQ270Q*365),
         days_heroin=case_when(DUQ310U=='Days' ~ DUQ310Q,
                               DUQ310U=='Weeks' ~ DUQ310Q*7,
                               DUQ310U=='Months' ~ DUQ310Q*30,
                               DUQ310U=='Year' ~ DUQ310Q*365),
         days_meth=case_when(DUQ350U=='Days' ~ DUQ350Q,
                             DUQ350U=='Weeks' ~ DUQ350Q*7,
                             DUQ350U=='Months' ~ DUQ350Q*30,
                             DUQ350U=='Year' ~ DUQ350Q*365),
         days_iv=case_when(DUQ400U=='Days' ~ DUQ400Q,
                           DUQ400U=='Weeks' ~ DUQ400Q*7,
                           DUQ400U=='Months' ~ DUQ400Q*30,
                           DUQ400U=='Year' ~ DUQ400Q*365)) %>%
  rowwise() %>%
  mutate(drugdur=min(days_cocaine,days_heroin,days_meth,days_iv,na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(drugdurcat=case_when(drugdur <= 365 ~ '<1y',
                              drugdur > 365 & drugdur < Inf ~ '>1y',
                              DUQ240 == 'No' | (DUQ250 == 'No' & DUQ290 == 'No' & DUQ330 == 'No' & DUQ370 == 'No') ~ 'none')) %>%
  # filter(!is.na(drugdurcat)) %>%
  mutate(druguse=ifelse(drugdurcat == 'none','No','Yes')) %>%
  select(SEQN,year,druguse,drugdur,drugdurcat,starts_with('days'))


opioidusers <- dat$druguse %>% 
  select(SEQN,year,rxcount=RXD295,rxuse=RXDUSE,rx=RXDDRUG,rxid=RXDDRGID,rxdur=RXDDAYS) %>% 
  # filter(opioiduse %in% c('Yes','No')) %>%
  mutate(opioiduse=case_when((rxid %in% opioids$id) ~ 'Yes',
                             !(rxid %in% opioids$id) | rxuse == 'No' ~ 'No'),
         opioiddur=ifelse(opioiduse=='Yes',as.integer(rxdur),NA),
         opioiddurcat=case_when(opioiduse == 'No' ~ 'none',
                                opioiddur < 90 ~ 'short',
                                opioiddur >= 90 & opioiddur <= 25550 ~ 'long')) %>%
  mutate(rxcount=as.integer(rxcount),
         rxcount_listed=str_count(as.character.factor(rx),';')) 


opioidusers <- dat$druguse %>% 
  select(SEQN,year,rxcount=RXD295,rxuse=RXDUSE,rx=RXDDRUG,rxid=RXDDRGID,rxdur=RXDDAYS) %>% 
  # filter(opioiduse %in% c('Yes','No')) %>%
  mutate(opioiduse=case_when((rxid %in% opioids$id) ~ 'Yes',
                             !(rxid %in% opioids$id) | rxuse == 'No' ~ 'No'),
         opioiddur=ifelse(opioiduse=='Yes',as.integer(rxdur),NA),
         rxcount=as.integer(rxcount),
         rxcount_listed=str_count(as.character.factor(rx),';'),
         rxcount_listed=ifelse(nchar(as.character.factor(rx)) > 0,rxcount_listed + 1,rxcount_listed)) %>%
  group_by(SEQN,year) %>%
  mutate(rxcount_listed=sum(rxcount_listed),
         opioiduse=ifelse(any(opioiduse == 'Yes'),'Yes',opioiduse),
         opioiddur=ifelse(opioiddur == 99999 | opioiddur == 77777,NA,opioiddur),
         opioiddur=max(opioiddur,na.rm=TRUE),
         opioiddur=ifelse(is.infinite(opioiddur),NA,opioiddur)) %>%
  ungroup() %>%
  nest(rx=c(rx,rxdur,rxid)) %>%
  mutate(opioiddurcat=case_when(opioiduse == 'No' ~ 'none',
                                opioiddur < 90 ~ 'short',
                                opioiddur >= 90 & opioiddur <= 25550 ~ 'long'))


users <- full_join(drugusers,opioidusers,by=c('SEQN','year')) %>%
  filter((druguse %in% c('Yes','No')) | (opioiduse %in% c('Yes','No')),
         !is.na(druguse) | !is.na(opioiduse)) %>%
  left_join(dat$smoke %>% # might be wrong
              select(SEQN,year,smokecurrent=SMQ040,smokelife=SMQ020) %>%
              mutate(smoke=case_when(smokecurrent %in% c('Every day,','Some days, or') ~ 'Current',
                                     smokelife == 'Yes' ~ 'Former',
                                     smokecurrent == 'Not at all?' | smokelife == 'No' ~ 'Never')) %>%
              select(SEQN,year,smoke),
            by=c('SEQN','year')) %>%
  left_join(dat$ins %>% 
              select(SEQN,year,ins1=HID010,ins2=HIQ011) %>%
              mutate(ins=case_when(ins1 == 'Yes' | ins2 == 'Yes' ~ 'Yes',
                                   ins1 == 'No' | ins2 == 'No' ~ 'No')) %>%
              select(SEQN,year,ins),
            by=c('SEQN','year')) %>%
  left_join(dat$demo %>% select(SEQN,year,age=RIDAGEYR,sex=RIAGENDR,race=RIDRETH1,edu=DMDEDUC2) %>% 
              mutate(age=as.integer(age),
                     edu=case_when(edu %in% c('Less Than 9th Grade',
                                              '9-11th Grade (Includes 12th grad') ~ 'Less than high school' ,
                                   edu == 'High School Grad/GED or Equivale' ~ 'High school or equivalent',
                                   edu == 'Some College or AA degree' ~ 'Some college',
                                   edu == 'College Graduate or above' ~ 'College or higher'),
                     race=case_when(race %in% c('Mexican American','Other Hispanic') ~ 'Hispanic',
                                    race == 'Other Race - Including Multi-Rac' ~ 'Non-Hispanic Other',
                                    race == 'Non-Hispanic White' ~ 'Non-Hispanic White',
                                    race == 'Non-Hispanic Black' ~ 'Non-Hispanic Black')),
            by=c('SEQN','year')) %>%
  left_join(dat$demo_raw %>% select(SEQN,cycle=SDDSRVYR,year,
                                    wtint2yr=WTINT2YR,wtmec2yr=WTMEC2YR,
                                    cluster=SDMVPSU,strata=SDMVSTRA),by=c('year','SEQN')) %>%
  filter(age >= 18 & age <= 100) %>%
  left_join(dat$body %>% select(SEQN,year,bmi=BMXBMI),by=c('SEQN','year')) %>%
  mutate(bmicat=case_when(bmi < 20 ~ 'underweight',
                          bmi >= 20 & bmi < 25 ~ 'normal',
                          bmi >= 25 & bmi < 30 ~ 'overweight',
                          bmi >= 30 & bmi < 35 ~ 'obese I',
                          bmi >= 35 & bmi < 40 ~ 'obese II',
                          bmi >= 40 & bmi <= 80~ 'obese III'),
         agecat=case_when(age >= 18 & age < 25 ~ '18-24',
                          age >= 25 & age < 35 ~ '25-34',
                          age >= 35 & age < 45 ~ '35-44',
                          age >= 45 & age < 55 ~ '45-54',
                          age >= 55 & age < 65 ~ '55-64',
                          age >= 65 & age < 75 ~ '65-74',
                          age >= 75 & age < 100 ~ '75-100')) %>%
  left_join(dat$preg %>% select(SEQN,year,preg=RHD143),by=c('SEQN','year')) %>%
  mutate(preg=case_when(preg == 'Yes' ~ 'Yes',TRUE ~ 'No')) %>%
  left_join(cancerdx,by=c('SEQN','year')) %>%
  filter(preg != 'Yes'|is.na(preg)) %>%
  filter(is.na(cancer1yr)) %>%
  filter(!is.na(smoke)) %>%
  filter(!is.na(edu)) %>%
  filter(ins %in% c('Yes','No')) %>%
  mutate(sex=relevel(as.factor(sex),ref='Male'),
         opioiddurcat=relevel(as.factor(opioiddurcat),ref='none'),
         drugdurcat=relevel(as.factor(drugdurcat),ref='none'),
         bmicat=relevel(as.factor(bmicat),ref='normal'),
         agecat=relevel(as.factor(agecat),ref='35-44'),
         race=relevel(as.factor(race),ref='Non-Hispanic White'),
         edu=relevel(as.factor(edu),ref='High school or equivalent'),
         ins=relevel(as.factor(ins),ref='No'),
         smoke=relevel(as.factor(smoke),ref='Never'),
         yearcat=relevel(as.factor(year),ref='2003')) %>%
  mutate(strata=dense_rank(strata),cluster=dense_rank(cluster),
         opioiddurrank=case_when(opioiddurcat == 'none' ~ 1,
                                 opioiddurcat == 'short' ~ 2,
                                 opioiddurcat == 'long' ~ 3)) %>% 
  select(SEQN,year,yearcat,bmi,bmicat,age,agecat,sex,race,edu,ins,smoke,
         strata,cluster,
         rx,rxcount,rxcount_listed,
         drugdur,drugdurcat,druguse,starts_with('days'),
         opioiddur,opioiddurcat,opioiddurrank,opioiduse,wtint2yr,wtmec2yr) %>% 
  mutate(wtintadj=wtint2yr/length(unique(year)),wtmecadj=wtmec2yr/length(unique(year))) 


dat_everything <- readRDS('/Research/Tiffany/dat_everything.rds')
dat_everything <- dat_everything %>%
  select(union(
    dat_everything %>%
      filter(SEQN %in% (users %>% filter(opioiduse == 'Yes') %>% select(SEQN) %>% distinct() %>% unlist())) %>%
      select_if(~sum(!is.na(.)) > 100) %>% # consider changing
      colnames(),
    dat_everything %>%
      filter(SEQN %in% (users %>% filter(druguse == 'Yes') %>% select(SEQN) %>% distinct() %>% unlist())) %>%
      select_if(~sum(!is.na(.)) > 100) %>% # consider changing
      colnames()
  )) %>%
  select(-BMXBMI) %>%
  mutate_if(is.factor,as.character.factor) %>%
  mutate_all(~ifelse(. %in% c('Refused',"Don't know"),NA,.))

users_everything <- users %>%
  left_join(dat_everything,by=c('SEQN','year')) %>%
  rename(normal_otoscopy_L=AUXOTSPL,normal_otoscopy_R=AUXROTSP,
         middle_ear_pres_R=AUXTMEPR,middle_ear_pres_L=AUXTMEPL,
         tympanogram_qual_R=AUAREQC,tympanogram_qual_L=AUALEQC,
         hear_gen=AUQ130,hear_last_test=AUQ140,hear_aid=AUQ150,hear_recent_sx=AUQ190,
         pulse=BPXPLS,
         bp_sys_1=BPXSY1,bp_sys_2=BPXSY2,bp_sys_3=BPXSY3,bp_sys_4=BPXSY4,
         bp_di_1=BPXDI1,bp_di_2=BPXDI2,bp_di_3=BPXDI3,bp_di_4=BPXDI4,
         bp_high=BPQ020,bp_med=BPQ040A,
         chol_high=BPQ080,chol_med=BPQ090D,
         wt=BMXWT,ht=BMXHT,overweight1=MCQ080,overweight2=MCQ160J,
         feno=ENXMEAN,
         fev1_bl=SPXNFEV1,fet_bl=SPXNFET,fvc_bl=SPXNFVC,
         vis_R=VIDRVA,vis_L=VIDLVA,trouble_seeing=MCQ140,
         vis_acuity_R=VIDROVA,vis_acuity_L=VIDLOVA,
         gen_health1=HSD010,days_in_month_poor_health=HSQ470,days_in_month_poor_mental=HSQ480,days_in_month_forced_inactive=HSQ490,
         cp=CDQ001,cp_exert=CDQ002,cp_noexert=CDQ003,cp_severe30=CDQ008,sob_stairs=CDQ010,
         dm=DIQ010,dm_insulin=DIQ050,dm_med=DIQ070,dm_eye=DIQ080,dm_ulcer=DIQ090,
         dm_numb_recent_handfeet=DIQ100,dm_numb_loc=DIQ110,dm_paintingle_recent_handfeet=DIQ120,dm_paintingle_loc=DIQ130,
         dm_pain_leg=DIQ140,dm_pain_calf=DIQ150,
         run_out_food=FSD032A,emo_support=SSQ011,
         hepa_vax=IMD010,hepb_vax=IMQ020,liverdz1=MCQ160L,liverdz2=MCQ170L,liverdz_age=MCQ180L,liver_cap_med=LUXCAPM,liver_cap_iqr=LUXCPIQR,
         pneumo_vax=IMQ030,
         asthma=MCQ010,asthma_recent_treatment=MCQ053,asthma_recent_attack=MCQ040,asthma_recent_er=MCQ050,
         emphysema=MCQ160G,emphysema_age=MCQ180G,chr_bronch1=MCQ160K,chr_bronch2=MCQ170K,chr_bronch_age=MCQ180K,copd=MCQ160O,
         cough_chronic=RDQ031,cough_prod=RDQ050,wheeze=RDQ070,wheeze_med=RDQ134,cough_nocturnal_dry=RDQ140,
         blood_trans=MCQ092,
         arthritis=MCQ160A,arthritis_age=MCQ180A,gout=MCQ160N,gout_age=MCQ180N,
         fx_hip=OSQ010A,fx_wrist=OSQ010B,fx_spine=OSQ010C,op=OSQ060,op_tx=OSQ070,
         chf=MCQ160B,chf_age=MCQ180B,chd=MCQ160C,chd_age=MCQ180C,angina=MCQ160D,angina_age=MCQ180D,mi=MCQ160E,mi_age=MCQ180E,
         stroke=MCQ160F,stroke_age=MCQ180F,
         thyroiddz1=MCQ160M,thyroiddz2=MCQ160I,thyroiddz_current=MCQ170M,thyroiddz_age=MCQ180M,goiter=MCQ160H,
         cancer=MCQ220,cancer_type=MCQ230A,
         dentist_last_visit=OHQ030,insecure_mouth_freq=OHQ630,insecure_mouth_freq_recent=OHQ680,dentist_could_not_get_recent=OHQ770,
         herpes=SXQ260,warts=SXQ265,gono=SXQ270,chlam=SXQ272,hysterectomy=RHD280,
         gen_health2=HUQ010,gen_health_comp=HUQ020,hc_loc=HUQ030,hc_type=HUQ040,hc_count=HUQ050,hc_last=HUQ060,
         hosp_nights=HUQ071,hosp_nights_count=HUD080,mental_visit_count=HUQ090,
         activity_walkbike=PAD020,activity_avg_amount=PAQ180,
         smoke_age_reg=SMD055,smoke_age=SMD030,smoke_num_cigs_per_day_now=SMD070,smoke_num_cigs_per_day_now_length=SMD075,
         smoke_num_cigs_per_day=SMD057,smoke_avg_cigs_per_day_recent=SMD650,smoke_100cigs=SMQ020) %>%
  mutate(gen_health1=gsub(',.*|\\?.*','',gen_health1,ignore.case=TRUE),
         gen_health2=gsub(',.*|\\?.*','',gen_health2,ignore.case=TRUE),
         dm=relevel(factor(dm,ordered=FALSE),ref='No'),
         gen_health1=factor(gen_health1,ordered=TRUE,levels=rev(c('Excellent','Very good','Good','Fair','Poor'))),
         gen_health2=factor(gen_health2,ordered=TRUE,levels=rev(c('Excellent','Very good','Good','Fair','Poor'))),
         bmicat=factor(bmicat,ordered=TRUE,levels=c('underweight','normal','overweight','obese I','obese II','obese III')),
         agecat=factor(agecat,ordered=TRUE,levels=c('18-24','25-34','35-44','45-54','55-64','65-74','75-100')),
         edu=factor(edu,ordered=TRUE,levels=c('Less than high school','High school or equivalent','Some college','College or higher')),
         yearrank=dense_rank(year),
         bmirank=dense_rank(bmicat),
         edurank=dense_rank(edu),
         agerank=dense_rank(agecat)) %>%
  mutate_at(vars(starts_with('bp_sys')),as.integer) %>%
  mutate_at(vars(starts_with('bp_di')),as.integer) %>%
  rowwise() %>%
  mutate(map_mean=mean(c((bp_sys_1 + 2*bp_di_1)/3,(bp_sys_2 + 2*bp_di_2)/3,(bp_sys_3 + 2*bp_di_3)/3),na.rm=TRUE),
         map_min=min(c((bp_sys_1 + 2*bp_di_1)/3,(bp_sys_2 + 2*bp_di_2)/3,(bp_sys_3 + 2*bp_di_3)/3),na.rm=TRUE)) %>%
  ungroup() %>%
  mutate_at(vars(starts_with('map')),~ifelse(is.infinite(.)|is.nan(.),NA,.))

sample_dat <- users_everything %>%
  select(seqn=SEQN,cp,sex,race,age=agerank,edu=edurank,ins,smoke,year=yearrank,bmi=bmirank,opioid=opioiddurcat,n=wtmecadj) %>%
  filter(!is.na(cp),!is.na(bmi),!is.na(opioid)) %>%
  mutate(cp=as.integer(ifelse(cp=='Yes',1,0)),
         ins=ifelse(ins=='Yes',1,0),
         sex=ifelse(sex=='Female',0,1),
         race=case_when(race=='Hispanic' ~ 1,race=='Non-Hispanic Black' ~ 2,race=='Non-Hispanic Other' ~ 3,race=='Non-Hispanic White' ~ 4),
         smoke=dense_rank(factor(smoke,ordered=TRUE,levels=c('Never','Former','Current'))),
         opioid=dense_rank(factor(opioid,ordered=TRUE,levels=c('none','short','long'))),
         age=dense_rank(age),
         edu=dense_rank(edu),
         year=dense_rank(year),
         bmi=dense_rank(bmi)) %>%
  distinct()

pop_dat <- sample_dat %>% 
  group_by(cp,sex,race,age,edu,ins,smoke,year,bmi,opioid) %>%
  dplyr::summarize(n=sum(n)) %>%
  ungroup() %>%
  complete(cp,sex,race,age,edu,ins,smoke,year,bmi,opioid,fill=list(n=0))

stan_data <- list(
  N = nrow(sample_dat),
  N_race = max(sample_dat$race),
  N_age = max(sample_dat$age),
  N_edu = max(sample_dat$edu),
  N_smoke = max(sample_dat$smoke),
  N_year = max(sample_dat$year),
  N_bmi = max(sample_dat$bmi),
  N_opioid = max(sample_dat$opioid),
  
  y = sample_dat$cp,
  
  sex = sample_dat$sex,
  ins = sample_dat$ins,
  
  race = sample_dat$race,
  age = sample_dat$age,
  edu = sample_dat$edu,
  smoke = sample_dat$smoke,
  year = sample_dat$year,
  bmi = sample_dat$bmi,
  opioid = sample_dat$opioid
)

stan_model <- '
  data {
    int<lower=0> N;
    int<lower=0> N_race;
    int<lower=0> N_age;
    int<lower=0> N_edu;
    int<lower=0> N_smoke;
    int<lower=0> N_year;
    int<lower=0> N_bmi;
    int<lower=0> N_opioid;
    
    int<lower=0,upper=1> y[N];
    int<lower=0,upper=1> sex[N];
    int<lower=0,upper=1> ins[N];
    int<lower=0,upper=N_race> race[N];
    int<lower=0,upper=N_age> age[N];
    int<lower=0,upper=N_edu> edu[N];
    int<lower=0,upper=N_smoke> smoke[N];
    int<lower=0,upper=N_year> year[N];
    int<lower=0,upper=N_bmi> bmi[N];
    int<lower=0,upper=N_opioid> opioid[N];
  }
  parameters {                                     
    real b0;                                          
    real b_sex;
    real b_ins;
    vector[N_race] a_race;
    vector[N_age] a_age;
    vector[N_edu] a_edu;
    vector[N_smoke] a_smoke;
    vector[N_year] a_year;
    vector[N_bmi] a_bmi;
    vector[N_opioid] a_opioid;
    real<lower=0> sigma_race;
    real<lower=0> sigma_age;
    real<lower=0> sigma_edu;
    real<lower=0> sigma_smoke;
    real<lower=0> sigma_year;
    real<lower=0> sigma_bmi;
    real<lower=0> sigma_opioid;
  }
  transformed parameters {                               
    vector[N] theta;
    
    for (n in 1:N)
      theta[n] = b0 + 
        b_sex*sex[n] + b_ins*ins[n] +
        a_race[race[n]] + a_age[age[n]] + a_edu[edu[n]] + a_smoke[smoke[n]] + a_year[year[n]] + a_bmi[bmi[n]] + a_opioid[opioid[n]];
  }
  model {                                                
    b0 ~ student_t(3,0,2.5);
    b_sex ~ student_t(3,0,2.5);
    b_ins ~ student_t(3,0,2.5);
    
    a_race ~ student_t(7,0,sigma_race);
    a_age ~ student_t(7,0,sigma_age);
    a_edu ~ student_t(7,0,sigma_edu);
    a_smoke ~ student_t(7,0,sigma_smoke);
    a_year ~ student_t(7,0,sigma_year);
    a_bmi ~ student_t(7,0,sigma_bmi);
    a_opioid ~ student_t(7,0,sigma_opioid);

    sigma_race ~ student_t(3,0,2.5);
    sigma_age ~ student_t(3,0,2.5);
    sigma_edu ~ student_t(3,0,2.5);
    sigma_smoke ~ student_t(3,0,2.5);
    sigma_year ~ student_t(3,0,2.5);
    sigma_bmi ~ student_t(3,0,2.5);
    sigma_opioid ~ student_t(3,0,2.5);

    y ~ bernoulli_logit(theta);
  }
  generated quantities{
    int<lower=0,upper=1> y_hat[N];
    vector[N] log_lik;
  
    y_hat = bernoulli_logit_rng(theta);
    
    for (n in 1:N)
      log_lik[n] = bernoulli_logit_lpmf(y[n]|theta[n]);
  }
'

chains <- 4
iters <- 750
n_samp <- chains*iters/2
system.time(fit <- stan(model_code=stan_model,data=stan_data,iter=iters,chains=chains,
                        # pars=c("theta","p"),include=FALSE,
                        control = list(max_treedepth=12, adapt_delta=0.95),
                        verbose=TRUE
))


s <- summary(fit)
pars <- extract(fit)

n_cell <- nrow(pop_dat)
y_hat <- matrix(NA,n_samp,cells)
for (n in seq_len(n_cell)){
  theta <- pars$b0 + 
    pars$b_sex*pop_dat$sex[n] + pars$b_ins*pop_dat$ins[n] + 
    pars$a_race[,pop_dat$race[n]] + pars$a_age[,pop_dat$age[n]] + pars$a_edu[,pop_dat$edu[n]] + 
    pars$a_smoke[,pop_dat$smoke[n]] + pars$a_year[,pop_dat$year[n]] + pars$a_bmi[,pop_dat$bmi[n]] +
    pars$a_opioid[,pop_dat$opioid[n]]
  y_hat[,n] <- invlogit(theta)
}

opioid_hat <- matrix(NA,n_samp,stan_data$N_opioid)
for(s in seq_len(n_samp)){
  for(o in seq_len(stan_data$N_opioid)){  
    opioid_hat[s,o] <- sum(pop_dat$n[pop_dat$opioid == o]*y_hat[s,pop_dat$opioid == o])/sum(pop_dat$n[pop_dat$opioid == o])
  }
}

as_tibble(opioid_hat,.name_repair = 'unique') %>% 
  dplyr::rename('none'=1,'short'=2,'long'=3) %>%
  gather(use,p) %>%
  mutate(use=factor(use,ordered=TRUE,levels=c('none','short','long'))) %>%
  ggplot(aes(x=use,y=p)) +
  geom_violin() +
  ylim(0,1)


svy_dat <- users_everything %>%
  select(seqn=SEQN,cp,sex,race,age=agerank,edu=edurank,ins,smoke,year=yearrank,bmi=bmirank,opioid=opioiddurcat,n=wtmecadj,strata,cluster) %>%
  filter(!is.na(cp),!is.na(bmi),!is.na(opioid)) %>%
  mutate(cp=as.integer(ifelse(cp=='Yes',1,0))) %>%
  distinct()

svy <- svydesign(id=~cluster,strata=~strata,weights=~n,data=svy_dat,nest=TRUE)
m <- svyglm(cp ~ sex + race + age + edu + ins + smoke + year + bmi + opioid,
            design=svy,
            family=binomial(link='logit'),
            rescale=TRUE)
summary(m)
tibble(predict(m,type='response'),use=svy$variables$opioid) %>%
  rename('p'=1) %>%
  mutate(use=factor(use,ordered=TRUE,levels=c('none','short','long'))) %>%
  ggplot(aes(x=use,y=p)) +
  geom_violin() +
  ylim(0,1)


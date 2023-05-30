#######################################################
############  ANALYSIS OF CLint DISTRIBUTION ##########
#######################################################

# Hepatocytes and enterocytes

# Libraries 

# Libraries 
library(reshape2)
library(gridExtra)
library(ggplot2)
library(stats)
library(dplyr)
library(forecast)
library(grid)
library(data.table)   # raname columns from input file
library(AID)
library(onewaytests)
library(diptest)

# Get working directory 
folder = getwd()

## Preliminary data for the model, here the file with all activity data is read 
df <- data.frame(read.csv(paste0(folder,"/Human_Hepatocytes_data.csv"), na.strings = c("NT", "BQL", "TBD")))


# Hepatic activity distribution. In this section the activity of several enzymes/isoforms are considered and statistically analysed 

stat_analysis = function(df){
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  
  setnames (df, skip_absent=TRUE, old= c("MTT","ECOD","X7.HCG","X7.HCS","PHEN", "COUM", "BUP",    "AMO",  "TOLB",   "MEPH",   "DEX",   "CZX",    "TEST",    "MID",  "Description"), 
            new= c( "MTT","ECOD","UGTs",   "ST","PHEN","CYP2A6","CYP2B6","CYP2C8","CYP2C9","CYP2C19","CYP2D6", "CYP2E1", "CYP3A4", "CYP3A4.5", "Sex"))
  enzyme <- c("MTT","ECOD","UGTs","ST","PHEN","CYP2A6","CYP2B6","CYP2C8","CYP2C9","CYP2C19","CYP2D6", "CYP2E1", "CYP3A4", "CYP3A4.5","AO","UGT1A1")
  
  # Transformation in log domain   
  df.log <- cbind(df[,c(1:6)], log(df[,c(7:22)]),df[,c(23:38)])
  
  
  # Preliminary analysis. How many values are not available? The values are printed 
  apply(df, 2, function(x) length(which(!is.na(x))))
  
  # Generation of a dataframe including as describtion the gender: male, female, and all
  res <- setNames(data.frame(matrix(ncol =  length(enzyme)+1, nrow = 3)), c("Description",enzyme))
  res$Description <- c("ALL", "FEMALE","MALE")
  
  
  
  enz.act.data <- data.frame(matrix(ncol=20,nrow=0, dimnames=list(NULL, c("ISOFORM", "MEDIAN", "MEAN","SD","p.NORM" ,"BARTLETT","p_UNIMODAL",
                                                                          "MEAN.LOG", "SD.LOG", 
                                                                          "MEDIAN.BC","MEAN.BC", "SD.BC", "LAMBDA","p.NORM.BC","NORM.BC.male", "NORM.BC.female" ,"p.ANOVA.BC", "BARTLETT.BC",
                                                                          "N.MALE", "N.FEMALE"))))
  
  
  for (k in 1 : length(enzyme)){
    # kth is refers to isoform
    isoform <- enzyme[k]
    
    # linear distribution
    all.lin <- df[!df[,isoform] == 0,]
    all.lin <- all.lin  %>% filter(all.lin[,isoform]!="NA") 
    
    # Get the mode
    mode_lin = getmode(all.lin[,isoform])
    
    # perform normality assesment
    shapiro.lin <- shapiro.test(all.lin[,isoform])
    shapiro.lin$p.value
    
    # Log distibution
    all.log <- df.log[!df.log[,isoform] == -Inf,]
    all.log <- all.log  %>% filter(all.log[,isoform]!="NA") 
    
    # Count number of female and male with available data
    female.log <-  all.log %>% filter(Sex=="Female")
    male.log <- all.log %>% filter(Sex=="Male")
    female.lin <- all.lin %>% filter(Sex=="Female")
    male.lin <- all.lin%>% filter(Sex=="Male")
    n.female <- nrow(female.lin)
    n.male <- nrow(male.lin) 
    n.all <- n.female + n.male
    
    
    # unimodal evaluation in linear domain
    test_ks=dip.test(all.lin[,isoform])
    p_unimod=test_ks$p.value
    
    ### BOX-COX tranformation ###
    
    # Box-Cox transformation and lambda determination lin = linear distr. and log = log distr.from forcast (approx. lambda +/- 0.05)
    lambda.lin <- BoxCox.lambda(all.lin[,isoform], method = c("loglik") )
    
    
    # transformation with lambda 
    bc.transf.lin <- BoxCox( all.lin[,isoform], lambda.lin)
    
    # Mean and SD calculation after Box-Cox transformation
    box.cox.transf.mean <- mean(bc.transf.lin)
    box.cox.transf.sd <- sqrt(var(bc.transf.lin))
    
    #Get mode of BC trans. data
    mode_BC = getmode(bc.transf.lin)
    anova_df= all.lin
    anova_df$tra <- bc.transf.lin
    
    # Shapiro and Bartlett's test and  both sexes to asses normality and homogeneity test, respectively  with the lambda lambda.lin
    # test performed with AID library
    out_test <- boxcoxfr(all.lin[,isoform], all.lin$Sex, lambda = lambda.lin, tau=0)
    out_test_shapiro_male = out_test$shapiro$Normality[2]
    out_test_shapiro_female = out_test$shapiro$Normality[1]
    out_test_Bartlett_all_tra = out_test$bartlett$Homogeneity
    out_test_Bartlett_all_lin = homog.test(get(isoform) ~ Sex, data = all.lin, method = "Bartlett")
    p_lin_Bartlett = out_test_Bartlett_all_lin$p.value
    
    # Shapiro test of the Box-Cox trans. data
    shapiro.tra <- shapiro.test(out_test$tf.data)
    
    if (p_lin_Bartlett <0.05){
      p_lin_Bartlett_out = "NO HOMOGENEOUS"
      
    } else{
      p_lin_Bartlett_out = "YES HOMOGENEOUS"
    }
    
    
    ### Anova test ###
    
    # Box-Cox transformed data
    res.anova_tra <- aov(tra ~ Sex, data = anova_df)
    
    # Box-Cox linear data
    res.anova_lin <- aov(get(isoform) ~ Sex, data = anova_df)
    
    
    # p value transformed  
    p_tra <- summary(res.anova_tra)
    p_tra <- p_tra[[1]]$`Pr(>F)`[1]
    
    # p value linear  
    p_lin <- summary(res.anova_lin)
    p_lin <- p_lin[[1]]$`Pr(>F)`[1]
    
    
    
    
    
    # Mean and SD for LOG dist
    all.stat <- c( isoform, median(all.lin[,isoform]),  mean(all.lin[,isoform]), sqrt(var(all.lin[,isoform])) , shapiro.lin$p.value, p_lin_Bartlett_out,p_unimod,
                   mean(all.log[,isoform]), sqrt(var(all.log[,isoform])),
                   median(bc.transf.lin), mean(bc.transf.lin), sqrt(var(bc.transf.lin)),lambda.lin, shapiro.tra$p.value, out_test_shapiro_male,out_test_shapiro_female,p_tra, out_test_Bartlett_all_tra,
                   n.male, n.female)
    
    for (i in 1 : length(all.stat)){
      enz.act.data[k,i] <- all.stat[i]
    } 
    all.stat <- as.numeric(all.stat[-1])
    
    
    
  }
  
  # hep.enz.act.<- reshape2::melt(enz.act.data,  id.vars = c( "ISOFORM"), variable.name = c("GENDER"), value.name = "VALUE")
  write.csv(file = paste0(folder,"/out_analysis.csv"), enz.act.data)
  return(enz.act.data)
}


data_set=stat_analysis(df)



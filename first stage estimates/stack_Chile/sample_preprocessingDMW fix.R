rm(list = ls(all = TRUE))

#Install and load necessary packages
#install.packages(c('openxlsx', 'rio', 'date', 'zoo'), quiet = TRUE)
#library('openxlsx', quietly = TRUE)
#library('rio', quietly = TRUE)
library('date', quietly = TRUE)
library('zoo', quietly = TRUE)
library(lubridate)

#File information
country <- 'Chile_state'
#for(country in c('Brazil','Ecuador', 'Chile', 'Mexico', 'US' )){

if (country == 'US') {
  ICD9 = TRUE
} else {
  ICD9 = FALSE
}
file_directory <- paste('C:/Users/dmw63/Dropbox (Personal)/Meta analysis datasets/', country, '/', sep = '')
  file_input_name <- 'Chile region-level series master 16 pah sch sy.csv'
  file_output_name <- paste(country, '_processed_data.csv', sep = '')

#Convert file types
#rio::convert(paste(file_directory, file_input_name, sep = ''), paste(file_directory, file_output_name, sep = ''))

#Load converted CSV
data_to_process <- read.csv(paste(file_directory, file_input_name, sep = ''))
#data_to_process$agegroup <- paste(data_to_process$agegroup, "_", sep='')
#data_to_process$agegroup <- paste(data_to_process$agegroup, data_to_process$region10, sep='')

##Change column names
names(data_to_process)[grep('age.*group', names(data_to_process), ignore.case = TRUE)] <- 'age_group'
names(data_to_process) <- gsub(paste('month(.*)', sep = ''), 'month\\1', names(data_to_process), ignore.case = TRUE)
names(data_to_process)[grep('date.*in', names(data_to_process), ignore.case = TRUE)] <- 'date_in'
names(data_to_process) <- gsub(paste('ach_', country, '(.*)', sep = ''), 'ach\\1', names(data_to_process), ignore.case = TRUE)
names(data_to_process)[grep('population', names(data_to_process), ignore.case = TRUE)] <- 'population'
names(data_to_process) <- gsub('(.*)SY', '\\1SY', names(data_to_process), ignore.case = TRUE)
if (country == 'Ecuador') {
  names(data_to_process) <- gsub('ach_NO_FO', 'ach', names(data_to_process), ignore.case = TRUE)
}
##Format dates
###Delete empty dates
data_to_process <- data_to_process[data_to_process$date_in != '',]

data_to_process$month<-as.numeric(as.character(substr(data_to_process$date_in,4,5)))
data_to_process$year<-as.numeric(as.character(substr(data_to_process$date_in,1,2)))+2000
date_char<-as.Date(paste(data_to_process$year,data_to_process$month,rep("01", nrow(data_to_process) )),format='%Y %m %d')
data_to_process$date<-date_char
  
###Get desired start date
start_year <- data_to_process$date[1] #Default; otherwise set manually using the form: '2001-01-01'
###Keep dates after start date, inclusive
data_to_process <- data_to_process[which(as.numeric(format(as.Date(data_to_process$date), '%Y')) >= as.numeric(format(as.Date(start_year), '%Y'))),]
##Add new variable columns

  #list.A10.B99 <- c(7,8,9,10,11,12,13,14,15,16,22,23,24,25,26)
  #data_to_process$A10_B99 = rowSums(data_to_process[,list.A10.B99])
  data_to_process$A10_B99_nopneumo = data_to_process$A10_B99 - data_to_process$B95 - data_to_process$A403
  data_to_process$ach_noj = data_to_process$ACH_Chile - data_to_process$J00_99
  data_to_process$cJ20_J22 = data_to_process$J20 + data_to_process$J21 + data_to_process$J22
  data_to_process$FLU = data_to_process$J09 + data_to_process$J10 + data_to_process$J11

  
  sub.ds<-data_to_process[data_to_process$region10==1 & data_to_process$age_group==2,]
  plot(sub.ds$date, sub.ds$J12_18)


##AGGREGATE FLU ACROSS Age groups
#flu.agg<- aggregate(data_to_process$FLU, by=list(data_to_process$date), FUN=sum, na.rm=TRUE)
#names(flu.agg)<-c("date","flu.agg")
#data_to_process<- merge(data_to_process, flu.agg, by="date")


##Remove unwanted variable columns; take out diahrea codes from ICD10 bc of rota vax
if (country == 'Brazil' | country=="Chile" | country=="Ecuador" | country=="Mexico" | country=="Chile_state" | country=="Mexico_state") {
  keep <- c('age_group', 'region10', 'STATE', 'date', 'J12_18',  'A10_B99_nopneumo','A17', 'A18', 'A19','A39','A41', 'B20_24','B34', 'B96', 'B97', 'B99',   'C00_D48', 'D50_89', 'E00_99', 'E10_14', 'E40_46', 'G00_99_SY', 'H00_99_SY', 'I00_99','I60_64',  'cJ20_J22', 'K00_99', 'K35', 'K80', 'L00_99', 'M00_99', 'N00_99', 'N39', 'P00_99', 'P05_07', 'Q00_99', 'S00_T99','U00_99', 'V00_Y99', 'Z00_99',  'ach_noj' )
} else if (country == 'US') {
  keep <- c('age_group', 'date', 'm_ACP','X_001_139_SY', 'X_240_279', 'X_280_289', 'X_320_359_SY', 'X_360_389_SY', 'X_390_459', 'X_520_579', 'X_580_629', 'X_680_709', 'X_740_759', 'X_760_779', 'X_780_799', 'X_800_999', 'X_V00_V91', 'ach_sid_noresp', 'm_00861', 'm_0088', 'm_430_434', 'm_466', 'm_599', 'm_78791', 'X_710_739', 'm_540', 'X_630_679', 'm_009', 'm_574')
}


keep2<-match(keep,names(data_to_process))
keep2<-keep2[!is.na(keep2)]
data_to_process <- cbind(data_to_process[,keep2])

##Sum values with same age_group and date (eliminates regional and hospital type specificity)
#data_to_process <- aggregate(data_to_process[, 3:length(data_to_process)], by = list(age_group = data_to_process$age_group, date = data_to_process$date), FUN = sum)
##Reorder by age_group
data_to_process <- data_to_process[order(data_to_process$age_group),]
##Modify age groups
data_to_process <- data_to_process[!is.na(data_to_process$age_group),]
  data_to_process$date.merge <- paste(data_to_process$date, "_", sep='')
  data_to_process$date.merge <- paste(data_to_process$date.merge, data_to_process$region10, sep='')
  
  data_to_process2s <- subset(data_to_process, data_to_process$age_group==2)
  data_to_process9s <- subset(data_to_process, data_to_process$age_group==9)
  
  data_to_process92s_2 <- aggregate(. ~ date.merge + region10, rbind(data_to_process2s,data_to_process9s), sum)
  data_to_process92s_2$date <- data_to_process2s$date
  data_to_process92s_2$region10 <- data_to_process2s$region10
  
  #data_to_process[data_to_process$age_group == 9, (match('date.merge', names(data_to_process)) + 1):length(names(data_to_process))] <- data_to_process[data_to_process$age_group == 2, (match('date.merge', names(data_to_process)) + 1):length(names(data_to_process))] + data_to_process[data_to_process$age_group == 9, (match('date.merge', names(data_to_process)) + 1):length(names(data_to_process))]
  data_to_process92s_2$age_group[data_to_process92s_2$age_group == 11] <- 92
  data_to_process <- data_to_process[data_to_process$age_group != 2,]
  data_to_process <- data_to_process[data_to_process$age_group != 9,]
  data_to_process <- rbind(data_to_process, data_to_process92s_2)
  
  data_to_process$age_group <- paste(data_to_process$age_group, "_", sep='')
  data_to_process$age_group <- paste(data_to_process$age_group, data_to_process$region10, sep='')

#data_to_process$age_group <- paste(data_to_process$age_group, "_", sep='')
#data_to_process$age_group <- paste(data_to_process$age_group, data_to_process$region10, sep='')

#data_to_process <- data_to_process[!is.na(data_to_process$age_group),]
#if (country == 'Mexico' || country == 'US') { #Weird thing for Mexico and US - data file is probably duplicated
#  data_to_process[, 3:ncol(data_to_process)] <- data_to_process[, 3:ncol(data_to_process)]/2
#}
##Remove data for age_group with too few counts
###Function to remove data with too few counts
#frequencyCheck <- Vectorize(function(age_group, variable) {
#  if (mean(data_to_process[variable][data_to_process$age_group == age_group,]) < 10) {
#    data_to_process[variable][data_to_process$age_group == age_group,] <<- rep(NA, times = length(data_to_process['age_group'][data_to_process$age_group == age_group,]))
#    return(FALSE)
#  } else {
#    return(TRUE)
#  }
#})
###Apply function to remove data with too few counts
#supress_output <- outer(unique(data_to_process$age_group), names(data_to_process)[(match('date', names(data_to_process)) + 1):length(names(data_to_process))], FUN = frequencyCheck)
##Reset row indices
#rownames(data_to_process) <- NULL

#Save processed data to file
  data_to_process<-data_to_process[substr(data_to_process$age_group,1,2)=='92',]
  data_to_process$region10<-NULL
  data_to_process$date.merge<-NULL
  
  write.csv(data_to_process,paste0('C:/Users/dmw63/Dropbox (Personal)/Meta analysis datasets/', country, '/',"prelog_Chile_state_u2.csv"), row.names = FALSE)

###END LOOP
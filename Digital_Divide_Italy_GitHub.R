#required package
#install.packages("poseticDataAnalysis")
library(poseticDataAnalysis)

#file containing the "POsetDR.plot" function
source("OBE_plot_v4.R")

#auxiliary function to built the frequency distribution of observed profiles
profiles_freq_dist_weighted<-function(matrix, weights){
  #matrix is the dataframe containing the list of observed profiles
  #weights is the vector of sample weights of each observed profiles
  k<-ncol(matrix)
  #Transforming profiles into integer numbers
  dec_profiles<-as.matrix(matrix)%*%(2^((k:1)-1))
  #building the frequency distribution of observed profiles (frequencies obtained from sample weights)
  observed_profiles_list<-sort(unique(dec_profiles))
  freq<-freq_W<-vector(length=length(observed_profiles_list))
  for(i in 1:length(observed_profiles_list)){
    freq_W[i]<-sum(weights*(dec_profiles==observed_profiles_list[i]))
    freq[i]<-sum(dec_profiles==observed_profiles_list[i])
  }
  obs_Profiles<-t(sapply(observed_profiles_list ,function(x){ as.integer(intToBits(x))})[k:1,])
  return(list(obs_Profiles,freq,freq_W))
}


#Loading original data
df<-read.table(file="Multiscopo-2023.txt", header=T, sep="\t")

#selecting the variables of interest
all<-c("PCTEMPO","FREQPC12",
       "INTTEMPO","FREQIN12","INCOMU5" ,"INCOMU6" ,"INTATT8","INTATT11","PCOPEWO" ,"PCOPE_FILE",
       "PCOPEX"  ,"PCOPEPH" ,"PCOPECO","PA_CO" ,"PA_BA"  ,"PA_SE"  ,"PA_DWN" ,"PA_PRE" ,"PA_RED",
       "PA_CER","PA_PPS" ,"PA_SC","PA_alt","INTCOM" ,"REGMf","ISTRMi","ETAMi","COEFIN")
Data<-df[,all]
#removing record with NA values
Data<-na.omit(Data)


#Binarization of non binary varaibles
Data$PCTEMPO<-(Data$PCTEMPO<=2)*1
Data$INTTEMPO<-(Data$INTTEMPO<=2)*1
Data$FREQPC12<-(Data$FREQPC12<=2)*1
Data$FREQIN12<-(Data$FREQIN12<=2)*1
Data$PA_RED<-(Data$PA_RED<=1)*1
Data$INTCOM<-(Data$INTCOM<=1)*1
Data$INCOMU5<-(Data$INCOMU5==2)*1
Data$INCOMU6<-(Data$INCOMU6==8)*1
Data$INTATT11<-(Data$INTATT11==2)*1
Data$INTATT8<-(Data$INTATT8==2)*1
Data$PCOPEWO<-(Data$PCOPEWO==2)*1
Data$PA_SE<-(Data$PA_SE==6)*1


#########################
#Final Variables selection
#########################
mantenute<-c("INCOMU5" ,"INCOMU6" ,"INTATT8","INTATT11","PA_SE","INTCOM","REGMf","ISTRMi","ETAMi","COEFIN")
Data3<-Data[, (names(Data) %in% mantenute)]
Data_def<-Data3[, !(names(Data3) %in% c("REGMf","ISTRMi","ETAMi","COEFIN"))]



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
######  GLOBAL ANALISYS
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#defining the weights vector
weights<-Data[,"COEFIN"]/sum(Data[,"COEFIN"])

B<-profiles_freq_dist_weighted(Data_def,weights)

#computing the optimal map
res<-OptimalBidimensionalEmbedding(B[[1]], B[[3]])

#global loss value
res$bestLossValue

names<-c("Email","Snet","Info","Bank","Padm","Ecomm")
N<-10000
#visulization of optimal embedding
x11()
POsetDR.plot(res, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = T, Profiles_col = "black", 
             #varNames = c("Fr_PC", "Fr_Int", "E_Comm","PA_Inf","PA_Serv", "PC_Soft", "Int_Serv", "PA_Doc"),
             varNames = names,
             Profiles_expand=c(1,1), title = "Global Map", type= "bubbles_and_approximation", expandY = c(0.1,0.1),expandX = c(0.1,0.1),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = "viridis") #plasma, turbo


########################################################################################################################
########################################################################################################################
########################################################################################################################
######  Decomposed analisys: subdivision by age
########################################################################################################################
########################################################################################################################
########################################################################################################################

###############################################################
### age classes
###############################################################
E1<-which(Data3$ETAMi %in% c(5,6)) #14-17
E2<-which(Data3$ETAMi %in% c(7,8)) #18-24
E3<-which(Data3$ETAMi %in% c(9,10)) #25-44
E4<-which(Data3$ETAMi %in% c(11,12)) #45-59
E5<-which(Data3$ETAMi %in%c(13,14)) #60-74
E6<-which(Data3$ETAMi %in%c(15)) #75+

nomi1<-c("Global loss", "Profile's Appr. error","Units' err. dist.","Prof.'s Appr. err.or vs Pop. shares ", "Map")

#14-15
eta<-"from 14 to 17 years"
B<-profiles_freq_dist_weighted(Data_def[E1,],weights[E1]/sum(weights[E1]))
resE1<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#18-24
eta<-"from 18 to 24 years"
B<-profiles_freq_dist_weighted(Data_def[E2,],weights[E2]/sum(weights[E2]))
resE2<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#25-44 
eta<-"from 25 to 44 years" 
B<-profiles_freq_dist_weighted(Data_def[E3,],weights[E3]/sum(weights[E3]))
resE3<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#45-59
eta<-"from 45 to 59 years"
B<-profiles_freq_dist_weighted(Data_def[E4,],weights[E4]/sum(weights[E4]))
resE4<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#60-74
eta<-"from 60 to 74 years"
B<-profiles_freq_dist_weighted(Data_def[E5,],weights[E5]/sum(weights[E5]))
resE5<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#75 e più
eta<-"over 75 years"
B<-profiles_freq_dist_weighted(Data_def[E6,],weights[E6]/sum(weights[E6]))
resE6<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

###############################################################
### plotting the maps for all age classes
###############################################################
titoli<-c("14-17","18-24","25-44","45-59","60-74","75+")
rampa<-"viridis"
par1<-par()
x11()
par(mfrow = c(3, 2), mai=c(0.1,0.1,0.1,0.1))
POsetDR.plot(resE1, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames = names,
             Profiles_expand=c(1,1), title = titoli[1], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resE2, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames = names,
             Profiles_expand=c(1,1), title = titoli[2], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resE3, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames = names,
             Profiles_expand=c(1,1), title = titoli[3], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resE4, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames = names,
             Profiles_expand=c(1,1), title = titoli[4], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resE5, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames = names,
             Profiles_expand=c(1,1), title = titoli[5], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resE6, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames = names,
             Profiles_expand=c(1,1), title = titoli[6], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
par<-par1


#########################################################################################################################
########################################################################################################################
########################################################################################################################
######  Decomposed analisys: subdivision by school degree
########################################################################################################################
########################################################################################################################
########################################################################################################################


T1<-which(Data3$ISTRMi %in% c(1)) #Bachelor's or Master's degree
T2<-which(Data3$ISTRMi %in%c(7)) #High-school degree
T3<-which(Data3$ISTRMi %in%c(9)) #Middle-school degree
T4<-which(Data3$ISTRMi %in%c(10)) #Elementary-school degree or nothing

#Bachelor's or Master's degree
B<-profiles_freq_dist_weighted(Data_def[T1,],weights[T1]/sum(weights[T1]))
resT1<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#High-school degree
B<-profiles_freq_dist_weighted(Data_def[T2,],weights[T2]/sum(weights[T2]))
resT2<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#Middle-school degree
B<-profiles_freq_dist_weighted(Data_def[T3,],weights[T3]/sum(weights[T3]))
resT3<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

#Elementary-school degree or nothing
B<-profiles_freq_dist_weighted(Data_def[T4,],weights[T4]/sum(weights[T4]))
resT4<-BidimentionalPosetRepresentation(as.matrix(B[[1]]), B[[3]], res$bestVariablesPriority)

###############################################################
### plotting the maps for all degrees
###############################################################
titoli<-c("Bachelor's or Master's degree","High-school degree","Middle-school degree","Elementary-school degree or nothing")
rampa<-"viridis"
x11()
par(mfrow = c(2, 2), mai=c(0.2,0.2,0.2,0.2))
POsetDR.plot(resT1, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames =  names,
             Profiles_expand=c(1,1), title = titoli[1], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resT2, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames =  names,
             Profiles_expand=c(1,1), title = titoli[2], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resT3, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames =  names,
             Profiles_expand=c(1,1), title = titoli[3], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
POsetDR.plot(resT4, N=N, k=6, Grid = T, Grid_col = "darkgray", Profile = F, Profiles_col = "black", 
             varNames =  names,
             Profiles_expand=c(1,1), title = titoli[4], type= "bubbles", expandY = c(9.5,9.5),expandX = c(9.5,9.5),
             points_col="black", center_col= "red", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 30.5, rampPalette = rampa) 
par<-par1

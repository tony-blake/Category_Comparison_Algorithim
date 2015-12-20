# algorithim for pairing indivuduals from generated data based on simliarity of chromosomes
# Created by Tony Blake
# 14/12/15

############ Install Package and load dependencies  #################################################

source('http://www.bioconductor.org/biocLite.R')
biocLite('gplots')
biocLite('stats')
biocLite('plyr')
biocLite('multtest')
biocLite('xlsx')
biocLite('bioDist')
library(gplots)
library(stats)
library(plyr)
library(multtest)
library(xlsx)
library(reshape2)
library(bioDist)





#####################################################################################################
#############################    Generate data and dataframe  #######################################
####################################################################################################

a=c(20:40)
b=c("m","f")
c=c("Af", "As", "Eu")
c_names <- c("Age", "Gender", "Chromo1", "Chromo2", "Chromo3", "Chromo4", "Chromo5", "Chromo6", "Chromo7", "Chromo8")
df <- data.frame()
df = data.frame(matrix( ncol = 10, nrow = 1000))
colnames(df) <- c_names

for(i in 1:1000){
  
  for(j in 3:10){
    
    df[i,1] <- sample(a, size=1, replace=T)
    
    df[i,2] <- sample(b, size=1, replace=T)
    
    df[i,j] <- sample(c, size=1, replace=T)
    
  }
}

rownames(df) <- paste("Human", rownames(df), sep = "_")

#####################################################################################################  
################### Create table for all matched pairs ##############################################
####################################################################################################

dff <- df[,1]
dff <- data.frame(dff)
test4 <- data.matrix(dist(dff)) #create distance matrix
CM <- test4
CM[upper.tri(CM, diag = TRUE)] <- NA # upper tri and diag set to NA
chromopairs <- subset(melt(CM, na.rm = TRUE), value <= 5)   # melt and subset as before

####################################################################################################
######################### create function to calculate score of matched rows ######################
###################################################################################################

list1 <- list()
list2 <- list()
biglist <- list()
score <- function(i,j){
  
  
  list1[[i]] <- df[i, 3:10]
  list2[[j]] <- df[j, 3:10]
  
  
  biglist[[1]] <- list1
  biglist[[2]] <- list2
  
  biglist[[1]] <- unlist(biglist[[1]])
  biglist[[2]] <- unlist(biglist[[2]])
  
  
  
  length(which(biglist[[1]]==biglist[[2]]))
  
}





######################################################################################################
############### make matrix  for pairwise comparisons based on score  ####################################################################### 
######################################################################################################


# make 1000 x 1000 matrix of where age difference are greater than 5
CA <- test4
for(i in 1:1000){
  
  CA[,i][CA[,i]<=5] <- 0
  
}

colnames(CA) <- paste("Human", colnames(CA), sep = "_")
rownames(CA) <- paste("Human", rownames(CA), sep = "_")

# make  1000 x 1000 matrix of all scores between pairs
ar <- array()
ar <- array(ar, c(1000,1000), dimnames=list(rownames(CA), rownames(CA)))
for(i in 1:1000) {
  for(j in 1:1000) {
    
    # intersect / union
    ar[i,j] <- score(i,j)/8
    
  }
}


# Add the two 1000 x 1000 matrices together and remove all entries greater than 1 
ARC <- ar + CA

for(i in 1:1000){
  
  ARC[,i][ARC[,i] > 1] <- -5/8
  
}




#####################################################################################################
################################# Find  unique pairs ###############################

#convert matrix of scores back to unnormalised form

CC <-  ARC
CC[upper.tri(CC, diag = TRUE)] <- NA 

CCC <- CC*8 

for(i in 1:1000){
  
  CCC[,i][CCC[,i]==8] <- 10
  
}

#for(i in 1:1000){
  
  #CCC[,i][CCC[,i]==0] <- -5
  
#}

#######################################################################################################################
################################# collect all pairs for 10, 7, 6, 5, 4, 3, 2, 1, -5 and exclude repeats ###############
#######################################################################################################################

tens <- subset(melt(CCC, na.rm = TRUE), value == 10) #find all pairs with score 10
tens <- tens[!duplicated(tens[,1]),] #remove duplicates in column 1
tens <- tens[!duplicated(tens[,2]),] #remove duplicates in column 2

#dataframe for sevens
sevens <- subset(melt(CCC, na.rm = TRUE), value == 7) #find all pairs with score 7
sevens <- sevens[!duplicated(sevens[,2]),] #remove duplicates in column 1
sevens <- sevens[!duplicated(sevens[,1]),] #remove duplicates in column 2
exclude<-!((sevens$Var2%in%tens$Var2)) # excluding individuals already included in "tens" that are also in sevens 
sven <- sevens[exclude,] #removing individuals from seven that are included in ten
exclude2 <- !((sven$Var2%in%sven$Var1)) # excluding individuals already included in column1 of sven that are also in column2 
sven <- sven[exclude2,] #removing individuals from column1 of sven that are included in column 2


# combine dataframes for tens and sevens
super <- rbind(tens, sven) 
super <- super[!duplicated(super[,1]),] #remove duplicates 
exclude3 <- !((super$Var1%in%super$Var2)) #excluding individuals already included in column 2 of super  that are also in column 1 
super <- super[exclude3,] # remove excluded individuals
intv <- intersect(super$Var1, super$Var2) #check to members of pairs are not repeated


# datframe for sixs
sixs <- subset(melt(CCC, na.rm = TRUE), value == 6) #find all pairs with score 7
sixs <- sixs[!duplicated(sixs[,1]),]  #remove duplicates in column 1
sixs <- sixs[!duplicated(sixs[,2]),]   #remove duplicates in column 1
exclude4 <-!((sixs$Var2%in%super$Var2)) # excluding individuals already included in "super" that are also in "sixs" 
six <- sixs[exclude4,] #removing excluded indivduals


#combine  dataframes for sixs with tens and sevens
supermax <- rbind(super, six) 
supermax <- supermax[!duplicated(supermax[,1]),] #remove duplicated individuals 
supermax <- supermax[!duplicated(supermax[,2]),] # remove duplicated individuals
exclude5 <- !((supermax$Var1%in%supermax$Var2)) # exclude individulas already included in column 2 of supermax that are also in column1
supermax <- supermax[exclude5,]


#need to check other way for  exclude5 <- !((supermax$Var2%in%supermax$Var1)) 

#dataframe for fives
fives <- subset(melt(CCC, na.rm = TRUE), value == 5) #comments follow same format as above with relevant index and name 
fives <- fives[!duplicated(fives[,1]),]
fives <- fives[!duplicated(fives[,2]),]
exclude6 <-!((fives$Var2%in%supermax$Var2))
fives <- fives[exclude6,]

#combine dataframe for fives with sixs, sevens and tens
minimax <- rbind(supermax, fives)
minimax <- minimax[!duplicated(minimax[,1]),]
minimax <- minimax[!duplicated(minimax[,2]),]
exclude5 <- !((minimax$Var1%in%minimax$Var2))
minimax <- minimax[exclude5,]

#datafame for fours
fours <- subset(melt(CCC, na.rm = TRUE), value == 4)
fours <- fours[!duplicated(fours[,1]),]
fours <- fours[!duplicated(fours[,2]),]
exclude8 <-!((fours$Var2%in%minimax$Var2))
fours <- fours[exclude8,]

#combine dataframe for fours with fives, sixs, sevens and tens
maxmax <- rbind(minimax, fours)
maxmax <- maxmax[!duplicated(maxmax[,1]),]
maxmax <- maxmax[!duplicated(maxmax[,2]),]
exclude9 <- !((maxmax$Var1%in%maxmax$Var2))
maxmax <- maxmax[exclude9,]

#dataframe for threes
trees <- subset(melt(CCC, na.rm = TRUE), value == 3)
trees <- trees[!duplicated(trees[,1]),]
trees <- trees[!duplicated(trees[,2]),]
exclude10 <-!((trees$Var2%in%maxmax$Var2))
trees <- trees[exclude10,]

#combine dataframe for threes withs fours, fives, sixs, sevens and tens 
verymax <- rbind(maxmax, trees)
verymax <- verymax[!duplicated(verymax[,1]),]
verymax <- verymax[!duplicated(verymax[,2]),]
exclude11 <- !((verymax$Var1%in%verymax$Var2))
verymax <- verymax[exclude11,]


#dataframe for twos
twos <- subset(melt(CCC, na.rm = TRUE), value == 2)
twos <- twos[!duplicated(twos[,1]),]
twos <- twos[!duplicated(twos[,2]),]
exclude12 <-!((twos$Var2%in%verymax$Var2))
twos <- twos[exclude12,]

#combine dataframe for twos with threes, fours, fives, six
moremax <- rbind(verymax, twos)
moremax <- moremax[!duplicated(moremax[,1]),]
moremax <- moremax[!duplicated(moremax[,2]),]
exclude13 <- !((moremax$Var1%in%moremax$Var2))
moremax <- moremax[exclude13,]


#dataframe for ones 
ones <- subset(melt(CCC, na.rm = TRUE), value == 1)
ones <- ones[!duplicated(ones[,1]),]
ones <- ones[!duplicated(ones[,2]),]
exclude14 <-!((ones$Var2%in%moremax$Var2))
ones <- ones[exclude14,]

#combine dataframe for ones with twos, threes, fours, fives, sixs, sevens, tens
almostmax <- rbind(moremax, ones)
almostmax <- almostmax[!duplicated(almostmax[,1]),]
almostmax <- almostmax[!duplicated(almostmax[,2]),]
exclude15 <- !((almostmax$Var1%in%almostmax$Var2))
almostmax <- almostmax[exclude15,]

#dataframe for ones 
zeros <- subset(melt(CCC, na.rm = TRUE), value == 0)
zeros <- zeros[!duplicated(zeros[,1]),]
zeros <- zeros[!duplicated(zeros[,2]),]
exclude16 <-!((zeros$Var2%in%almostmax$Var2))
zeros <- zeros[exclude16,]

#combine dataframe for ones with twos, threes, fours, fives, sixs, sevens, tens
penultimatemax <- rbind(moremax, zeros)
penultimatemax <- penultimatemax[!duplicated(penultimatemax[,1]),]
penultimatemax <- penultimatemax[!duplicated(penultimatemax[,2]),]
exclude17 <- !((penultimatemax$Var1%in%penultimatemax$Var2))
penultimatemax <- penultimatemax[exclude17,]


#dataframe for zeros








#dataframe for minus fives
minus <- subset(melt(CCC, na.rm = TRUE), value == -5)
minus <- minus[!duplicated(minus[,1]),]
minus <- minus[!duplicated(minus[,2]),]
exclude18 <-!((minus$Var2%in%penultimatemax$Var2))
minus <- minus[exclude18,]
exclude18a <-!((minus$Var1%in%penultimatemax$Var1))
minus <- minus[exclude18a,]



#combine dataframe for minus with ones, twos, threes, fours, fives, sixs, sevens, tens
max <- rbind(penultimatemax, minus)
max <- max[!duplicated(max[,1]),]
max <- max[!duplicated(max[,2]),]
exclude19 <- !((max$Var1%in%max$Var2))
max <- max[exclude19,]


intv <- intersect(max$Var1, max$Var2)

# create remainder dataset
N <- 500 - length(rownames(max))
rem <- data.frame()
rem = data.frame(matrix( ncol = 3, nrow = N))

#Populate dataset
col_names <- c("Var1", "Var2", "value")
colnames(rem) <- col_names
for(i in 1:N){
  rem[i,3] <- -5
}

for(i in 1:N){
  rem[i,2] <-  "_"
}

for(i in 1:N){
  rem[i,1] <-  "_"
}





total <- rbind(max, rem)




hist(total$value)
sum(total$value)
median(total$value)
mean(total$value)

write.xlsx(max, file="score_distrib.xlsx")
write.xlsx(df, file="simulated_data.xlsx")


########################################  end  ####################################################

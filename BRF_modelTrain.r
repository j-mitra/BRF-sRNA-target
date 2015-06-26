#Author: Joydeep Mitra
#Email: jmitra@uga.edu
#Last modified: 1st June, 2015

#Load randomForest library
library(randomForest)


#Read in Training data
positives.train<-read.table("./posTrain.txt",header=T,sep="\t")
negatives.train<-read.table("./negTrain.txt",header=T,sep="\t")

#Tune the mtry
tuned<-tuneRF(rbind(positives.train[,-ncol(positives.train)],negatives.train[,-ncol(positives.train)]), as.factor(c(positives.train$labels,negatives.train$labels)),plot=F,trace=F)
tuned.mtry<-tuned[,1][which.min(tuned[,2])][[1]]


#Set the number of samples to draw from both classes at 
#90% of positive class size
no.pos.samps<-round(0.9*nrow(positives.train))

#Set the number of trees to use per sampling
trees.per.forest<-5	
no.forests<-1000

#Create training labels for training
train.labels<-as.factor(c(rep(1,times=no.pos.samps),
						rep(-1,times=no.pos.samps)))

#Draw equal number of random samples from both classes
new.pos.sample<-sample(1:nrow(positives.train),no.pos.samps)
new.neg.sample<-sample(1:nrow(negatives.train),no.pos.samps)

#Initiate the BRF with mini-RF from sampled data
balanced.forest<-randomForest(rbind(positives.train[new.pos.sample,], 
				 negatives.train[new.neg.sample,]),train.labels,
				 ntree=trees.per.forest,mtry=tuned.mtry,importance=T) 

#The out-of-bag confusion matrix is used for evaluation on training set
oob.confusion<-balanced.forest$confusion[,1:2]


for(n in 2:no.forests)
{	
	#Draw new samples at each iteration
	new.pos.sample<-sample(1:nrow(positives.train),no.pos.samps)
	new.neg.sample<-sample(1:nrow(negatives.train),no.pos.samps)
	
	#Construct mini-RF to add to the BRF model
	add.forest<-randomForest(rbind(positives.train[new.pos.sample,], 
				negatives.train[new.neg.sample,]),train.labels,
				ntree=trees.per.forest,mtry=tuned.mtry,importance=T) 

	balanced.forest<-combine(balanced.forest,add.forest)

	#Cumulate out-of-bag confusion matrices
	oob.confusion<-oob.confusion + add.forest$confusion[,1:2]

}

#Calculate the OOB performance measures
oob.accuracy<-(oob.confusion[1]+oob.confusion[4])/sum(oob.confusion)
oob.posclass.acc<-oob.confusion[4]/(oob.confusion[2]+oob.confusion[4])
oob.negclass.acc<-oob.confusion[1]/(oob.confusion[1]+oob.confusion[3])
mcc.nume<-(oob.confusion[4]*oob.confusion[1]-oob.confusion[3]*oob.confusion[2])
mcc.denome_a<-(oob.confusion[4]+oob.confusion[3])*(oob.confusion[2]+oob.confusion[4])
mcc.denome_b<-(oob.confusion[1]+oob.confusion[3])*(oob.confusion[1]+oob.confusion[2])
mcc.denome<-sqrt(mcc.denome_a*mcc.denome_b)
oob.mcc<-mcc.nume/mcc.denome

#Consolidate results in a single variable
oob.results<-c("Accuracy"=oob.accuracy,"Sensitivity"=oob.posclass.acc,
			 "Specificity"=oob.negclass.acc, "MCC"=oob.mcc)

#Print the results for viewing
print(oob.results)

#Save the model for later use
save(balanced.forest,"./BRF_model.Rdata")
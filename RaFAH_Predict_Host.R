library(ranger)

args = commandArgs(trailingOnly=TRUE)

model_file<-args[1]
TestFileTable<-args[2]
out_name<-args[3]
threads<-as.numeric(args[4])


print(paste("Loading Model from ",model_file))
load(model_file)

print(paste("Reading input file ",TestFileTable))
TestData<-read.table(file=TestFileTable,sep="\t",header=TRUE,quote="",comment="")

TestSet<-TestData

gc()

print(paste("Passing data to Random Forest using ",threads," threads"))
TestDataPred <- predict(ranger_model_3, TestSet, type = "response", num.threads=threads)

TestDataPredResults<-rbind()
for (i in 1:nrow(TestDataPred$predictions)) {
	Winner_Score<-0
	Predicted_Host<-NA
	for (j in 1:ncol(TestDataPred$predictions)) {
		Taxon<-colnames(TestDataPred$predictions)[j]
		Score<-TestDataPred$predictions[i,j]
		if (Score > Winner_Score) {
			Winner_Score<-Score
			Predicted_Host<-Taxon
		}
	}
	TestDataPredResults<-rbind(TestDataPredResults,c(as.vector(TestData[i,1]),Predicted_Host,Winner_Score))
}

colnames(TestDataPredResults)<-c("Sequence","Predicted_Host","Winner_Score")
TestDataPredResults<-cbind(TestDataPredResults,TestDataPred$predictions)
print(paste("Writing predictions to ",out_name))
write.table(TestDataPredResults,file=out_name,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#
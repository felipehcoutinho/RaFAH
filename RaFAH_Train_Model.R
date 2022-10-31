#This code builds random forest models for the predictions of hosts of prokaryote viruses based on a given set of predictors.
#In this case the predictors are computed hmmsearch scores of genomic proteins versus a custom database of HMM models 
#It should also work with any other set of numerical predictors (e.g. number of occurrence of pFAM od pVOG domains in genomes)
args = commandArgs(trailingOnly=TRUE)

out_model_file<-args[1]
TrainFileTable<-args[2]
threads<-as.numeric(args[3]) #Number of threads to use when growing trees
confusion_file<-args[4]
importance_file<-args[5]


library(ranger)
options(expressions = 5e5)

#out_model_file<-"Test.RData"
#TrainFileTable<-"RaFAH_Genome_to_Domain_Score_Min_Score_50-Max_evalue_1e-05_Training.tsv"
#threads<-72
#confusion_file<-"Test_Confusion.tsv"
#importance_file<-"Test_Importancen.tsv"

tree_num<-1000 #Number of trees to grow
var_tries<-5000 #Variables to try at each split when growing trees
node_size<-1 #Node size for the probabilistic forest

TrainSet<-read.table(file=TrainFileTable,sep="\t",header=TRUE,quote="",comment="")
TrainSet<-subset(TrainSet, select=-c(Variable))
TrainSet<-TrainSet[!is.na(TrainSet$Host),]
TrainSet$Host<-factor(as.character(TrainSet$Host))

var_count<-(ncol(TrainSet) - 1)
if (var_count < var_tries) {
	var_tries<-var_count
}

start_time <- Sys.time()
set.seed(100)
ranger_model_3<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=TRUE,importance="impurity",min.node.size=node_size,num.threads=threads)
end_time <- Sys.time()

end_time - start_time

save.image(file=out_model_file)

#Print table of confusion matrix
#write.table(ranger_model_3$confusion.matrix,file=confusion_file,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

#Print Global importance matrix
#write.table(ranger_model_3$variable.importance,file=importance_file,sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

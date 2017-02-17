FoldCrossValidation.V1.BayesB <- function(dataset, trait, genoID,M,nFolds, nRepeats){
  require(BGLR)  
  data<-dataset# [which(dataset[ ,genoID] %in% rownames(Klist[[1]])),]
  data<-data[!is.na(data[,trait]),]    
  # rownames(data)<-data$CLONE
  nInd <- dim(data)[1] 
  accuracies<-data.frame()
  blups<-data.frame()
  for (rep in 1:nRepeats){ 
    print(paste("Rep ",rep," of ",nRepeats,sep=""))
    folds <- sample(rep(1:nFolds, length.out=nInd))
    BLUPSthisRep<-data.frame()
    for (fold in 1:nFolds){
      print(paste("Fold ",fold," of ",nFolds,sep=""))
      indInFold <- which(folds == fold)
      indNotInFold <- which(folds != fold)
      ValidSet<-data[indInFold,genoID]
      TrainSet<-data[indNotInFold,genoID]
      data1 <- data[which(data$CLONE %in% union(TrainSet,ValidSet)),c("CLONE",trait,paste(trait,"ebv",sep="."),paste(trait,"wt",sep="."))]
      data1 <- data1[!is.na(data1[,trait]),]
      rownames(data1) <- data1$CLONE
      train = TrainSet[TrainSet %in% data1$CLONE]
      test = ValidSet[ValidSet%in% data1$CLONE]
      datatrain <- data1[train,]
      datatest <- data1[test,]
      M.trn=as.matrix(M[TrainSet,])
      Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=datatrain) 
      y = datatrain[,trait]
      X = model.matrix(~1,data=datatrain)
      weight<-sqrt(data1[train,paste(trait,"wt",sep=".")])
      Z = Z/weight
      X = X/weight
      y = y/weight
      ETA.BayesB = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BayesB"))
      fit.BayesB = BGLR(y=y,ETA=ETA.BayesB,nIter=10000,burnIn=1000,thin=5)
      pred.tst.BayesB=(as.matrix(M[ValidSet,])%*%as.matrix(fit.BayesB$ETA[[2]]$b))
      pred.tst.Bayes<-as.data.frame(pred.tst.BayesB)
      BLUPSthisFold<-data.frame(CLONE=rownames(pred.tst.BayesB), PRED=pred.tst.Bayes$V1)
      BLUPSthisRep<-rbind(BLUPSthisRep,BLUPSthisFold)
    } 
    BLUPSthisRep[,"Rep"]<-rep
    # Calc accuracy after predicting each fold to complete the rep
    BLUPSthisRep<-merge(BLUPSthisRep,data[,c("CLONE",paste(trait,".ebv",sep=""))],by="CLONE")
    accuracy.thisrep<-cor(BLUPSthisRep$PRED,BLUPSthisRep[,paste(trait,".ebv",sep="")], use="complete.obs")
    AccuracyThisRep<-data.frame(Trait=trait,Rep=rep,Accuracy=accuracy.thisrep)
    accuracies<-rbind(accuracies,AccuracyThisRep) 
    blups<-rbind(blups,BLUPSthisRep)
  }
  return(list(accuracies=accuracies,blups=blups))
}

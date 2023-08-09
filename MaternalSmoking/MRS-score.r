
source("/home/dengwq/pgpstore/Data-Wei/LassoSum/elasticNet_lassosum.R")
combined <- read.csv("maternal_smoking_score_weights_Joubert.csv")


inver_norm <- function(x){
    x <- qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))
    x <- (x-mean(x, na.rm=T))/sd(x, na.rm=T)
    x
    } 



##############
##############
##############
##############

generate_MS_score <- function(CPG_MAT = NULL, combined = NULL, Y_cont = NA)

commonX <- intersect(names(CPG_MAT), combined$CPG)

print(paste(length(commonX), "CpGs have summary statistics"))

Xin <- as.matrix(scale(CPG_MAT[,commonX]))

print(paste("Input data dimension is ", dim(Xin)) 

BETA_MAT <- combined[match(commonX, combined$CPG),c(1,2,3)]
rownames(BETA_MAT) <- BETA_MAT$CPG


#### elements of lasso sum: 

library("lassosum")

minP <- NA
nk <- 50
nr <- 21
pred_p <-list()

for (k in 1: nr){
	
lambda <- exp(seq(log(0.005), log(20), length.out= nk))
pred_p[[k]] <- NA
alpha <- seq(0,1,0.05)[k]

for (j in 1:nk){
	
	lam1 <- lambda[j]*alpha
	lam2 <- lambda[j]*(1-alpha)
	
elRes <- elnetR(lambda1 = lam1, lambda2=lam2, X = Xin, b = BETA_MAT$WEIGHT);
MScore <- scale(elRes$pred)

if (sum(elRes$pred==0, na.rm=T)==length(elRes$pred)){
pred_p[[k]][j] <- NA
	} else {
pred_p[[k]][j]  <- summary(lm(Y_cont ~ MScore))$coef[2,4]
}
}
minP[k] <- min(pred_p[[k]], na.rm=T)
}
minP


	alpha <- seq(0,1,0.05)[which.min(minP)]
	j = which.min(pred_p[[which.min(minP)]])
	lam1 <- lambda[j]*alpha
	lam2 <- lambda[j]*(1-alpha)
	print(c(alpha, lam1, lam2))

elRes <- elnetR(lambda1 =lam1, lambda2=lam2, X = Xin, b = BETA_MAT$WEIGHT);
MScore <- (elRes$pred)
MScore_target <- MScore

 
#### elements of lasso sum end 

combined_wt <- data.frame("CpGs" = BETA_MAT[,1], "WEIGHTS" = elRes$beta , "SCORE" = "MScore")

print(paste(sum(elRes$beta!=0), "CpGs are in the final score"))

print("A quick view of CpGs contributing to the final score")
print(head(BETA_MAT[elRes$beta!=0,]))

return(combined_wt)
}

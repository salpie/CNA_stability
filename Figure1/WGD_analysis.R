library(viridis)

##### margin


hist(dataf$probs_diff, breaks=10)

svmfit = svm(WGD ~ ., data=dataf, scale = F, kernel = "linear", cost = 5,probability=TRUE)

beta = t(svmfit$coefs)%*%as.matrix(dataf[svmfit$index,c(2,3)])


b0 = svmfit$rho
b1<-beta[[1]]
b2<-beta[[2]]

plot(dataf$noSegments[dataf$WGD==T],
     dataf$percGenomeAbberated[dataf$WGD==T],
     xlab="number segments", 
     ylab = "% genome altered",
     xlim=c(0,250),ylim=c(0,1),
     col="red",
     pch=17)
points(dataf$noSegments[dataf$WGD==F],
       dataf$percGenomeAbberated[dataf$WGD==F],
       col="black",
       pch=16)
legend(200,1,c("WGD","noWGD"),col=c("red","black"),pch=c(17,16))

abline(b=-0.00035,a=0.43)

intercept = b0/b1
gradient = -b2/b1


abline(intercept, gradient)
abline((b0 - 1) / b1, -b2 / b1, lty = 2)
abline((b0 + 1) / b1, -b2/ b1, lty = 2)

(b0 - 1) / b1
[1] 0.5382894
(b0 + 1) / b1
0.3262133

## probabilities

predictions <- predict(fit,dataf[,c(2,3)], probability=TRUE)
probs = attr(predictions,"probabilities")
probs = as.data.frame(probs)
dataf$mid = "FALSE"
dataf$mid[ probs[,1] <0.6 & probs[,1] > 0.4 | probs[,2] <0.6 & probs[,2] > 0.4] <- TRUE



plot(dataf$noSegments[dataf$WGD==T],
     dataf$percGenomeAbberated[dataf$WGD==T],
     xlab="number segments", 
     ylab = "% genome altered",
     xlim=c(0,250),ylim=c(0,1),
     col="red",
     pch=17)
points(dataf$noSegments[dataf$WGD==F],
       dataf$percGenomeAbberated[dataf$WGD==F],
       col="black",
       pch=16)
legend(200,1,c("WGD","noWGD"),col=c("red","black"),pch=c(17,16))

abline(b=-0.00035,a=0.43)
abline((b0 - 0.5) / b1, -b2 / b1, lty = 2)
abline((b0 + 0.5) / b1, -b2/ b1, lty = 2)

Close = dataf[,c(2,3)][dataf$mid==TRUE,]

symbols(x=c(Close$noSegments), y=c(Close$percGenomeAbberated), circles=rep(1:16), add=T, inches=F)


numbers = seq(0, 1, by=0.05)


for (i in 1:length(numbers)) {

dataf$WGD_class_lower = dataf$percGenomeAbberated >= (dataf$noSegments * +gradient) + (b0 + numbers[i]) / b1
dataf$WGD_class_upper = dataf$percGenomeAbberated <= (dataf$noSegments * +gradient) + (b0 - numbers[i]) / b1

predictionsOnMyData <- predict(svmfit,dataf[,c(2,3)])
predictionsOnMyData <- predictionsOnMyData[dataf$WGD_class_lower==dataf$WGD_class_upper]
Middles <- dataf[dataf$WGD_class_lower==dataf$WGD_class_upper, ]

data_pred <- cbind(Middles,predictionsOnMyData)
# What fraction of positives are actually true?
ppv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData)
print(ppv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData))

# What fraction of negatives are actually false?
npv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData)
print(npv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData))

if (i ==1) {
	final_row = cbind(numbers[i], ppv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData)[3], npv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData)[3])
}

	new_row = cbind(numbers[i], ppv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData)[3], npv(data_pred,truth=data_pred$WGD,estimate=data_pred$predictionsOnMyData)[3])
	final_row = rbind(final_row,new_row)

}

names(final_row) <- c("margin", "ppv", "npv")
write.table(final_row, file="ppv_npv_cointoss_WGD.txt", row.names=F, quote=F)

load("/Users/nowins01/Documents/Gland_project/Analysis/nSEG_PGA.RData")
boxplot(nSeg_PGA$pga[nSeg_PGA$type=="carcinomas"], dataf$percGenomeAbberated, names=c("stability_paper_CRCs", "ICGC"), col=c("orange","blue"))


plot(nSeg_PGA$numSeg,
     nSeg_PGA$pga,
     xlab="number segments", 
     ylab = "% genome altered",
     xlim=c(0,100),ylim=c(0,1),
     col="red",
     pch=17)

abline(b=-0.00035,a=0.43)
abline((b0 - 0.4) / b1, -b2 / b1, lty = 2)
abline((b0 + 0.4) / b1, -b2/ b1, lty = 2)



type = c("#ea7e5d", "#5deac5", "#5d83ea")
b <- ggplot(nSeg_PGA, aes(x = numSeg, y = pga)) + geom_point(aes(shape = type, color = type, size=3))+
scale_color_manual(values = type) + theme_cowplot(12) +
geom_abline(slope=-0.00035, intercept=0.43)+ geom_abline(slope=-b2 / b1, intercept=(b0 - 0.5) / b1)+
geom_abline(slope=-b2 / b1, intercept=(b0 + 0.5) / b1, colour=“red”)

#new intercept of 0.3792323


boxplot(nSeg_PGA$pga[nSeg_PGA$type=="ca-in-ads"], nSeg_PGA$pga[nSeg_PGA$type=="carcinomas"],dataf$percGenomeAbberated, names=c("stability_paper_ca-in-ads", "stability_paper_CRCs", "ICGC"), col=c("orange","purple","blue"), main="PGA")

b <- ggplot(nSeg_PGA, aes(x = numSeg, y = pga)) + geom_point(aes(shape = type, color = type, size=3))+
scale_color_manual(values = type) + theme_cowplot(12) +
geom_abline(slope=-0.00035, intercept=0.43)+ geom_abline(slope=-b2 / b1, intercept=(b0 - 0.5) / b1)+
geom_abline(slope=-b2 / b1, intercept=(b0 + 0.5) / b1)

intercept_low = 0.3792323
intercept_up = 0.4852704
slope = -0.00035


nSeg_PGA$WGD_class = nSeg_PGA$pga >= (nSeg_PGA$numSeg * +slope) + intercept_low
nSeg_PGA$WGD_class_other = nSeg_PGA$pga <= (nSeg_PGA$numSeg * +slope) + intercept_up
nSeg_PGA$WGD_class[nSeg_PGA$WGD_class == nSeg_PGA$WGD_class_other] <- "MIDDLE"


nSeg_PGA <- nSeg_PGA[order(nSeg_PGA$names),]
nSeg_PGA$names <- gsub("Utrecht_pT1s_.", "", nSeg_PGA$names, fixed = TRUE)

nSeg_PGA$patient <- sub("\\.[^.]*$", "", nSeg_PGA$names)

b <- ggplot(nSeg_PGA, aes(x = numSeg, y = pga)) + geom_point(aes(shape = type, color = type, size=3)) +
  scale_color_manual(values = type) + theme_cowplot(12) + geom_abline(slope=-b2 / b1, intercept=(b0 - 0.5) / b1, colour="red")+
geom_abline(slope=-b2 / b1, intercept=(b0 + 0.5) / b1, colour="red") 


sub = nSeg_PGA[,c(4,5)]
tab = table(sub)

tab = as.data.frame(tab)

plot = ggplot(tab, aes(fill=type, y=Freq, x=WGD_class)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    facet_wrap(~type) +
	theme_cowplot(12) +
    theme(legend.position="none") +
    xlab("")

#make new classification from data and remove "middle" group
nSeg_PGA$WGD_class_new = nSeg_PGA$WGD_class


#rules for deciding which sample is WGD
for (i in 1:length(unique(nSeg_PGA$patient))){
	sub = nSeg_PGA[nSeg_PGA$patient==unique(nSeg_PGA$patient)[i],]
	if (any(sub$WGD_class_new != "FALSE")) {
		if (any(sub$WGD_class_new == "TRUE") && any(sub$WGD_class_new == "MIDDLE")) {
			sub$WGD_class_new[sub$WGD_class_new=="MIDDLE"] <- "FALSE"
		}
		if (sum(sub$WGD_class_new == "MIDDLE") >= 2) {
			sub$WGD_class_new[sub$WGD_class_new=="MIDDLE"] <- "TRUE"
		}		
		if (any(sub$WGD_class_new == "FALSE") && any(sub$WGD_class_new == "MIDDLE")) {
			sub$WGD_class_new[sub$WGD_class_new=="MIDDLE"] <- "TRUE"
		}	
	}
	if (i == 1) {
		final_sub = sub
	} else {
		final_sub=rbind(final_sub,sub)
	}

}
## plot the new thing

b <- ggplot(final_sub, aes(x = numSeg, y = pga)) + geom_point(aes(shape = WGD_class_new, color = type, size=3))+
scale_color_manual(values = type) + theme_cowplot(12)+ geom_line(aes(group = patient))+ geom_line(aes(group = patient))+ geom_abline(slope=-b2 / b1, intercept=(b0 - 0.5) / b1)+
geom_abline(slope=-b2 / b1, intercept=(b0 + 0.5) / b1)

#now plot the fractions of types and samples WGD or not

sub = final_sub[,c(4,12)]
tab = table(sub)

tab = as.data.frame(tab)

plot = ggplot(tab, aes(fill=type, y=Freq, x=WGD_class_new)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    facet_wrap(~type) +
	theme_cowplot(12) +
    theme(legend.position="none") +
    xlab("")



##### redo same analysis as before:

nSeg_PGA=final_sub


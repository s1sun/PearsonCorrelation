#vcfdata=read.csv("C:\\Users\\sesun\\Documents\\pipeline\\somaticseq\\result\\mutstrvarSXY.txt")
vcfdata=read.delim("C:\\Users\\sesun\\Documents\\pipeline\\somaticseq\\result\\mutstrvarSXY.txt", head=TRUE, sep="\t", dec=".")
#remove first 5 columns
vcfdata=vcfdata[,-c(1:5)]

vcfname=names(vcfdata)
vcfname
#[1] "AC"  "AF"  "DP_N"  "DP_T" "ECNT" ... "VariantType" "set" "lable" "GroundT"
plot(AC~AF, data=vcfdata)
vcfmodel=lm(AC~AF, data=vcfdata)
vcfmodel
abline(vcfmodel,col="red")
summary(vcfmodel)

cor(vcfdata$AC, vcfdata$AF, use="complete")
cor(vcfdata[,1], vcfdata[,2], use="complete")
#Pcorr:0.06179868

cor(vcfdata$AF, vcfdata$ECNT, use="complete")
cor(vcfdata[,2], vcfdata[,3], use="complete")
#Pcorr:-0.0993995

dist=cor(vcfdata, use="complete")
myf=abs(dist)
myT=which(myf == max(myf[row(myf)!=col(myf)]), arr.ind = TRUE)
myT
#     	row col
#TQSI_NT	47   46
#TQSI		46   47

install.packages("corrplot")
library('corrplot')
corrplot.mixed(dist,tl.cex=0.5, number.cex=0.6) #lower.col="black", 
corrplot(dist, method="circle", tl.cex=0.7)

#Choose the circle with the largest diameter( corr:-0.98)

cc=as.formula(paste(colnames(vcfdata)[myT[2,2]], "~", paste(colnames(vcfdata)[myT[2,1]])))
maxlm=lm(cc, data=vcfdata)
maxlm
#Coefficients:
#(Intercept)         Var1  
#    3.083e-05     9.996e-01 
#equ: TQSI_NT=0.03268-1.06569(TQSI)

plot(cc, data=vcfdata, pch=20, type='p', las=1,  xlab=bquote( "X(" ~ .(vcfname[T[1,2]]) ~ ")" ),
         ylab=bquote( "Y(" ~ .(vcfname[T[2,2]]) ~ ")" ),main="Two most related columns")
abline(maxlm,col="red")
rp = vector('expression',2)
rp[1] = substitute(expression( Y == MYVALUE*X+VALUE), list( VALUE =        
                       format(maxlm$coefficients[1],dig=3), MYVALUE=format(maxlm$coefficients[2],dig=4) ) )[2]
pCORR=cor(vcfdata[T[2,2]], vcfdata[T[1,2]], use="complete")
rp[2] = substitute(expression(italic(Pcorr) == MYVALUE), list(MYVALUE = format(pCORR[1], digits =3)))[2]
legend('topright', legend = rp, bty = 'n')

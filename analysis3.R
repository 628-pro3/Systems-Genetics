library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(MASS)

geno <- read_csv("b6btbr_geno.csv")
pmap <- read_csv("b6btbr_pmap.csv", col_types=list(col_character(),col_character(),col_double()))
pheno <- read_csv("b6btbr_pheno.csv")
covar <- read_csv("b6btbr_covar.csv")

# sort the marker by chromosomes and positions
pmap_1_19 <- pmap %>% 
    filter(chr != "X") %>%
    mutate(chr=as.numeric(chr)) %>%
    group_by(chr) %>%
    arrange(pos)
pmap_X <- pmap %>% 
    filter(chr == "X") %>%
    arrange(pos)
pmap_sort <- rbind(pmap_1_19, pmap_X)

# sort the key phenotype
insulin <-  dplyr::select(pheno, MouseNum, log10_insulin_10wk) %>%
    arrange(desc(log10_insulin_10wk))

# reorder geno data sets
col_ind <- match(pmap_sort$marker, colnames(geno))
geno_col_sort <- geno[,c(1,col_ind)]

row_ind <- match(insulin$MouseNum,geno_col_sort$MouseNum)
geno_row_sort <- geno_col_sort[row_ind,]
add_sex <- covar[,c(1,2)]
data <- merge(geno_row_sort, add_sex)

dat <- stack(data[,-1])
head(dat)
colnames(dat)[2] <- "marker"
dat$marker <- as.character(dat$marker)
dat$insulin <- rep(insulin$log10_insulin_10wk,ncol(geno_row_sort))
geno_pheno <- left_join(dat,pmap)
write.csv(geno_pheno,"geno_pheno.csv",row.names=FALSE)
 
# plot
geno_pheno <- read.csv("geno_pheno.csv",colClasses = "character")
str(geno_pheno)
geno_pheno$insulin <- as.numeric(geno_pheno$insulin)
geno_pheno$pos <- as.numeric(geno_pheno$pos)
str(geno_pheno)

## Separate Chromosome 1-19 & Chromosome X
chrX <- geno_pheno %>%
    filter(marker != "Sex") %>%
    filter(chr=="X")
ggplot(chrX, aes(x=pos,y=insulin)) + 
    geom_point(aes(colour=values)) +
    xlab("position") +
    ylab("log10_insulin_10wk")
ggsave("chrX.png",width = 9, height = 9, dpi=80)

chr1_19 <- geno_pheno %>%
    filter(chr!="X" & marker!="Sex")
chr1_19$chr <- as.numeric(as.character(chr1_19$chr))
ggplot(chr1_19, aes(x=pos,y=insulin)) + 
    geom_point(aes(colour=values)) +
    facet_wrap(~chr) +
    xlab("position") +
    ylab("log10_insulin_10wk")
ggsave("chr1_19.png",width = 12, height = 12, dpi=80)


## Loop for ploting for each chromosome

for (i in 1:19) {
    chr_i <- chr1_19 %>%
        filter(chr == i)
    fn <- paste0("Chromosome ",i)
    ggplot(chr_i, aes(x=pos,y=insulin)) + 
        geom_point(aes(colour=values)) +
        facet_wrap(~chr) +
        xlab("position") +
        ylab("log10_insulin_10wk") +
        xlim(0,200) +
        ggtitle(fn)
    fn_address <- paste0("Chromosome_",i,".png")
    ggsave(fn_address,width = 9, height = 9, dpi=80)
}


## Modeling

# Chr1
chr1 <- chr1_19 %>%
    filter(chr == 1)
head(chr1)
str(chr1)

# AIC=-464.97, R2 = 0.2046  
chr1_marker1 <- chr1 %>%
    filter(insulin == min(insulin) & as.character(values)=="RR") # find the important markers
chr1_mod1 <- chr1 %>%
    filter(marker %in% chr1_marker1$marker)

chr1_mod1_wide <- chr1_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr1_mod1_wide)<-c(names(chr1_mod1_wide)[1], gsub("values.","",names(chr1_mod1_wide)[-1]))

na.ind <- unname(as.matrix(which(is.na(chr1_mod1_wide),arr.ind=TRUE)))
chr1_mod1_wide <- chr1_mod1_wide[-na.ind[1:(length(na.ind)/2)],]

fit1 <- stepAIC(lm(insulin~.,data=chr1_mod1_wide, na.action=na.exclude),direction="both")
summary(fit1)

r1 <- rstudent(fit1)
fitted1 <- fitted(fit1)
df1 = data.frame(residual=r1, fitted.value=fitted1)
ggplot(df1,aes(y=r1,x=fitted1)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)

## subset, R2=0.01696

fit1.1 <- lm(insulin~rs13475697+rs13475710+rs13475764+rs13475794+rs13475801+rs4222269+rs13476045+rs3685700,data=chr1_mod1_wide, na.action=na.exclude)
summary(fit1.1)

## subset, R2=0.01696

fit1.2 <- lm(insulin~rs13475801*rs13475794,data=chr1_mod1_wide, na.action=na.exclude)
summary(fit1.2)



# Chr4
chr4 <- chr1_19 %>%
    filter(chr == 4)
head(chr4)
str(chr4)

# AIC=-892.06, R2 = 0.01941, residual plot has patterns 
chr4_marker1 <- chr4 %>%
    filter(insulin == min(insulin) & pos > 80) # find the important markers
chr4_mod1 <- chr4 %>%
    filter(marker %in% chr4_marker1$marker)

chr4_mod1_wide <- chr4_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr4_mod1_wide)<-c(names(chr4_mod1_wide)[1], gsub("values.","",names(chr4_mod1_wide)[-1]))

na.ind <- unname(as.matrix(which(is.na(chr4_mod1_wide),arr.ind=TRUE)))
chr4_mod1_wide <- chr4_mod1_wide[-na.ind[1:(length(na.ind)/2)],]

fit2 <- stepAIC(lm(insulin~.,data=chr4_mod1_wide, na.action=na.exclude),direction="both")
summary(fit2)

r2 <- rstudent(fit2)
fitted2 <- fitted(fit2)
df2 = data.frame(residual=r2, fitted.value=fitted2)
ggplot(df2,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)



# Chr8
chr8 <- chr1_19 %>%
    filter(chr == 8)
head(chr8)
str(chr8)

# AIC=-866.8, R2 = 0.02164, residual plot has patterns
chr8_marker1 <- chr8 %>%
    filter(insulin == min(insulin) & as.character(values) == "RR") # find the important markers
chr8_mod1 <- chr8 %>%
    filter(marker %in% chr8_marker1$marker)

chr8_mod1_wide <- chr8_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr8_mod1_wide)<-c(names(chr8_mod1_wide)[1], gsub("values.","",names(chr8_mod1_wide)[-1]))

na.ind <- unname(as.matrix(which(is.na(chr8_mod1_wide),arr.ind=TRUE)))
chr8_mod1_wide <- chr8_mod1_wide[-na.ind[1:(length(na.ind)/2)],]

fit3 <- stepAIC(lm(insulin~.,data=chr8_mod1_wide, na.action=na.exclude),direction="backward")
summary(fit3)

r3 <- rstudent(fit3)
fitted3 <- fitted(fit3)
df3 = data.frame(residual=r3, fitted.value=fitted3)
ggplot(df3,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)


# Chr9
chr9 <- chr1_19 %>%
    filter(chr == 9)
head(chr9)
str(chr9)

# AIC=-912.04, R2 = 0.009757, residual plot has patterns
chr9_marker1 <- chr9 %>%
    filter(insulin <0 & as.character(values) == "BR") %>%
    filter(insulin == max(insulin))# find the important markers
chr9_mod1 <- chr9 %>%
    filter(marker %in% chr9_marker1$marker)

chr9_mod1_wide <- chr9_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr9_mod1_wide)<-c(names(chr9_mod1_wide)[1], gsub("values.","",names(chr9_mod1_wide)[-1]))

na.ind <- unname(as.matrix(which(is.na(chr9_mod1_wide),arr.ind=TRUE)))
chr9_mod1_wide <- chr9_mod1_wide[-na.ind[1:(length(na.ind)/2)],]

fit4 <- stepAIC(lm(insulin~.,data=chr9_mod1_wide, na.action=na.exclude),direction="backward")
summary(fit4)

r4 <- rstudent(fit4)
fitted4 <- fitted(fit4)
df4 = data.frame(residual=r4, fitted.value=fitted4)
ggplot(df4,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)



# Chr13
chr13 <- chr1_19 %>%
    filter(chr == 13)
head(chr13)
str(chr13)

# AIC=-872.91, R2 = 0.006594, residual plot has patterns
chr13_marker1 <- chr13 %>%
    filter(insulin > 0) %>%
    filter(insulin == min(insulin) & as.character(values)=="BR")
chr13_marker2 <- chr13 %>%
    filter(marker %in% chr13_marker1$marker) %>%
    filter(insulin > 0.1) %>%
    filter(insulin == min(insulin) & as.character(values)=="BR")# find the important markers
chr13_mod1 <- chr13 %>%
    filter(marker %in% chr13_marker2$marker)

chr13_mod1_wide <- chr13_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr13_mod1_wide)<-c(names(chr13_mod1_wide)[1], gsub("values.","",names(chr13_mod1_wide)[-1]))

chr13_mod1_wide <- na.omit(chr13_mod1_wide)

fit5 <- stepAIC(lm(insulin~.,data=chr13_mod1_wide, na.action=na.exclude),direction="backward")
summary(fit5)

r5 <- rstudent(fit5)
fitted5 <- fitted(fit5)
df5 = data.frame(residual=r5, fitted.value=fitted5)
ggplot(df5,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)



# Chr14
chr14 <- chr1_19 %>%
    filter(chr == 14)
head(chr14)
str(chr14)

# not significant
chr14_marker1 <- chr14 %>%
    filter(insulin == max(insulin) & as.character(values)=="RR")# find the important markers
chr14_mod1 <- chr14 %>%
    filter(marker %in% chr14_marker1$marker)

chr14_mod1_wide <- chr14_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr14_mod1_wide)<-c(names(chr9_mod1_wide)[1], gsub("values.","",names(chr14_mod1_wide)[-1]))

na.ind <- unname(as.matrix(which(is.na(chr14_mod1_wide),arr.ind=TRUE)))
chr14_mod1_wide <- chr14_mod1_wide[-na.ind[1:(length(na.ind)/2)],]

fit6 <- stepAIC(lm(insulin~.,data=chr14_mod1_wide, na.action=na.exclude),direction="both")
summary(fit6)



# Chr16
chr16 <- chr1_19 %>%
    filter(chr == 16)
head(chr16)
str(chr16)

# duplicated name 'NA' in data frame using '.'?
chr16_marker1 <- chr16 %>%
    filter(insulin == min(insulin) & as.character(values) == "RR")# find the important markers
chr16_mod1 <- chr16 %>%
    filter(marker %in% chr16_marker1$marker)

chr16_mod1_wide <- chr16_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr16_mod1_wide)<-c(names(chr16_mod1_wide)[1], gsub("values.","",names(chr16_mod1_wide)[-1]))

chr16_mod1_wide <- na.omit(chr16_mod1_wide)

fit7 <- stepAIC(lm(insulin~.,data=chr16_mod1_wide),direction="backward")
summary(fit7)



# Chr18
chr18 <- chr1_19 %>%
    filter(chr == 18)
head(chr18)
str(chr18)

# problem
chr18_marker1 <- chr18 %>%
    filter(insulin == max(insulin) & pos < 75)# find the important markers
chr18_mod1 <- chr18 %>%
    filter(marker %in% chr18_marker1$marker)

chr18_mod1_wide <- chr18_mod1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr18_mod1_wide)<-c(names(chr18_mod1_wide)[1], gsub("values.","",names(chr18_mod1_wide)[-1]))

na.ind <- unname(as.matrix(which(is.na(chr18_mod1_wide),arr.ind=TRUE)))
chr18_mod1_wide <- chr18_mod1_wide[-na.ind[1:(length(na.ind)/2)],]

fit8 <- stepAIC(lm(insulin~.,data=chr18_mod1_wide, na.action=na.exclude),direction="backward")
summary(fit8)


# geno_pheno_wide
geno_pheno_wide <- geno_pheno %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(geno_pheno_wide)<-c(names(geno_pheno_wide)[1], gsub("values.","",names(geno_pheno_wide)[-1]))

write.csv(geno_pheno_wide,"geno_pheno_wide.csv",row.names=FALSE)
write.table(geno_pheno_wide,"geno_pheno_wide.txt",row.names=FALSE, col.names=FALSE, quote=FALSE)

geno_pheno_wide <- read.csv("geno_pheno_wide.csv", stringsAsFactors=FALSE)
str(geno_pheno_wide)

des <- paste0(1:length(colnames(geno_pheno_wide)), " ", colnames(geno_pheno_wide), " b")
write.table(des, "geno_pheno_wide_des.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


# Random Forest
library(randomForest)
rf <- randomForest(insulin~.,data=geno_pheno_wide,
                   importance=TRUE,na.action=na.omit)
print(rf)
varImpPlot(rf)
plot(rf)
PartialPlot

head(round(importance(rf,type=1),2))
str(round(importance(rf,type=1),2))

impVar <- attr(importance(rf,type=1),"dimnames")[[1]][order(importance(rf,type=1),decreasing=TRUE)]
write.csv(impVar,"impVarList.csv", row.names=FALSE)
inf_imp_marker <- pmap %>%
    filter(marker %in% impVar[1:8])
write.csv(inf_imp_marker,"inf_imp_marker.csv", row.names=FALSE)




## RF Modeling

# Chr1 R2=0.004484
chr1_wide <- chr1 %>%
    dplyr::select(insulin,values,marker) %>%
    reshape(timevar="marker", idvar="insulin", direction="wide")
names(chr1_wide)<-c(names(chr1_wide)[1], gsub("values.","",names(chr1_wide)[-1]))

# Check chr1

pmap_chr1 <- pmap %>%
    filter(chr == "1")
ind.marker.chr1 <- which(colnames(geno) %in% pmap_chr1$marker)
geno_chr1 <- geno[,c(1,ind.marker.chr1)]
pheno_2 <- pheno %>%
    dplyr::select(MouseNum,log10_insulin_10wk)
chr1_check <- merge(geno_chr1,pheno_2)
chr1_check <- chr1_check[,-1]
ind_not <- which(!(colnames(chr1_check) %in% colnames(chr1_wide)))
col_not <- colnames(chr1_check)[ind_not]

rf.fit1 <- lm(insulin~rs13475752+rs4222269+rs13475822+rs13475826+rs13475834,data=chr1_wide)
summary(rf.fit1)


rf.df1 = data.frame(residual=rstudent(rf.fit1), fitted.value=fitted(rf.fit1))
ggplot(rf.df1,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)


# Overall R2=0.001518
rf.fit2 <- lm(insulin~rs13475752+rs4222269+rs13475822+rs13475826+rs13475834+rs6234650+rs6244558+rs13483071,data=geno_pheno_wide)
summary(rf.fit2)

rf.df2 = data.frame(residual=rstudent(rf.fit2), fitted.value=fitted(rf.fit2))
ggplot(rf.df2,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)


# Overall+sex R2=0.00434
rf.fit3 <- lm(insulin~rs13475752+rs4222269+rs13475822+rs13475826+rs13475834+rs6234650+rs6244558+rs13483071+Sex,data=geno_pheno_wide)
summary(rf.fit3)

rf.df3 = data.frame(residual=rstudent(rf.fit3), fitted.value=fitted(rf.fit3))
ggplot(rf.df3,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)


# Overall+sex＋Chr1 R2=0.00434
rf.LM.fit4 <- lm(insulin~rs13475697+rs13475710+rs13475729+rs13475752+rs13475764+rs13475794+rs13475801+rs4222269+rs13475822+rs13475826+rs13475834+rs13475912+rs13476045+rs6234650+rs6244558+rs13483071+rs3685700+Sex,data=geno_pheno_wide)
summary(rf.LM.fit4)

rf.df4 = data.frame(residual=rstudent(rf.LM.fit4), fitted.value=fitted(rf.LM.fit4))
ggplot(rf.df3,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)



# Tree R2=0.00434
tree.fit5 <- lm(insulin~rs13478610+rs13478667+Sex,data=geno_pheno_wide)
summary(tree.fit5)

tree.df5 = data.frame(residual=rstudent(tree.fit5), fitted.value=fitted(tree.fit5))
ggplot(tree.df5,aes(y=residual,x=fitted.value)) +
    geom_point() +
    geom_hline(yintercept=0,colour="red",lty=2)





library(vegan); library(spdep)

##############################################################
### Pripravne kroky
##############################################################

# Nacti vsech sest promennych zdravotniho stavu lesa 
def <- readXL("F:/IFER/clanek_Beskydy/#Resubmission/revize/propojeni_uprava.xls", rownames=FALSE, header=TRUE, na="", sheet="def", stringsAsFactors=TRUE)
def.nestand <- readXL("F:/IFER/clanek_Beskydy/#Resubmission/revize/propojeni_uprava.xls", rownames=FALSE, header=TRUE, na="", sheet="def_unst", stringsAsFactors=TRUE)

ord <- rda(def[2:6], scaled=TRUE)
summary(ord) # Podil variability vysvetlene prvnimi dvema osami, komponenty
set.seed(1)
biplot (ord, display = 'species', scaling=0, xlim=c(-0,0.6), ylim=c(-0.4,0.8))


# Nacti CI4, PC1 a PC2 (obe komponenty jsem si vypocital v Excelu podle vystupu z predchozich funkci)
ind <- readXL("F:/IFER/clanek_Beskydy/#Resubmission/revize/propojeni_uprava.xls", rownames=FALSE, header=TRUE, na="", sheet="ind", stringsAsFactors=TRUE)

rcorr.adjust(ind[,c("PC1","PC2")], type="pearson", use="complete") # velmi silna korelace mezi CI4 a PC1
scatterplot(PC1~PC2, reg.line=FALSE, smooth=FALSE, spread=FALSE, boxplots=FALSE, span=0.5, ellipse=FALSE, levels=c(.5, .9), data=ind) # Zavislost graficky

###############################################################
### GLM na PC1 a PC2
###############################################################
###############################################################

Dataset <- readXL("F:/IFER/clanek_Beskydy/#Resubmission/revize/propojeni_uprava.xls", rownames=FALSE, header=TRUE, na="", sheet="Dataset", stringsAsFactors=TRUE) # Nactu si cely dataset environmentalnich promennych

###############################################################
Dataset <- local({
  .Z <- scale(Dataset[,c("AGE_5","GEOLOGY","N","S","SN","SR","T", "ALTITUDE_M")])
  within(Dataset, {
    ALTITUDE_M <- .Z[,8]
    T <- .Z[,7]
    SR <- .Z[,6]
    SN <- .Z[,5]
    S <- .Z[,4]
    N <- .Z[,3]
    GEOLOGY <- .Z[,2]
    AGE_5 <- .Z[,1] 
  })
})
###############################################################

cbind <- cbind(ind, Dataset) # propojeni tabulek pro nasledujici vypocty

### PC1

PC1_spatial <- lagsarlm(PC1 ~ ALTITUDE_M + SN + GEOLOGY + AGE_5, cbind, lw1); summary(PC1_spatial)


### PC2

PC2_spatial <- lagsarlm(PC2 ~ ALTITUDE_M + AGE_5, cbind, lw1); summary(PC2_spatial)


#####################################################################
### Vypocet prostorove autokorelace (Moranovo I) vysvetlovane promenne a residualu linearniho modelu
#####################################################################

souradnice <- readXL("F:/IFER/clanek_Beskydy/#Resubmission/revize/propojeni_uprava.xls", rownames=FALSE, header=TRUE, na="", sheet="souradnice", stringsAsFactors=TRUE) # Nactu si cely dataset environmentalnich promennych

library(ape)
### Vytvorim inverzni matici vzdalenosti, na jejiz diagonalu manualne doplnim nuly
dists <- as.matrix(dist(souradnice))
dists.inv <- 1/dists
diag(dists.inv) <- 0
dists.bin <- (dists > 0 & dists <= 85)
lw1 <- mat2listw(dists.inv, style="W")


### LM se zohlednenou prostorovou autokorelaci
PC1_lm_spat <- lm(PC1 ~ ALTITUDE_M + SN + GEOLOGY + AGE_5 + colSums(dists.inv*cbind$PC1), cbind); summary(PC1_lm_spat)
PC2_lm_spat <- lm(PC2 ~ ALTITUDE_M + SN + GEOLOGY + AGE_5 + colSums(dists.inv*cbind$PC2), cbind); summary(PC2_lm_spat)

### Vypocet residualu linearniho modelu
PC1_lm <- lm(PC1 ~ ALTITUDE_M + SN + GEOLOGY + AGE_5, cbind); summary(PC1_lm)
PC2_lm <- lm(PC2 ~ ALTITUDE_M + AGE_5, cbind); summary(PC2_lm)
residual <- cbind(data.frame(residuals(PC1_lm)), data.frame(residuals(PC2_lm)), data.frame(residuals(PC1_lm_spat)), data.frame(residuals(PC2_lm_spat)))

### Klasicke moranovo kriterium
Moran.I(ind$PC1, dists.inv)
Moran.I(ind$PC2, dists.inv)
Moran.I(residual$residuals.PC1_lm., dists.inv)
Moran.I(residual$residuals.PC2_lm., dists.inv)
Moran.I(residual$residuals.PC1_lm_spat., dists.inv)
Moran.I(residual$residuals.PC2_lm_spat., dists.inv)

### Max a prumerna vzdalenost od nejblizsiho souseda
library(FNN)
row.names(dists) <- (souradnice$IDPLOT) # Prepisu ID v tabulce distance
colMeans(knn.dist(dists))


#####################################################################
### Zjisteni hodnoty prumerne PC1 a PC2 nekolika nejblizsich sousedu
#####################################################################

distance <- read.table("clipboard", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) # Nactu jenom sloupce X_m a Y_m
id_plots <- read.table("clipboard", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) # Nactu ID ploch
row.names(distance) <- (id_plots$IDPLOT) # Prepisu ID v tabulce distance
 
knn <- get.knn(distance, 7)
# Rucne prekopirovano do Excelu

komplet <- read.table("clipboard", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) # Nactu vse z listu knn
pc <- read.table("clipboard", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) # Nactu pomocna ID a PC1 a PC2

names(pc)[c(1,2,3)] <- c("ID_nn6","PC1_nn6","PC2_nn6") # Postupne prejmenovavam promenne v datasetu "pc"

komplet <- merge(komplet, pc, by="ID_nn1")
# ... menim hodnotu by=, dokud nepripojim vechno
write.table(komplet, "C:/Users/JT_2/Desktop/ifer/clanek/Beskydy/oprava/komplet_knn.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA") # Ulozim si to

#####################################################################
### Vykresleni histogramu
#####################################################################

layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
layout(matrix(c(1,2), 2, 1, byrow = TRUE))

### PC komponenty
hist(ind$PC1, main="", xlab="", ylab="",col="gray", xlim=c(-10,40), breaks=seq(-10,40,by=2.5), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5);hist(ind$PC2, main="", xlab="", ylab="", col="gray",xlim=c(-10,40), breaks=seq(-10,40,by=2.5), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

### Environmentalni promenne
hist(Dataset$AGE_5, main="Stand age", xlab="year", ylab="",col="white", xlim=c(0,70), breaks=seq(0,70,by=2.5), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(Dataset$ALTITUDE_M, main="Elevation", xlab="m a.s.l.", ylab="",col="white", xlim=c(200,1200), breaks=seq(200,1200,by=25), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(Dataset$N, main="N-deposition", xlab=parse(text="g.m^-2*.year^-1"), ylab="",col="white", xlim=c(0.5,2), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(Dataset$S, main="S-deposition", xlab=parse(text="g.m^-2*.year^-1"), ylab="",col="white", xlim=c(0.5,2), breaks=seq(0.5,2,by=0.05), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(Dataset$GEOLOGY, main="Geology", xlab="-", ylab="",col="white", xlim=c(0,1), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist((Dataset$SR/1000), main="Solar radiation", xlab=parse(text="kWh.m^-2"), ylab="",col="white", xlim=c(800,1300), breaks=seq(800,1300,by=12.5), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)

### Promenne zdravotniho stavu lesa
hist(def.nestand$R_DRYTREETOP, main="Dry tree top", xlab="%", ylab="",col="white", xlim=c(0,100), breaks=seq(0,100,by=5), ylim=c(0,200), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(def.nestand$R_REDUCEDINCREMENT, main="Reduced increment", xlab="%", ylab="",col="white", xlim=c(0,100), breaks=seq(0,100,by=5),ylim=c(0,200), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(def.nestand$R_RESIN, main="Resin exudation", xlab="%", ylab="",col="white", xlim=c(0,100), breaks=seq(0,100,by=5),ylim=c(0,200), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(def.nestand$INDEX_COLOR, main="Discoloration", xlab="%", ylab="",col="white", xlim=c(0,100), breaks=seq(0,100,by=5),ylim=c(0,200), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(def.nestand$R_PEELING, main="Peeling", xlab="%", ylab="",col="white", xlim=c(0,100), breaks=seq(0,100,by=5),ylim=c(0,200), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)
hist(def.nestand$R_BREAK, main="Crown breaks", xlab="%", ylab="",col="white", xlim=c(0,100), breaks=seq(0,100,by=5),ylim=c(0,200), cex.lab=1.5, cex.axis=2, cex.main=2.0, cex.sub=2.5)



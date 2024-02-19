
#Variables:
#Fipsno        (State fips*1000)+County fips
#South         Counties in the southern region, scored 1,
#hr            Homicide Rate per 100,000 (numerator is a 3 year
#                              average centered on the decennial census year,
#                              e.g., 1959, 1960, and 1961). Individual deaths
#                              are aggregated to the county level according to
#                              the decendent's county of residence.
#hc            Homicide count (3 year average centered on
#                             decennial census year, e.g., 1959, 1960, and
#                             1961). Individual deaths are aggregated to the
#                             county level according to the decendent's county
#                             of residence.
#po            Population of each county for the decennial census year
#rd            Resource Deprivation/Affluence Component
#                       (prinicipal component composed of percent black,
#                        log of median family income, gini index of
#                        family income inequality, percent of families
#                        female headed (percent of families single parent
#                 for 1960), and percent of families below poverty
#               (percent of families below $3,000 for 1960)) (See Land et al., 1990)
#ps              Population Structure Component(prinicipal
#                                  component composed of the log of population and
#                                  the log of population density) (See Land et al., 1990)
#ue              Percent of civilian labor force that is unemployed
#dv              Percent of males 14 and over who are divorced (aged 15 and over for 1980 and 1990)
#ma              Median age
#pol             Population logged
#dnl             Population density logged
#mfil            Median family income logged
#fp              Percent of families below poverty (percent of families below $3,000 for 1960)
#blk             Percent black
#gi              Gini index of family income inequality
#fh              Percent of families female headed (percent of families single parent for 1960)


library(sf)
library(sp)
library(spdep)
library(maptools)

# Percorso diretto al file shapefile
file_shapefile <- "C:/Users/manet/OneDrive/Desktop/Spatial Data/progetto/natregimes/Missisipi.shp"

# Leggere lo shapefile
data <- st_read(file_shapefile)

data = data[,-c(1, 3, 4, 5, 6, 7, 8, 9, 10)]

plot(data$geometry, col = data[, 1])

# Calcola i centroidi delle geometrie
centroidi <- st_centroid(data$geometry)

# Plotto
plot(data$geometry, main = "Map of Mississippi Counties")
text(st_coordinates(centroidi), labels = data$NAME, cex = 0.5, check.overlap = TRUE)

# creo vicinato (condivisione di un confine in comune)
data.nb <- poly2nb(data)
unlist(data.nb)
summary(data.nb)

plot(data$geometry)
plot(data.nb, centroidi, add=TRUE)

data.nb.lag<-nblag(data.nb,maxlag=5)

plot(data$geometry,col="gray")
plot(data.nb.lag[[1]],centroidi, add=TRUE,col=1)

plot(data$geometry,col="gray")
plot(data.nb.lag[[2]],centroidi, add=TRUE,col=2)

plot(data$geometry,col="gray")
plot(data.nb.lag[[3]],centroidi, add=TRUE,col=3)

plot(data$geometry,col="gray")
plot(data.nb.lag[[4]],centroidi, add=TRUE,col=4)

plot(data$geometry,col="gray")
plot(data.nb.lag[[5]],centroidi, add=TRUE,col=5)

#vicinato con knn
library(RANN)

col.knn <- knearneigh(centroidi, k=2)
plot(data$geometry, border="grey")
data.nb<-knn2nb(col.knn)
plot(data.nb, centroidi, add=TRUE)


col.knn <- knearneigh(centroidi, k=4)
plot(data$geometry, border="grey")
data.nb<-knn2nb(col.knn)
plot(data.nb, centroidi, add=TRUE)

# vicinato per griglie regolari

data.nb = cell2nb(7,7,type="rook",torus=TRUE)

data.nb = cell2nb(7,7,type="queen",torus=FALSE)

# Spatial Weights (dal vicinato alla matrice di adiacenze)
# default standardizzata (tipo W, altrimenti va specificato tipo B)

lw_W <- nb2listw(data.nb) #s0, s1, s2 sono le somme dei pesi(vedi slide)
lw_W

lw_B <- nb2listw(data.nb, style = "B") #qui i pesi non sono standardizzati quindi s0 non sara piu la somma delle righe
lw_B

#qual è il range dei pesi?
1/rev(range(card(lw_W$neighbours))) #da 1/10 (10 numero massimo vicini)= 0.1

summary(unlist(lw_W$weights))
summary(sapply(lw_W$weights, sum))

summary(unlist(lw_B$weights))
summary(sapply(lw_B$weights, sum))

#studiamo la y
brks <- round(quantile(data$HR90, probs = seq(0, 1, by = 0.2), na.rm = TRUE) , digits = 2)
cols <- rev(heat.colors(5))

plot(data$geometry, col=cols[findInterval(data$HR90, brks, all.inside=TRUE)])
title(main="Homicide Rate per 100,000")

min_val <- round(min(data$HR90, na.rm = TRUE), digits = 2)
max_val <- round(max(data$HR90, na.rm = TRUE), digits = 2)

leg_labels <- c(paste(min_val, "-", brks[2]), 
                paste(brks[2], "-", brks[3]), 
                paste(brks[3], "-", brks[4]), 
                paste(brks[4], "-", brks[5]), 
                paste(brks[5], "-", max_val))


legend(x="bottomleft", legend=leg_labels, fill=cols, bty="n")


moran.test(data$HR90, nb2listw(data.nb)) #test con le permutazioni sui dati originali

moran.test(data$HR90, nb2listw(data.nb),randomisation=FALSE) #distribuzione asintotica del moran

geary.test(data$HR90, nb2listw(data.nb))

lA <- lag.listw(nb2listw(data.nb), data$HR90) 
summary(lA)
cor(lA, data$HR90)

moran.plot(data$HR90, nb2listw(data.nb),labels=data$NAME) #in alto a destra hotspot in basso a sinistra coldspot
lmii<-localmoran(data$HR90,nb2listw(data.nb))
lmii

brks <- c(-1,0,0.5,1,1.5,2,3)
cols <- grey((length(brks):2)/length(brks))

plot(data$geometry, col=cols[findInterval(lmii, brks, all.inside=TRUE)])
title(main="Local Moran")
legend(x="bottomleft", leglabs(brks), fill=cols, bty="n")

#togliamo il trend
HR.lm <- lm(HR90 ~ BLK90 + DV90 + FP89, data=data) #modello lineare
summary(HR.lm)
res <- residuals(HR.lm)

brks <- c(min(res),-2,-1,0,1,2,max(res))
cols <- rev(cm.colors(6))

plot(data$geometry, col=cols[findInterval(res, brks, all.inside=TRUE)])
title(main="Regression residuals")
legend(x="bottomleft", legend=leglabs(brks), fill=cols,bty="n")

#provo diverse matrici W
data.nb <- poly2nb(data)
data.nb<- nblag(data.nb, 2)[[2]]
data.nb<-knn2nb(knearneigh(centroidi, k=2))
data.nb<-knn2nb(knearneigh(centroidi, k=4))
data.nb<-knn2nb(knearneigh(centroidi, k=6))

moran.test(res,nb2listw(data.nb)) #fare moran sui residui di un modello di regressione, in genere migliore
lm.morantest(HR.lm, nb2listw(data.nb)) #moran sul modello

moran.plot(res,nb2listw(data.nb), labels=data$NAME) #senza trend ora ho meno elementi significativi


#sar
library(spatialreg)
library(expm)
sar_lm=spautolm(HR90 ~ BLK90 + DV90 + FP89, data=data, listw=nb2listw(data.nb))
summary(sar_lm)

#spatial error model (non significativo prendo spatial lag)
m1=errorsarlm(HR90 ~ BLK90 + DV90 + FP89, data=data, listw=nb2listw(data.nb))
summary(m1)

#spatial lag model 
m2=lagsarlm(HR90 ~ BLK90 + DV90 + FP89, data=data, listw=nb2listw(data.nb))
summary(m2)

#spatial durbin model (non significativo prendo spatial lag)
m3=lagsarlm(HR90 ~ BLK90 + DV90 + FP89, data=data, listw=nb2listw(data.nb), type="mixed")
summary(m3)

AIC(m1);AIC(m2);AIC(m3);AIC(sar_lm) #lag migliore 

#test moltiplicatore di lagrange
mod_lin=lm(HR90 ~ BLK90 + DV90 + FP89, data=data)
lm.LMtests(mod_lin, listw=nb2listw(data.nb), test="all")
#differenza non significativa con spatial error, mentre signif. per spatial lag
#i primi due sono quelli non robusti, basterebbe guardare quello, ma anche i robusti confermano
#la scelta dello spatial lag rispetto a uno spatial error

# Estrazione dei residui dal modello M2
residui <- residuals(m2)

# Calcolo di Moran's I per i residui
moran_test <- moran.test(residui, nb2listw(data.nb), randomisation = FALSE)
print(moran_test)

moran.plot(residui,nb2listw(data.nb), labels=data$NAME)

#inserendo queste covariate e la struttura autoregressiva sulla risposta 
#siamo riusciti a spiegare l'intera struttura spaziale di questo fenomeno



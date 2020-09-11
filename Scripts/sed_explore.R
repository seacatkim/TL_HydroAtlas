## Explore sediment plume data 
## interpolated and no interpoloation
## land is 0 not null
## CKim
## Sept 2020

# read sedientation data without interpoloated
sed <- read.csv('data/sedplume_trans.csv')
sed$Trans <- as.factor(sed$Trans)
sed[which(sed$RASTERVALU ==0),]

# read in interpolated sedimentation data
sed_int <- read.csv('data/sedplume_int_trans.csv')
sed_int$Trans <- as.factor(sed_int$Trans)
sed_int[which(sed_int$RASTERVALU ==0),]


par(mfrow = c(2,1))
plot(RASTERVALU ~ id_group, data = sed)

plot(RASTERVALU ~ id_group, data = sed_int)

library(ggplot2)
ggplot(sed, aes(Trans, RASTERVALU)) +
   geom_boxplot()

ggplot(sed_int, aes(Trans, RASTERVALU)) +
   geom_boxplot()

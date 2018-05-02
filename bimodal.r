#######################
#
# Full Code for Analysis of Bimodality vs Metallicity - M87 / NGC4486
#
#######################
require(mixtools)
require(ggplot2)
require(mclust)
source('/home/emilio/HAWKI/products-inter-combine/functions_to_bimodal.r')

corr_cat = read.table('/home/emilio/HAWKI/products-inter-combine/oldham_chefs.dat')

u = corr_cat$V3
uerr = corr_cat$V4
g = corr_cat$V5
gerr = corr_cat$V6
r = corr_cat$V7
rerr = corr_cat$V8
i_band = corr_cat$V9
i_band_err = corr_cat$V10
z = corr_cat$V11
zerr = corr_cat$V12
k = corr_cat$V13
kerr = corr_cat$V14

#removing outliers in the k band:

bad_index = c()
bad_index = which(kerr >= 1.000)

df_cols = data.frame(u, g, r, i_band, z, k, uerr, gerr, rerr, i_band_err, zerr, kerr)

cols = remove_out(df_cols, bad_index)

u = cols$u
uerr = cols$uerr
g = cols$g
gerr = cols$gerr
r = cols$r
rerr = cols$rerr
i_band = cols$i_band
i_band_err = cols$i_band_err
z = cols$z
zerr = cols$zerr
k = cols$k
kerr = cols$kerr

#bad_index = c()
#for(item in kerr){
#  if(item > 1.0){
#    bad_index = append(bad_index, match(item, kerr))
#  }
#}

#cols = remove_out(df_cols, bad_index)

######### probaibility of being on the MW or being a star, checks

pMW = corr_cat$V21
pStar = corr_cat$V22

mw = c()
star = c()

for(i in 1:length(pMW)){
  if(pMW[i] > 0.0000){
    mw = append(mw, i)
  }
}

#### colors

zk = z-k
gk = g-k
uk = u-k       #with outliers
gz = g-z

zkpure = zk[-mw]
gkpure = gk[-mw]
ukpure = uk[-mw]  #removed outliers
gzpure = gz[-mw]

kpure = k[-mw]
kerrpure = kerr[-mw]

########### error in colors

errzk = sqrt(zerr^2+kerr^2)
errgk = sqrt(gerr^2+kerr^2)
erruk = sqrt(uerr^2+kerr^2)
errgz = sqrt(gerr^2+zerr^2)

errzkpure = errzk[-mw]
errgkpure = errgk[-mw]
errukpure = erruk[-mw]
errgzpure = errgz[-mw]

######### bimodality tests

estgz = normalmixEM(gz)
estzk = normalmixEM(zk)
estgk = normalmixEM(gk)
estuk = normalmixEM(uk)

plot(estzk, which=2)
plot(estgk, which=2)
plot(estuk, which=2)
plot(estgz, which=2)

##### using ggplot to plot bimodality tests

gg.mixEM(estzk, '(z-k)')
gg.mixEM(estgk, '(g-k)')
gg.mixEM(estuk, '(u-k)')
gg.mixEM(estgz, '(g-z)')

estzkpure = normalmixEM(zkpure)
gg.mixEM(estzkpure, '(z-k)')
print(estzkpure$loglik)
estgkpure = normalmixEM(gkpure)
gg.mixEM(estgkpure, '(g-k)')
print(estgkpure$loglik)
estukpure = normalmixEM(ukpure)
gg.mixEM(estukpure, '(u-k)')
print(estukpure$loglik)
estgzpure = normalmixEM(gzpure)
gg.mixEM(estgzpure, '(g-z)')
print(estgzpure$loglik)

dat = data.frame(zk, gk, uk, gz)
datpure = data.frame(zkpure, gkpure, ukpure, gzpure)
ggplot(datpure, aes(x=zkpure)) + geom_histogram(aes(y=..density..), fill='darkgrey')+geom_density()+labs(title='M87 GC (z-k) color distribution', x='(z-k)', y='Density')

##### histograms according to k-band error 

Kerr_hist(gzpure, kerrpure, 0.03, 0.06, '(g-z)')
Kerr_hist(zkpure, kerrpure, 0.03, 0.06, '(z-k)')
Kerr_hist(gkpure, kerrpure, 0.03, 0.06, '(g-k)')
Kerr_hist(ukpure, kerrpure, 0.03, 0.06, '(u-k)')

###############################
#
# CaT to [Fe/H] conversion & Analysis for M87/NGC4486
#
################################

catalog = read.table('/home/emilio/HAWKI/products-inter-combine/corr_romanowski_oldham.dat')

df = data.frame(catalog)

names(df) = c('index', 'RA', 'DEC', 'i', '(g-i)', 'err_col', 'CaT', 'err_cat', 'cross', 'indexOld', 'RA2', 'DEC2','u', 'uerr', 'g','gerr', 'r','rerr', 'i', 'ierr', 'z', 'zerr', 'sep' )

FeH = -3.641+0.438*df$CaT                  #Foster et. al, 2010
FeHerr = sqrt((0.438^2)*(df$err_cat^2))
#FeHerr = df$err_cat

gi_rom_old = df$g-df$i

errgi = sqrt(df$gerr^2+df$ierr^2)

plot(gi_rom_old,df$CaT, pch=5, ylab='CaT (A)', xlab='(g-i)')
arrows(gi_rom_old, df$CaT-df$err_cat, gi_rom_old, df$CaT+df$err_cat, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(gi_rom_old-errgi, df$CaT, gi_rom_old+errgi, df$CaT, length=0.05, angle=90, code=3, lwd= 0.3)

        #removing strange outliers in (g-i)

out_indexes = c(3,25)  #less than 0.0 [25] and more than 1.0 [3]

gi_rom_new = gi_rom_old[-out_indexes]
CaTnew = df$CaT[-out_indexes]
errginew = errgi[-out_indexes]
errcatnew = df$err_cat[-out_indexes]
fehnew = FeH[-out_indexes]
errfehnew = FeHerr[-out_indexes]

plot(gi_rom_new,CaTnew, pch=20, ylab='CaT (A)', xlab='(g-i)', xlim=c(0.6, 1.1), ylim=c(1, 9))
arrows(gi_rom_new, CaTnew-errcatnew, gi_rom_new, CaTnew+errcatnew, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(gi_rom_new-errginew, CaTnew, gi_rom_new+errginew, CaTnew, length=0.05, angle=90, code=3, lwd= 0.3)


plot(gi_rom_new,fehnew, pch=20, ylab='[Fe/H]', xlab='(g-i)', xlim=c(0.65, 1.0),ylim=c(-3.0, 0.4))
arrows(gi_rom_new, fehnew-errfehnew, gi_rom_new,fehnew+errfehnew, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(gi_rom_new-errginew, fehnew, gi_rom_new+errginew, fehnew, length=0.05, angle=90, code=3, lwd= 0.3)

##### plot K vs Kerr
kfit = lm(kerrpure ~ kpure+I(kpure^2)+I(kpure^3))

plot(kpure, kerrpure, xlab='k', ylab='Kerr')
lines(sort(kpure), fitted(kfit)[order(kpure)], col='red', lwd=2) 
legend('topleft', c('Cubic fit'), lty=c(1), col=c('red'))

k_sim_err <- function(mag){
  kse = kfit$coefficients[1]+kfit$coefficients[2]*mag+(kfit$coefficients[3]*mag^2)+(kfit$coefficients[4]*mag^3)
  return(kse)
}

gse = k_sim_err(g)
zse = k_sim_err(z)

gsepure = gse[-mw]
zsepure = zse[-mw]

##### from g-z, obtain simulations of (g-k) and (z-k). See Chies-Santos+2012a, sec. 5.1

gk_sim = (gzpure + 0.349)/0.465
zk_sim = (gzpure + 0.387)/0.746

zk_err_sim = get_err_vecs(zk_sim, zse, 0.3, 0.5)
gk_err_sim = get_err_vecs(gk_sim, gse, 0.5, 0.8)

zk_sim_less_1 = zk_err_sim[[1]]
zk_sim_less_half = zk_err_sim[[2]]
gk_sim_less_1 = gk_err_sim[[1]]
gk_sim_less_half = gk_err_sim[[2]]

Kerr_hist(zk_sim,zse, 0.3, 0.5, '(z-k) (simulated)')
Kerr_hist(gk_sim,gse, 0.5, 0.8, '(g-k) (simulated)')

estzksim = normalmixEM(zk_sim)
gg.mixEM(estzksim, '(z-k) (simulated)')
print(estzksim$loglik)
estzksim.1 = normalmixEM(zk_sim_less_1)
gg.mixEM(estzksim.1, '(z-k) (simulated, for Kerr < 0.5)')
print(estzksim.1$loglik)
estzksim.2 = normalmixEM(zk_sim_less_half)
gg.mixEM(estzksim.2, '(z-k) (simulated, for Kerr < 0.3)')
print(estzksim.2$loglik)

estgksim = normalmixEM(gk_sim)
gg.mixEM(estgksim, '(g-k) (simulated)')
print(estgksim$loglik)
estgksim.1 = normalmixEM(gk_sim_less_1)
gg.mixEM(estgksim.1, '(g-k) (simulated, for Kerr < 0.8)')
print(estgksim.1$loglik)
estgksim.2 = normalmixEM(gk_sim_less_half)
gg.mixEM(estgksim.2, '(g-k) (simulated, for Kerr < 0.5)')
print(estgksim.2$loglik)


########### ploting (g-z) vs (z-k) 

color_color(zk, gz, errzk, errgz, '(z-k)', '(g-z)')
color_color(zk_sim, gzpure, zsepure, errgzpure, '(z-k)(simulated)', '(g-z)(simulated)')   ####not meaningful!

## perform Parametric Bootstrap Cross-Fitting Method (PBCM)

#pbcm(-enter a color distribution-)

pbcm(g-i_band)
pbcm(z-k)

## running GMM (externaly)

setwd('/home/emilio/Downloads/gmm/')

zk_err_vecs = get_err_vecs(zkpure,errzkpure, 0.05, 0.03)
gk_err_vecs = get_err_vecs(gkpure,errgkpure, 0.05, 0.03)
gz_err_vecs = get_err_vecs(gzpure,errgzpure, 0.05, 0.03)

zk_less_1 = zk_err_vecs[[1]]
zk_less_half = zk_err_vecs[[2]]
gk_less_1 = gk_err_vecs[[1]]
gk_less_half = gk_err_vecs[[2]]

colors_errors = list(zkpure, zk_less_1, zk_less_half, gkpure, gk_less_1, gk_less_half)

file_names = c('zk_gmm_pure.dat', 'zk_gmm_less1.dat', 'zk_gmm_lesshalf.dat', 'gk_gmm_pure.dat', 'gk_gmm_less1.dat', 'gk_gmm_lesshalf.dat')

# the following loop takes about 30 min. 

for(i in 1:length(file_names)){
    ptm = proc.time()
    write(colors_errors[[i]], file=file_names[i], ncolumns=1)
    set.seed(7803)
    estaux = normalmixEM(colors_errors[[i]])
    mean1 = estaux$mu[1]
    mean2 = estaux$mu[2]
    if(as.integer(mean1)==as.integer(mean2)){mean2 = mean2+1}
    system2('./gmm', args=c(file_names[i], '0', as.character(as.integer(mean1)), as.character(as.integer(mean2)), paste0('> gmm_output_', file_names[i])))
    system2('mv',  args=c('peakprob.out', paste0('peakprob_',file_names[i])))
    system2('./dip', args=c(as.character(length(colors_errors[[i]])), file_names[i], paste0('> dipout_',file_names[i])))
    print(proc.time()-ptm)/60
}

#### bin GC by radius then plot histograms and bimodality tests on each bin - update: now generates files for GMM for each bin

dummybins = list()

RA = corr_cat$V1
DEC = corr_cat$V2

RA = RA[-c(mw, bad_index)]
DEC = DEC[-c(mw, bad_index)]

hg = radialbins(RA, DEC, zkpure, kerrpure, 1, '(z-k)')
gzbins(hg$RA_less_half, hg$DEC_less_half, 6, hg$est_less_half, set_bins=FALSE, dummybins)

#gzbins(RA, DEC, 5, gzpure, set_bins=FALSE, dummybins)
#gzbins(RA, DEC, 5, zkpure, set_bins=FALSE, dummybins)

### delete the last line of gmm files and...

rungmm('gmm-gzbin5', 'gzbin5')   #system2 not recognizing gmm, don't know why

#### do the same analysis here for LIRIS sample

RAliris = read.table('/home/emilio/HAWKI/products-inter-combine/ana-complete-m87.cat')$V4
DECliris = read.table('/home/emilio/HAWKI/products-inter-combine/ana-complete-m87.cat')$V5
gana = read.table('/home/emilio/HAWKI/products-inter-combine/ana-complete-m87.cat')$V6
zana = read.table('/home/emilio/HAWKI/products-inter-combine/ana-complete-m87.cat')$V8

gzana = gana - zana

gzbins(RAliris, DECliris, 5, gzana, set_bins=FALSE, dummybins)

## do the same for the whole Oldham&Auger catalog

RAold = read.table('/home/emilio/HAWKI/products-inter-combine/oldham_RADECgz.cat', sep='\t')$V1
DECold = read.table('/home/emilio/HAWKI/products-inter-combine/oldham_RADECgz.cat', sep='\t')$V2
gold = read.table('/home/emilio/HAWKI/products-inter-combine/oldham_RADECgz.cat', sep='\t')$V3
zold = read.table('/home/emilio/HAWKI/products-inter-combine/oldham_RADECgz.cat', sep='\t')$V5

gzold = gold - zold

NAindex = which(gzold %in% c(NA))
RAold = RAold[-NAindex]
DECold = DECold[-NAindex]
gzold = gzold[-NAindex]

rbins = read.table('temp_bins_Rpy')$V1

gzbins(RAold, DECold, 5, gzold, set_bins=TRUE, rbins)

###### spatial plot 

plot(RA, DEC)

r = read.table('/home/emilio/HAWKI/products-inter-combine/temp_radius_Rpy')$V1

RAbim = c()
DECbim = c()

for(j in 1:length(r)){
  if((19.34>=r[j])&(r[j]>22.36)){
     RAbim = append(RAbim, RA[j]) #to do: not working?
     DECbim = append(DECbim, DEC[j])
  }
}

points(RAbim, DECbim, color='red', pch=20)

RAgal = 187.7042
DECgal = 12.3911

plot(x=RAgal, y=DECgal, pch=4)
points(RA, DEC)

require(plotrix)
draw.circle(RAgal, DECgal, 2.91, border='black')
draw.circle(RAgal, DECgal, 10.24, border='black')
draw.circle(RAgal, DECgal, 13.76, border='black')
draw.circle(RAgal, DECgal, 16.94, border='black')
draw.circle(RAgal, DECgal, 19.34, border='red')
draw.circle(RAgal, DECgal, 22.36, border='red')
draw.circle(RAgal, DECgal, 30.54, border='black')

legend('topright', c('region with low bimodality'), col=('red'), pch=1)

##################### Kerr (hawki) vs Kerr (liris)

RAmatch = read.table('~/HAWKI/products-inter-combine/match_liris_hawki_new.dat')$V1
DECmatch = read.table('~/HAWKI/products-inter-combine/match_liris_hawki_new.dat')$V2
Kerr_hawki = read.table('~/HAWKI/products-inter-combine/match_liris_hawki_new.dat')$V16 
Kerr_liris = read.table('~/HAWKI/products-inter-combine/match_liris_hawki_new.dat')$V4

Kerr_hawki_op = c()
Kerr_liris_op = c()

for(i in 1:length(Kerr_hawki)){
  if(Kerr_hawki[i] < 0.06){
    Kerr_hawki_op = append(Kerr_hawki[i], Kerr_hawki_op)
    Kerr_liris_op = append(Kerr_liris[i], Kerr_liris_op)
  }
}

plot(Kerr_liris, Kerr_hawki, xlab='Kerr (LIRIS)', ylab='Kerr (HAWKI)', ylim=c(0,0.1))
kerrfit = lm(Kerr_hawki ~ Kerr_liris)
abline(kerrfit, col='red', lty=2)
kerrfitbest = lm(Kerr_hawki_op ~ Kerr_liris_op)
abline(kerrfitbest, col='blue', lty=1)
legend('topleft', c('best fit for full sample', 'best fit excluding Kerr(HAWKI) > 0.06'), col=c('red', 'blue'), lty=c(2, 1))


######## gerr vs g

gfit = lm(gerr ~ g+I(g^2)+I(g^3))

plot(g, gerr, xlab='g', ylab='g_err')
lines(sort(g), fitted(gfit)[order(g)], col='red', lwd=2) 
legend('topleft', c('Cubic fit'), lty=c(1), col=c('red'))

####### CMDs

plot(gk, g, ylim=c(24,19),xlab='(g-k)', ylab='g ', pch=5)
arrows(gk, g-gerr, gk, g+gerr, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(gk-errgk, g, gk+errgk, g, length=0.05, angle=90, code=3, lwd= 0.3)

zpure = z[-mw]
zin = match(zkpure, zk)
zkout = zk[-zin]
zout = z[-zin]
plot(zkpure, zpure, ylim=c(23.3,18), xlab='(z-k)', ylab='z ', pch=20, col='black')
points(zkout, zout, ylim=c(23.3,18), xlab='(z-k)', ylab='z ', pch=20, col='red')
legend('topleft', c('p_MW > 0'), col=c('red'), pch=20)
#pca test script

require(stats)
require(ggplot2)
source('/home/emilio/HAWKI/products-inter-combine/functions_to_bimodal.r')

corr_cat = read.table('/home/emilio/HAWKI/products-inter-combine/chefs_FULL_catalog.dat')

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
k = corr_cat$V23
kerr = corr_cat$V24

#removing outliers on the k band:

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

pMW = corr_cat$V17
pStar = corr_cat$V19

mw = c()
star = c()

for(i in 1:1000){
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

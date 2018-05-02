###############################################
###  Correlate catalogs to find zp!
###  Author = Emilio Zanatta
###############################################

import numpy as np
import matplotlib.pyplot as plt
import rpy2.robjects as robjects

robjects.r('''
        lin_model <- function(x, y){
        
        library(stats)
        
        xx = unlist(x)
        yy = unlist(y) 
        
        lin = lm(yy ~ xx)
        
        print(summary(lin))
        
        png(file = '~/HAWKI/products-inter-combine/fit_plot.png')
        plot(xx, yy, xlab='HAWKI', ylab='LIRIS')
        abline(lin, color='red')
        
        png(file = '~/HAWKI/products-inter-combine/rzp_plot.png')
        plot(x, y, xlab='HAWKI', ylab='LIRIS', main='Hawki vs Liris Ks Mag for M87')
        abline(lin)
        dev.off()

        b1 = coef(lin)[1]             
        b2 = coef(lin)[2]

        ret = list(b1, b2)            
            
        return(ret)
        }
''')

def read_catalog(cat):
    ra = np.loadtxt(cat, usecols=(0,))
    dec = np.loadtxt(cat, usecols=(1,))
    k = np.loadtxt(cat, usecols=(12,))
    kerr = np.loadtxt(cat, usecols=(13,))
    k2 = np.loadtxt(cat, usecols=(29,))
    k2_err = np.loadtxt(cat, usecols=(30,))
    
    return ra, dec, k, kerr, k2, k2_err
    
ra, dec, kh, kh_err, kl, kl_err = read_catalog('low_errors_zp.dat')

hawki = kh.tolist()
liris = kl.tolist()

#indexes = list()
#
#for i in range(0, len(hawki)):
#    if hawki[i] > 0.0:
#       indexes.append(i)
#      
#hawki = [i for j, i in enumerate(hawki) if j not in indexes]
#liris = [i for j, i in enumerate(liris) if j not in indexes]

lin_model = robjects.r('lin_model')

coef = lin_model(hawki, liris)

b1 = coef[0]
b2 = coef[1]

print 'zero point estimated:', b1

raf = np.loadtxt('chefs_oldham_m87.dat', usecols=(0,))
decf = np.loadtxt('chefs_oldham_m87.dat', usecols=(1,))
kf = np.loadtxt('chefs_oldham_m87.dat', usecols=(12,))

k_adj = b2*kf+b1

with open('new_catalog.dat', 'w') as f:
    with open('correlation.dat', 'w') as g:

        for i in range(0,len(ra)):
            print >>f, ra[i], dec[i], k_adj[i]
        for i in range(0, len(raf)):
            print >>g, raf[i], decf[i], k_adj[i]   

        
    
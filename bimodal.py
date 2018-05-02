import numpy as np
import rpy2.robjects as r

r.r('''

read_catalog <- function(){

  cat1 = 'oldham_and_auger_m87.cat'

  catalog = read.table(cat1, sep="\t")
  
  ra = catalog$V2
  dec = catalog$V3
  u = catalog$V4
  u_err = catalog$V5
  g = catalog$V6
  g_err = catalog$V7
  r = catalog$V8
  r_err = catalog$V9
  i = catalog$V10
  i_err = catalog$V11
  z = catalog$V12    
  z_err = catalog$V13
  
  cat2 = 'new_catalog.dat'  

  catalog2 = read.table(cat2)

  ra_h = catalog2$V1
  dec_h = catalog2$V2
  k = catalog2$V3
  k_err = catalog2$V4

  cat1 = data.frame(RAold = ra, DECold = dec, u, u_err, g, g_err, r, r_err,z, z_err, i, i_err)      
        
  return(cat1)
}''')

r.r('''

read_catalog_hawki <- function(){
  
  cat2 = 'new_catalog.dat'  

  catalog2 = read.table(cat2)

  ra_h = catalog2$V1
  dec_h = catalog2$V2
  k = catalog2$V3
  k_err = catalog2$V4

  cat2 = data.frame(RAnew=ra_h, DECnew=dec_h, k, k_err)      
        
  return(cat2)
}''')

r.r('''

plot_bimodality <- function(cat1, cat2){
    
    


}


''')

read_cat = r.r('read_catalog')
read_catalog_hawki = r.r('read_catalog_hawki')
    
oldham = read_cat()
hawki = read_catalog_hawki()

robj = list()
dobj = list()
u = list()
g = list()
r = list()
i = list()
z = list()
k = list()

thres = .0001

for i in range(0, ra.size):
    for j in range(0, ra_c.size):
        if (abs(ra[i] - ra_c[j]) < thres and abs(dec[i]-dec_c[j]) < thres):       
            robj.append(ra_c[j])
            dobj.append(dec_c[j])
            u.append(u[j])
            g.append(g[j])
            r.append(r[j])
            i.append(i[j])
            z.append(z[j])
            k.append(k[j])


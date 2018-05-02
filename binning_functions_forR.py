def binning(RA, DEC, RAgal, DECgal, d, nbins):
  import numpy as np
  from astropy import units as u
  from astropy.coordinates import SkyCoord
  
  RA = RA*u.deg
  DEC = DEC*u.deg
  n = len(RA)
  
  gcs = list()
  dist = list()
  d = d*u.pc
  g = SkyCoord(RAgal, DECgal, distance=d)
  
  x1l, y1l, z1l = [],[],[]
  
  for i in range(0, n):
    x1 = d.value*np.cos(RA[i])*np.cos(DEC[i])
    y1 = d.value*np.sin(RA[i])*np.cos(DEC[i])
    z1 = d.value*np.sin(DEC[i])
    
    RAg = g.ra.degree*u.degree
    DECg = g.dec.degree*u.degree
    
    xg = d.value*np.cos(RAg)*np.cos(DECg)
    yg = d.value*np.sin(RAg)*np.cos(DECg)
    zg = d.value*np.sin(DECg)
    
    dx = (x1-xg)
    dy = (y1-yg)
    dz = (z1-zg)
    
    #print x1, y1, z1
    #print xg, yg, zg
    
    aux = np.sqrt(dx*dx+dy*dy+dz*dz)*1000
    
    x1l.append(x1)
    y1l.append(y1)
    z1l.append(z1)
    dist.append(aux)
  
  ###binning ####
  r = sorted(dist)
  n = len(r)   #lenght of r --- number of objects to iterate on the following loops
  
  rgal = list()      #rgal contains the limits of each bin
  for h in range(1, nbins):
    rgal.append(r[int(n*h/nbins)])     #binning
  rgal.insert(nbins, r[n-1])
  rgal.insert(0, r[0])
  
  with open('temp_radius_Rpy', 'w') as f:
    for item in r:
      print >>f, item
  with open('temp_bins_Rpy', 'w') as f:
    for item in rgal:
      print >>f, item
      
  with open('temp_radius_Rpy', 'w') as f:
    for item in x1l:
      print >>f, item
  with open('temp_radius_Rpy', 'w') as f:
    for item in y1l:
      print >>f, item
  with open('temp_radius_Rpy', 'w') as f:
    for item in z1l:
      print >>f, item  
      
def density(nbin, r, rgal):
    import numpy as np
    
    n = len(r)
    NGC = np.zeros(nbin)    #number of GC in each bin
    rmed = [[] for i in range(0,nbin)]
    for h in range(0, nbin):                              
         for i in range(0, n):
             if (rgal[h+1] >= r[i] and r[i] > rgal[h]):
                 NGC[h] = NGC[h] + 1
                 rmed[h].append(r[i])

    area = []
    for h in range(0, nbin):
        area.append((np.pi*rgal[h+1]**2)-(np.pi*rgal[h]**2)) #area per bin
    
    binsize = np.linspace(0,0,nbin)    
    for h in range(0,nbin):
        binsize[h] = rgal[h+1]-rgal[h]  
    binsize = binsize/2
    
    return area, NGC, median, binsize, poi_err 
    

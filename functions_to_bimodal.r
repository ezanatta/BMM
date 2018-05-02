remove_out <- function(df, bad_index){

  u = df$u[-c(bad_index)]
  g = df$g[-c(bad_index)]
  r = df$r[-c(bad_index)]
  i_band = df$i_band[-c(bad_index)]
  z = df$z[-c(bad_index)]
  k = df$k[-c(bad_index)]
  uerr = df$uerr[-c(bad_index)]
  gerr = df$gerr[-c(bad_index)]
  rerr = df$rerr[-c(bad_index)]
  i_band_err = df$i_band_err[-c(bad_index)]
  zerr = df$zerr[-c(bad_index)]
  kerr = df$kerr[-c(bad_index)]
  
  df = data.frame(u, g, r, i_band, z, k, uerr, gerr, rerr, i_band_err, zerr, kerr)
  
  return(df)
}


get_err_vecs <- function(est, err, err_thres1, err_thres2){
  est_less_half = c()
  est_less_1 = c()
  
  for(i in 1:length(est)){
    if(err[i] < err_thres1){
      est_less_half = append(est_less_half, est[i])
    }
    if(err[i] < err_thres2){
      est_less_1 = append(est_less_1, est[i])
    }
  }
  
  outp = list(est_less_1, est_less_half)
  
  return(outp)
}

Kerr_hist <- function(est, err, err_thres1, err_thres2, title){
  
  #this function only works if Kerr and K_band magnitudes are already loaded in the main code
  
  est_less_half = c()
  est_less_1 = c()
  
  for(i in 1:length(est)){
    if(err[i] < err_thres1){
      est_less_half = append(est_less_half, est[i])
    }
    if(err[i] < err_thres2){
      est_less_1 = append(est_less_1, est[i])
    }
  }
  
  hist(est, breaks='Sturges', main=NA, xlab=title)
  hist(est_less_1, add=TRUE, color='grey', density=10, breaks='Sturges')
  hist(est_less_half, add=TRUE, color='black', density=100, breaks='Sturges')
  legend('topleft', c('full sample', paste0('K_err < ',err_thres1), paste0('K_err < ', err_thres2)), fill=c('white', 'black', 'black'), density=c(NA, 10, 100))

  #ests = data.frame(est_less_half, est_less_1)
  return(est_less_half)
}

radialbins <- function(RA, DEC, est, err, err_thres, title){
  
  #this function only works if Kerr and K_band magnitudes are already loaded in the main code
  
  est_less_half = c()
  est_less_1 = c()
  RA_less_half = c()
  DEC_less_half = c()
  
  for(i in 1:length(est)){
    if(err[i] < err_thres){
      est_less_half = append(est_less_half, est[i])
      RA_less_half = append(RA_less_half, RA[i])
      DEC_less_half = append(DEC_less_half, DEC[i])
    }
  }
  
  hist(est, breaks='Sturges', main=NA, xlab=title)
  hist(est_less_half, add=TRUE, color='black', density=100, breaks='Sturges')
  legend('topleft', c('full sample', paste0('K_err < ',err_thres)), fill=c('white', 'black'), density=c(NA, 10))
  
  halfgc = data.frame(RA_less_half, DEC_less_half, est_less_half)
  return(halfgc)
}


gg.mixEM <- function(EM, color) {
  require(ggplot2)
  x       <- with(EM,seq(min(x),max(x),len=1000))
  pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))
  em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
  em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
  fname = paste0(color, '-bimodaltest_m87.png')
  
  ggplot(data.frame(x=EM$x),aes(x,y=..density..)) + 
    geom_histogram(fill='darkgrey')+
    #geom_density()+
    geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5)+
    scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))+
    theme_bw()+labs(title=' ', x=color, y='Density')
  
  ggsave(fname)
}

pbcm <- function(est1){
  require(mclust)
  
  est1.gmm = Mclust(est1, G=2, modelNames = 'V')
  est1.gmm.1 = Mclust(est1, G=1, modelNames = 'V')
  
  g1 = 0
  g2 = 0
  
  for(item in est1.gmm$class){
    if(item == 1){
      g1 = g1+1
    }else{
      g2 = g2+1
    }
  }
  
  m1.1 = est1.gmm$parameters$mean[1]
  m1.2 = est1.gmm$parameters$mean[2]
  v1.1 = est1.gmm$parameters$variance$sigmasq[1]
  v1.2 = est1.gmm$parameters$variance$sigmasq[2]
  
  m2.1 = est1.gmm.1$parameters$mean[1]
  v2.1 = est1.gmm.1$parameters$variance$sigmasq
  
  print(summary(est1.gmm))
  print(summary(est1.gmm.1))
  
  
  set.seed(7809)
  B = 1000;    x2.d = vector(length=B);    x1.d = vector(length=B)
  for(i in 1:B){
    x2      = c(rnorm(g1, mean=m1.1, sd=sqrt(v1.1)), rnorm(g2, mean=m1.2, sd=sqrt(v1.2)))
    x1      = rnorm(length(est1), mean=m2.1, sd=sqrt(v2.1))
    x2.d[i] = Mclust(x2, G=2)$loglik - Mclust(x2, G=1)$loglik
    x1.d[i] = Mclust(x1, G=2)$loglik - Mclust(x1, G=1)$loglik
  }
  
  x2.d = sort(x2.d);  x1.d = sort(x1.d)

  print(summary(x1.d))
  print(summary(x2.d))

  plot(density(x1.d), xlab='difference in log-likelihood', xlim=c(min(x1.d), max(x2.d)))
  lines(density(x2.d), col='green')
  legend('topright', c('1 component', '2 component'), lty=c(1, 1), col=c('black', 'green'))
}

color_color <- function(col1, col2, errcol1, errcol2, xlabel, ylabel){
  
  plot(col1, col2, xlab=xlabel, ylab=ylabel, pch=5)
  arrows(col1, col2-errcol2, col1, col2+errcol2, length=0.05, angle=90, code=3, lwd= 0.3)
  arrows(col1-errcol1, col2, col1+errcol1, col2, length=0.05, angle=90, code=3, lwd= 0.3)
  #### binning and median plot
  
  number_of_obj_per_bin = 50
  
  col1_or = sort(col1, decreasing=FALSE)
  col2_or = sort(col2, decreasing=FALSE)
  
  col1_bins = split(col1_or, ceiling(seq_along(col1_or)/number_of_obj_per_bin))
  col2_bins = split(col2_or, ceiling(seq_along(col2_or)/number_of_obj_per_bin))
  
  col1_bins = data.frame(col1_bins)
  col2_bins = data.frame(col2_bins)
  
  col1avg = list()
  col2avg = list()
  
  for(i in 1:(1000/number_of_obj_per_bin)){
    col1avg = append(col1avg, median(col1_bins[[i]]))
    col2avg = append(col2avg, median(col2_bins[[i]]))
  }
  
  lines(col1avg, col2avg, type='l', col='green', lwd=2)
  leg.txt = c('M87 GCs', 'Outliers','Median')
  legend('bottomleft',leg.txt, pch=c(5, 5, NA), lty=c(NA, NA, 1), col=c('black', 'red', 'green') )
}

gzbins <- function(RA, DEC, nbins, gz, set_bins=FALSE, bins){
  
  require(rPython)
  python.load('/home/emilio/HAWKI/products-inter-combine/binning_functions_forR.py')
  pa = 153.0
  incl = 31.0
  RAgal = '12h30m49.4s'
  DECgal = '+12d23m28s'
  d = 22.2
  python.call("binning", RA, DEC, RAgal, DECgal, d, nbins)
  
  r = read.table('temp_radius_Rpy')$V1
  
  if(set_bins==TRUE){rbins = bins}else{rbins = read.table('temp_bins_Rpy')$V1}
  
  ibins = replicate(nbins, c())

  for(i in 1:nbins){
    for(j in 1:length(r)){
      if((rbins[i+1]>=r[j])&(r[j]>rbins[i])){
        ibins[[i]] = append(ibins[[i]], j)
      }
    }
  }
   ########## verifying if (g-z) is bimodal in each bin
 
  gzbins = replicate(nbins, c())
  estbins = replicate(nbins, list())
  
  require(propagate)
  for(i in 1:nbins){
    #pbcm(gzbins[[i]])
    
    set.seed(5678)
    gzbins[[i]] = gz[ibins[[i]]]
    est1 = normalmixEM(gzbins[[i]])
    estbins[[i]] = append(estbins[[i]], est1)
    titles = paste0(round(rbins[i], 2), ' < r (kpc) < ', round(rbins[i+1], 2))
    #plot(est1, which=2, xlab2='(z-k)', main2=title, col2=c('blue', 'red'), breaks=20)
    xrand=c(rnorm(10000,8,2),rnorm(10000,17,4))
    sdnorm =
      function(xrand, mean=0, sd=1, lambda=1, n=length(gzbins[[i]]), bw=0.2){lambda*dnorm(xrand, mean=mean, sd=sd)*n*bw}
    
    a = ggplot(data.frame(x=gzbins[[i]])) + 
      geom_histogram(aes(x=gzbins[[i]],y=..count..),fill="white",color="black", binwidth = 0.2) +
      stat_function(fun=sdnorm,
                    args=list(mean=est1$mu[2],
                             sd=est1$sigma[2],
                             lambda=est1$lambda[2]),
                    colour="red",geom="line", size=2) +
      stat_function(fun=sdnorm,
                    args=list(mean=est1$mu[1],
                             sd=est1$sigma[1],
                             lambda=est1$lambda[1]),
                    colour="blue",geom="line", size=2)+
      labs(title=titles, x='(z-k)', y='NGC')+theme(plot.title = element_text(hjust = 0.5))+
      xlim(2, 6)
    
    plot(a)
    
    print(i)
    #to aid gmm usage:
    gmmfile = paste0('/home/emilio/HAWKI/products-inter-combine/gmm-files-bins/gmm-gzbin', i)
    write.table(as.double(format(gzbins[[i]], digits=4)), file=gmmfile, row.names = FALSE, col.names = FALSE)
    print(paste0('kurtosis=', kurtosis(gzbins[[i]])))
    print(paste0('Number of GC = ', length(gzbins[[i]])))
  }
}

rungmm <- function(gmmfile, title){
  gmmfile = paste0('/home/emilio/HAWKI/products-inter-combine/gmm-files-bins/hawki/', gmmfile)
  f = read.table(gmmfile)$V1
  est = normalmixEM(f)
  m1 = est$mu[1]
  m2 = est$mu[2]
  print(m1)
  print(m2)
  system2('./gmm', args=c(gmmfile, '0', as.character(as.integer(m1)), as.character(as.integer(m2)), paste0('> /home/emilio/Downloads/gmm/gmm_output_', title)))
  system2('mv',  args=c('/home/emilio/Downloads/gmm/peakprob.out', paste0('/home/emilio/Downloads/gmm/peakprob_',title)))
  system2('./dip', args=c(as.character(length(f)), gmmfile, paste0('> /home/emilio/Downloads/gmm/dipout_',title)))
}

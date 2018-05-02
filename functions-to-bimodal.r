gg.mixEM <- function(EM, color) {
  require(ggplot2)
  x       <- with(EM,seq(min(x),max(x),len=1000))
  pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))
  em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
  em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
  
  ggplot(data.frame(x=EM$x),aes(x,y=..density..)) + 
    geom_histogram(fill='darkgrey')+
    #geom_density()+
    geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5)+
    scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))+
    theme_bw()+labs(title=color, x=color, y='Density')
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
  
  print(c('n1', g1, 'n2', g2, 'm1.1', m1.1, 'm1.2', m1.2, 'v1.1', v1.1, 'v1.2', v1.2, 'm2.1', m2.1, 'v2.1', v2.1))
  
  set.seed(7809)
  B = 100;    x2.d = vector(length=B);    x1.d = vector(length=B)
  for(i in 1:B){
    x2      = c(rnorm(g1, mean=m1.1, sd=sqrt(v1.1)), rnorm(g2, mean=m1.2, sd=sqrt(v1.2)))
    x1      = rnorm(length(est1), mean=m2.1, sd=sqrt(v2.1))
    x2.d[i] = Mclust(x2, G=2)$loglik - Mclust(x2, G=1)$loglik
    x1.d[i] = Mclust(x1, G=2)$loglik - Mclust(x1, G=1)$loglik
  }
  
  x2.d = sort(x2.d);  x1.d = sort(x1.d)
  
  plot(density(x1.d), xlab='difference in log-likelihood', xlim=c(0, 300))
  lines(density(x2.d), col='green')
  legend('topright', c('1 component', '2 component'), lty=c(1, 1), col=c('black', 'green'))
}
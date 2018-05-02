###############################3
#
# CaT to [Fe/H] conversion & Analysis for M87/NGC4486
#
################################

catalog = read.table('romanowskyCaT.dat')

df = data.frame(catalog)

names(df) = c('index', 'RA', 'DEC', 'i', '(g-i)', 'err_col', 'CaT', 'err_cat', 'cross')

FeH = -3.641+0.438*df$CaT                  #Foster et. al, 2010
FeHerr = -3.641+0.438*df$err_cat

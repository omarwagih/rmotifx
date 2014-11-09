setwd('~/Development/motifx/package/')

remove.packages('motifx')
version = '1.0'

targz = sprintf('motifx/build/motifx_%s.tar.gz', version)

system(sprintf('rm -rf %s', targz))
system('R CMD build motifx')

system(sprintf('mv motifx_%s.tar.gz %s', version, targz))

install.packages(targz, repos = NULL, type="source")
#detach("package:motifx", unload=TRUE)

setwd('~/Development/motifx/')

remove.packages('rmotifx')
version = '1.0'

targz = sprintf('rmotifx/build/rmotifx_%s.tar.gz', version)

system(sprintf('rm -rf %s', targz))
system('R CMD build rmotifx')

system(sprintf('mv rmotifx_%s.tar.gz %s', version, targz))

install.packages(targz, repos = NULL, type="source")
#detach("package:motifx", unload=TRUE)

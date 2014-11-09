# motif-x implementation in R
# Created on : 09-11-14 
# Author     : omarwagih
# Todo       : Implement degenerate enrichment, add DNA support, add non-centered site support

# All amino acids
.GetAA <- function(){
  AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  AA
}

# Constructing amino acids
.BuildPWM <- function(seqs, relative.freq=T){
  
  # Ensure same length characters 
  num.pos = seq.len = nchar(seqs[1])
  if(F){
    # Slows things down for large number of sequences
    seq.len = sapply(seqs, nchar)
    num.pos = seq.len[1]
    if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  }
  
  # List of valid amino acids, sorted
  namespace = .GetAA()
  
  # Make matrix of letters
  pwm.matrix = sapply(1:num.pos, function(col){
    i = col
    # Get frequencies 
    t = table(substr(seqs, col,col))
    # Match to aa
    ind = match(namespace, names(t))
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    names(col) = namespace
    
    
    # Do relative frequencies
    if(relative.freq)
      col = col / sum(col)
    
    col
  })
  colnames(pwm.matrix) = 1:num.pos
  
  return(pwm.matrix)
}

.MotifRegex <- function(motifs, cent.regex='[ST]', n=15){
  if(class(motifs) == 'data.frame') motifs = list(motifs)
  
  # Central residue position
  cent.ind = ceiling(n/2)
  # Must have at least one row
  motifs = motifs[sapply(motifs, nrow) > 0]
  regex = sapply(motifs, function(mt){
    aa = mt[[1]]
    col = mt[[2]]
    z = rep('.', n)
    z[col] = aa
    z[cent.ind] = cent.regex
    paste0(z, collapse='')
  })
  regex
}

.FindMotif <- function(fg.seqs, bg.seqs, min.seqs=20, pval.cutoff=0.05, cent.regex='[ST]', verbose=F, perl.impl=F){
  
  AA = .GetAA()
  
  .fg.seqs = fg.seqs
  .bg.seqs = bg.seqs
  
  fg.size = length(fg.seqs)
  bg.size = length(bg.seqs)
  
  kmer.length = nchar(fg.seqs[1])
  cent.index = ceiling(kmer.length/2)
  
  if(verbose) writeLines('Step 1: recursive motif building')
  motif = matrix(0, 0, 2)
  motif.parts = matrix(0, 0, 2)
  pvals = c()
  iter = 1
  while(TRUE){
    # Break if we dont have any sequences remaining
    if(length(bg.seqs) == 0 | length(fg.seqs) == 0) break
    
    # Construct PFM and PWM
    ppwm = .BuildPWM(fg.seqs, relative.freq=F)
    npwm = .BuildPWM(bg.seqs)
    
    # Compute binomial matrix
    bin = 1 - pbinom(ppwm-1, length(fg.seqs), npwm, lower.tail=TRUE)
    
    rownames(bin) = rownames(ppwm)
    colnames(bin) = colnames(ppwm)
    
    # Lowest p-value computed by pbinom is 10e-16, set any zero values to that
    if(perl.impl){
      if(verbose) writeLines('Using perl implementation.')
      bin[bin == 0] = 1e-16 # this is what weblogo uses in perl, different for R
    }else{
      bin[bin == 0] = .Machine$double.xmin # Usually ~= 2.2e-308
    }
    
    
    # Don't enrich on central residue or previous motifs
    bin[,cent.index] = 1
    bin[motif] = 1
    bin[motif.parts] = 1
    
    # bin[, c(1:2, 14:15)] = 1 # exclude outermost
    
    # Give anything with a frequency less than the threshold a p-value of 1
    bin[ppwm < min.seqs] = 1
    
    # Find the row/column of the min binomal p-value
    min.bin = min(bin)
    mbin = which(bin == min.bin & bin < pval.cutoff, arr.ind=T)
    
    # No more significant, breaking
    if(nrow(mbin) == 0) break
    
    # Case where there are more than one matches (likely pvalue=0)
    # Find match with greatest value in PFM
    if(nrow(mbin) > 1){
      r = which.max(ppwm[mbin])
      rn = rownames(mbin)
      mbin = mbin[r,]
      mbin = t(as.matrix(mbin))
      rownames(mbin) = rn[r]
    }
    
    aa = rownames(mbin)[1]
    row = mbin[1,1]
    col = mbin[1,2]
    
    aa. = aa
    
    # Extract sequences 
    ind.pos = substr(fg.seqs, col, col) %in% aa.
    ind.neg = substr(bg.seqs, col, col) %in% aa.
    
    
    if(verbose)
      writeLines(sprintf('\tIteration %s, aa=%s row=%s col=%s, pval=%s, fg=%s, bg=%s', 
                         iter, aa, row, col, signif(min.bin, 3), length(fg.seqs), length(bg.seqs)))
    
    fg.seqs = fg.seqs[ind.pos]
    bg.seqs = bg.seqs[ind.neg]
    motif = rbind(motif, c(row, col))
    
    
    pvals = append(pvals, min.bin)
    iter = iter + 1
  }
  
  # Motif data: data frame with amino acid and positions for the motif
  motif.data = as.data.frame(motif) 
  names(motif.data) = c('aa', 'col')
  
  AA. = AA
  motif.data$aa = AA.[motif.data$aa]
  
  # If we don't have any data, return NULL so that the parent function can deal with it
  if(nrow(motif.data) == 0) return(NULL)
  
  # Find the regex of the motif given the row/col and aa values
  motif.regex = .MotifRegex(motif.data, cent.regex, kmer.length)
  
  # Compute the score of the motif
  motif.score = ifelse( is.null(pvals), NA, sum( -log10(pvals) ) )
  
  # Compute fg/bg matches
  fg.matches = length(grep(motif.regex, .fg.seqs))
  bg.matches = length(grep(motif.regex, .bg.seqs))
  
  return(list(pos=fg.seqs, motif.data=motif.data, motif.score=motif.score, motif.regex=motif.regex,
              fg.matches=fg.matches, fg.size=fg.size, 
              bg.matches=bg.matches, bg.size=bg.size))
}



#' Find overrepresented sequence motifs
#' 
#' @export 
#' 
#' @param fg.seqs Foreground k-mer sequences in a pre-aligned format. All k-mers must have same lengths.
#' @param bg.seqs Background k-mer sequences in a pre-aligned format. All k-mers must have same lengths.
#' @param central.res Central amino acid of the k-mer. Sequences without this amino acid in the centre position are filtered out. This can be one or more letter. For example, 'S', 'ST', 'Y', or 'STY'.
#' @param min.seqs This threshold refers to the minimum number of times you wish each of your extracted motifs to occur in the data set. An occurrence threshold of 20 usually is appropriate, although this parameter may be adjusted to yield more specific or less specific motifs.
#' @param pval.cutoff The p-value threshold for the binomial probability. This is used for the selection of significant residue/position pairs in the motif. A threshold of 0.000001 is suggested to maintain a low false positive rate in standard protein motif analyses.
#' @param verbose If true, motifx will show textual details of the steps while running.
#' @param perl.impl The original implementation of motifx in perl, P-values below 1e-16 cannot be computed and are thus set to zero. Motifx therefore sets any P-values of zero to the minimal P-value of 1e-16. In R, the minimal P-value is much lower (depending on the machine). If this option is set to TRUE, P-values with a value of zero are set to 1e-16, as in perl. Otherwise, the R P-value minimum will be used. For results identical to that of the webserver implementation, set to TRUE.
#' @return Data frame with seven columns containing overrepresented motifs. Motifs are listed in the order in which they are extracted by the algorithm, not with regard to statistical significance. Thus it should not be assumed that a motif found at a higher position in the list is more statistically significant than a motif found at a lower position in the list. The columns are as follows:
#' \describe{
#'   \item{motif}{The overrepresented motif}
#'   \item{score}{The motif score, which is calculated by taking the sum of the negative log probabilities used to fix each position of the motif. Higher motif scores typically correspond to motifs that are more statistically significant as well as more specific }
#'   \item{fg.matches}{Frequency of sequences matching this motif in the foreground set}
#'   \item{fg.size}{Total number of foreground sequences}
#'   \item{bg.matches}{Frequency of sequences matching this motif in the background set}
#'   \item{bg.size}{Total number of background sequences}
#'   \item{fold.increase}{An indicator of the enrichment level of the extracted motifs. Specifically, it is calculated as (foreground matches/foreground size)/(background matches/background size).}
#' }
#' @examples
#' # Get paths to sample files
#' fg.path = system.file("extdata", "fg-data-ck2.txt", package="motifx")
#' bg.path = system.file("extdata", "bg-data-serine.txt", package="motifx")
#' 
#' # Read in sequences
#' fg.seqs = readLines(fg.path)
#' bg.seqs = readLines(bg.path)
#' 
#' # Find overrepresented motifs
#' mot = motifx(fg.seqs, bg.seqs, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
#' 
#' # View results
#' head(mot)
motifx <- function(fg.seqs, bg.seqs, central.res='ST', min.seqs=20, pval.cutoff=1e-6, verbose=F, perl.impl=F){
  
  AA = .GetAA()
  # Degenerate option not implemented yet.
  .CheckEmptySeqs <- function(){
    if(length(fg.seqs) == 0) stop('Could not find any foreground sequences!')
    if(length(bg.seqs) == 0) stop('Could not find any background sequences!')
  }
  
  .CheckEmptySeqs()
  
  cent.regex = ifelse(nchar(central.res) > 1, paste0('[',central.res,']'), central.res)
  c.res = strsplit(central.res, '')[[1]]
  c.res = intersect(AA, c.res)
  
  if(length(c.res) == 0) stop('Central residue must contain at least one amino acid character')
  
  # Check sequence widths 
  width = nchar(fg.seqs[1])
  if(width < 3 | width > 35) stop('Sequence width must be between 3 and 35!')
  if(width %% 2 == 0) stop('Sequence width must be an odd number!')
  if(width != nchar(bg.seqs[1])) 
    stop('Widths for foreground and background data must be equal!')
  
  # Ensure k-mers have same lengths
  nc.pos = nchar(fg.seqs)
  if( any(nc.pos[1] != nc.pos[-1]) ) stop('Foreground k-mers must be same lentth.')
  nc.bg = nchar(bg.seqs)
  if( any(nc.bg[1] != nc.bg[-1]) ) stop('Background k-mers must be same lentth.')
  
  # Get central index
  ci = ceiling(width/2)
  
  # Filter the central residue for allowed residues
  fg.seqs = fg.seqs[substr(fg.seqs, ci,ci) %in% c.res]
  bg.seqs = bg.seqs[substr(bg.seqs, ci,ci) %in% c.res]
  
  
  if(T){
    # Only set to true for exact match to motifx webserver 
    # Remove non-amino acid residues
    fg.seqs = fg.seqs[!grepl('\\-|\\*|[BJOUXZ]', fg.seqs)]
    bg.seqs = bg.seqs[!grepl('\\-|\\*|[BJOUXZ]', bg.seqs)]
    bg.seqs = unique(bg.seqs)
  }
  
  .CheckEmptySeqs()
  
  seqs = fg.seqs
  data = list()
  while(TRUE){
    # Find the motif
    mt = .FindMotif(fg.seqs = seqs, bg.seqs = bg.seqs, min.seqs = min.seqs, 
                   pval.cutoff = pval.cutoff, cent.regex = cent.regex, verbose = verbose, 
                   perl.impl=perl.impl)
    if(is.null(mt)) break
    
    # Append to list of data
    data = append(data, list(mt))
    # Remove stuff already apart of a motif
    if(verbose) writeLines('Step 2: positive and negative set reduction')
    seqs = seqs[! grepl(mt$motif.regex, seqs)]
    bg.seqs = bg.seqs[! grepl(mt$motif.regex, bg.seqs)]
    
    # No more sequences left to process, breaking.
    if(length(seqs) < min.seqs) break 
  }
  
  
  if(verbose) writeLines('Converged, no more enrichments!')
  
  df = data.frame(motif = sapply(data, function(x) x$motif.regex), 
                  score = sapply(data, function(x) x$motif.score), 
                  fg.matches = sapply(data, function(x) x$fg.matches), 
                  fg.size = sapply(data, function(x) x$fg.size), 
                  bg.matches = sapply(data, function(x) x$bg.matches), 
                  bg.size = sapply(data, function(x) x$bg.size),
                  stringsAsFactors=F)
  
  df$fold.increase = (df$fg.matches/df$fg.size)/(df$bg.matches/df$bg.size)
  rownames(df) = NULL
  if(nrow(df) == 0) return(NULL)
  
  return(df)
}

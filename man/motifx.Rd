\name{motifx}
\alias{motifx}
\title{Find overrepresented sequence motifs}
\usage{
motifx(fg.seqs, bg.seqs, central.res = "ST", min.seqs = 20,
  pval.cutoff = 1e-06, verbose = F, perl.impl = F)
}
\arguments{
  \item{fg.seqs}{Foreground k-mer sequences in a
  pre-aligned format. All k-mers must have same lengths.}

  \item{bg.seqs}{Background k-mer sequences in a
  pre-aligned format. All k-mers must have same lengths.}

  \item{central.res}{Central amino acid of the k-mer.
  Sequences without this amino acid in the centre position
  are filtered out. This can be one or more letter. For
  example, 'S', 'ST', 'Y', or 'STY'.}

  \item{min.seqs}{This threshold refers to the minimum
  number of times you wish each of your extracted motifs to
  occur in the data set. An occurrence threshold of 20
  usually is appropriate, although this parameter may be
  adjusted to yield more specific or less specific motifs.}

  \item{pval.cutoff}{The p-value threshold for the binomial
  probability. This is used for the selection of
  significant residue/position pairs in the motif. A
  threshold of 0.000001 is suggested to maintain a low
  false positive rate in standard protein motif analyses.}

  \item{verbose}{If true, motifx will show textual details
  of the steps while running.}

  \item{perl.impl}{The original implementation of motifx in
  perl, P-values below 1e-16 cannot be computed and are
  thus set to zero. Motifx therefore sets any P-values of
  zero to the minimal P-value of 1e-16. In R, the minimal
  P-value is much lower (depending on the machine). If this
  option is set to TRUE, P-values with a value of zero are
  set to 1e-16, as in perl. Otherwise, the R P-value
  minimum will be used. For results identical to that of
  the webserver implementation, set to TRUE.}
}
\value{
Data frame with seven columns containing overrepresented
motifs. Motifs are listed in the order in which they are
extracted by the algorithm, not with regard to statistical
significance. Thus it should not be assumed that a motif
found at a higher position in the list is more
statistically significant than a motif found at a lower
position in the list. The columns are as follows:
\describe{ \item{motif}{The overrepresented motif}
\item{score}{The motif score, which is calculated by taking
the sum of the negative log probabilities used to fix each
position of the motif. Higher motif scores typically
correspond to motifs that are more statistically
significant as well as more specific }
\item{fg.matches}{Frequency of sequences matching this
motif in the foreground set} \item{fg.size}{Total number of
foreground sequences} \item{bg.matches}{Frequency of
sequences matching this motif in the background set}
\item{bg.size}{Total number of background sequences}
\item{fold.increase}{An indicator of the enrichment level
of the extracted motifs. Specifically, it is calculated as
(foreground matches/foreground size)/(background
matches/background size).} }
}
\description{
Find overrepresented sequence motifs
}
\examples{
# Get paths to sample files
fg.path = system.file("extdata", "fg-data-ck2.txt", package="rmotifx")
bg.path = system.file("extdata", "bg-data-serine.txt", package="rmotifx")

# Read in sequences
fg.seqs = readLines( fg.path )
bg.seqs = readLines( bg.path )

# Find overrepresented motifs
mot = motifx(fg.seqs, bg.seqs, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)

# View results
head(mot)
}


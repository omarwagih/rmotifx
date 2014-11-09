
<img style="display: block; margin-left: auto; margin-right: auto" src="inst/extdata/motifxlogowhite.gif">

------------

# Introduction
This package contains a useable implementation the motif-x tool in the R programming language. motif-x aims at extracting overrepresented sequence motifs. For more information, please refer to the original [motif-x resource](http://motif-x.med.harvard.edu/). The current implementation only supports sequences with a fixed length (i.e. pre-aligned) and have a fixed central residue. For example, phosphorylation sites. 

# How to install?
The motif-x R package can be directly installed from github. First ensure the `devotools` package is installed. Then run the following:

    require(devotools)
    install_github('motifx', 'omarwagih')
    
# How to use?
The package contains the function `motifx` which does everything. For a basic run, you will need a foreground and background set of sequences. 

We can go ahead and use the examples provided with the package: 
     
    # Get paths to sample files
    fg.path = system.file("extdata", "fg-data-ck2.txt", package="motifx")
    bg.path = system.file("extdata", "bg-data-serine.txt", package="motifx")
    
    # Read in sequences
    fg.seqs = readLines( fg.path )
    bg.seqs = readLines( bg.path )
    
    # You can take a look at the format of the sample data
    head(fg.seqs)
    head(bg.seqs)
    
The foreground data here represents phosphorylation binding sites of Casein Kinase 2. The negative data represents 10,000 random serine-centered sites.

To start the program, run the following:

    # mot = motifx(fg.seqs, bg.seqs, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
    # head(mot)

That's it!

For detailed explanations of all parameters and output, check out the documentation by typing `?motifx`. You can also refer to the original motif-x [resource](http://motif-x.med.harvard.edu/motif-x.html) or [paper](http://motif-x.med.harvard.edu/publications/Chou_Schwartz_motif-x_CPBI_2011.pdf). 


# Todo

- Add support for degenerate motifs
- Add support for DNA sequences. Currently, only protein supported.
- Allow motif discovery in non-centered k-mers




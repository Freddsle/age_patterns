install.packages("devtools")
library(devtools)
install_github('arendsee/rhmmer')
install_github("arendsee/phylostratr")

#########

library(phylostratr)
library(reshape2)
library(taxizedb)
library(dplyr)
library(readr)
library(magrittr)


focal_taxid <- '4932'

saccharomyces <- Strata(
  tree = ape::read.tree(system.file('extdata', 'yeast', 'tree', package='phylostratr')),
  data = list(faa=list(
    Saccharomyces_cerevisiae   = 'yeast/cerevisiae.faa',
    Saccharomyces_paradoxus    = 'yeast/paradoxus.faa',
    Saccharomyces_mikatae      = 'yeast/mikatae.faa',
    Saccharomyces_kudriavzevii = 'yeast/kudriavzevii.faa',
    Saccharomyces_arboricola   = 'yeast/arboricola.faa',
    Saccharomyces_eubayanus    = 'yeast/eubayanus.faa',
    Saccharomyces_uvarum       = 'yeast/uvarum.faa'
  )),
  focal_species = 'Saccharomyces_cerevisiae'
) %>% strata_convert(target='tip', to='id')


weights <- c(
  '1355160' = 0,   # these three are deprecated
  '245562'  = 0,
  '1094981' = 0,
  '284813'  = 1.1, # these are the reference replacements
  '766039'  = 1.1,
  '237561'  = 1.1
)


strata <- focal_taxid %>%
  # Get all UniProt proteomes
  uniprot_strata %>%
  # build a tree of all UniProt genomes
  # Select a diverse subset of 5 or fewer representatives from each stratum.
  # Only do this above the Saccharomyces genus, since we will later replace
  # Saccharomyces with out own tree.
  strata_apply(f=diverse_subtree, n=5, weights=weights) %>%
  # add prokaryote stratum
  use_recommended_prokaryotes %>%
  # download UniProt sequences (this may take 10+ minutes)
  uniprot_fill_strata


strata <- replace_branch(strata, y=saccharomyces, node='4930') 

# Next we run BLAST
strata <- strata %>%
  strata_blast(blast_args = list(nthreads=8)) %>%
  strata_besthits
# 

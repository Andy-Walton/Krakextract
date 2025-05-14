<h1 align="center">Krakextract</h1>

A little R package that extracts species of interest from Kraken reports and downloads good quality reference genomes for them from the NCBI refseq database.

<h2 align="center">Installation</h2>
Install the package using:

```
install.packages("devtools")
devtools::install()
devtools::install_github("Andy-Walton/Krakextract")
```

<h2 align="center">Instructions</h2>
The package contains 3 functions. They should be run in order:

<h3 align="left">1) get_metadata():</h3>
NB: this only needs to be run once for each 'kingdom' when package is first installed. Specify the kingdom (e.g. bacteria or fungi) that you are interested in and it will donwload the NCBI RefSeq metadata file for it.

<h3 align="left">2) parse_report():</h3>
This reads .report files produced by kraken. You specify families of interest (e.g. Plasmodiidae) and it will produce a .txt file with the scientific names of the most abundant species from that family in your sample. Note you can use other taxonomic levels over 'family' but this is experimental.

<h3 align="left">3) download_refs():</h3>
This takes the .txt files produced above and downloads the reference genome (or if not available, most recent complete genome) for each species from NCBI RefSeq.

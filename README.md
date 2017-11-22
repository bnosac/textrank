# Summarize Text by Ranking Sentences and Extracting Keywords

This repository contains an R package which handles summarizing text by using textrank. 

For ranking sentences, this algorithm basically consists of.

- Finding links between sentences by looking for overlapping terminology
- Using Google Pagerank on the sentence network to rank sentences in order of importance

For finding keywords, this algorithm basically consists of.

- Extract words following one another to construct a word network
- Using Google Pagerank on the word network to rank words in order of importance
- Constructing keywords - which are the combination of relevant words identified by the Pagerank algorithm which follow each other


## Installation & License

The package is available under the Mozilla Public License Version 2.0.
Installation can be done as follows. Please visit the package documentation and package vignette for further details.

```
install.packages("textrank")
vignette("textrank", package = "textrank")
```

For installing the development version of this package: `devtools::install_github("bnosac/textrank", build_vignettes = TRUE)`

## Support in text mining

Need support in text mining?
Contact BNOSAC: http://www.bnosac.be


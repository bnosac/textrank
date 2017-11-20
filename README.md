# Summarize Text by Ranking Sentences

This repository contains an R package which handles summarizing text by using textrank. This algorithm basically consists of.

- Finding links between sentences by looking for overlapping terminology
- Using Google Pagerank on the sentence network to rank sentences in order of importance

## Installation & License

The package is availabe under the Mozilla Public License Version 2.0.
Installation can be done as follows. Please visit the package documentation and package vignette for further details.

```
install.packages("textrank")
vignette("textrank", package = "textrank")
```

For installing the development version of this package: `devtools::install_github("bnosac/textrank", build_vignettes = TRUE)`

## Support in text mining

Need support in text mining?
Contact BNOSAC: http://www.bnosac.be


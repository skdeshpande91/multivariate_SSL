## mSSL: The Multivariate Spike-and-Slab LASSO

[Sameer K. Deshpande](http://people.csail.mit.edu/sameerd/)

An R package for simultaneous variable and covariance selection in multivariate linear regression with spike-and-slab LASSO (SSL) regularization.
For further details about the method, please see [Deshpande, Rockova, and George (2018)[https://arxiv.org/abs/1708.08911].


### Details

#### Installation

The package source files are contained in the sub-directory mSSL/.
To install, you can either download that directory and then build and install the package from the command line (e.g. `R CMD BUILD mSSL` followed by `R CMD INSTALL mSSL_1.0.tar.gz`).
You can also install the package using `devtools::install_github` as follows.

```r
library(devtools)
devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")
```

#### Examples

The sub-directory examples/ contains R scripts to replicate some of the examples in our paper.


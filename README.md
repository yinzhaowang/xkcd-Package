# xkcd-Package

This is a R package of xkcd distribution. The structure is similar to that of rnorm, dnorm, pnorm, qnorm in 'stats' package.

In this package, we have corresponding functions rxkcd, dxkcd, pxkcd, qxkcd.

Folder named xkcd is pure R. Folder named xkcd_c is combination of R and C.

## Notice: if you want to R CMD check xkcd_c, remember to change the folder name to 'xkcd' before checking, since the NAMESPACE file is written as useDynLib(xkcd)

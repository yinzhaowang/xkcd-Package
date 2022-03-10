pkgname <- "xkcd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('xkcd')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("xkcd")
### * xkcd

flush(stderr()); flush(stdout())

### Name: xkcd
### Title: xkcd
### Aliases: xkcd dxkcd pxkcd qxkcd rxkcd

### ** Examples

require(graphics)

dxkcd(0) == Inf
dxkcd(1/sqrt(2*pi)) == 0


## Using "log = TRUE" for an extended range :
par(mfrow = c(2,1))
plot(function(x) dxkcd(x, log = TRUE), 0.001, 1/sqrt(2*pi),
     main = "log { xkcd density }")
curve(log(dxkcd(x)), add = TRUE, col = "red", lwd = 2)
mtext("dxkcd(x, log=TRUE)", adj = 0)
mtext("log(dxkcd(x))", col = "red", adj = 1)


plot(function(x) pxkcd(x, log.p = TRUE), 0.001, 1/sqrt(2*pi),
     main = "log { xkcd Cumulative }")
curve(log(pxkcd(x)), add = TRUE, col = "red", lwd = 2)
mtext("pxkcd(x, log=TRUE)", adj = 0)
mtext("log(pxkcd(x))", col = "red", adj = 1)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

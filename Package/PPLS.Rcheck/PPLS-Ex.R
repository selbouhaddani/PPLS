pkgname <- "PPLS"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "PPLS-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('PPLS')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PPLS")
### * PPLS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PPLS
### Title: Performs PPLS fit.
### Aliases: PPLS

### ** Examples

exX = scale(matrix(rnorm(100*10),100,10))
exY = scale(matrix(rnorm(100*12),100,12))
PPLS(X = exX, Y = exY, nr_comp = 3, EMsteps = 1e4)
PPLS(X = exX, Y = exY, nr_comp = 3, EMsteps = 1e4, initialGuess = "random")

# devtools::install_github("selbouhaddani/O2PLS")
if(require(O2PLS)){
 exinitGuess = list(W = orth(1:10), C = orth(1:12), B = 0.1,
                     sigE = 1,sigF = 1,sigH = 1,sigT = 0.1)
 PPLS(X = exX, Y = exY, nr_comp = 1, EMsteps = 1e4,
       initialGuess = "custom", customGuess = exinitGuess)
}
exconstraints = list(fconstraint(list(B = 1)) , fconstraint(list(sigT = 1)))
PPLS(X = exX, Y = exY, nr_comp = 2, EMsteps = 1e4, constraints = exconstraints)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PPLS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fconstraint")
### * fconstraint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fconstraint
### Title: Defines an appropriate list with chosen constraints
### Aliases: fconstraint

### ** Examples

fconstraint(list(sigH=.1))
fconstraint(list(B=3,sigE=.5))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fconstraint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

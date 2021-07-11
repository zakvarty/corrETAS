# corrETAS
R package for simulating and fitting the Epidemic Type Aftershock Sequence (ETAS) model with correlated magnitudes.

## Installation
I haven't yet set things up nicely for you to be able to download the package directly from Github using `devtools`.

```
  devtools::intall_github("zakvarty/corrETAS")
```

Currently to install the package you have to: 
  1. download the source (this repo),
  2. open the project in RStudio 
  3. then build and install it yourself using (CMD + SHIFT + B)

Note that this will require the devtools Rpackage in addition to one of Rtools, Xcode or r-base-dev depending on your operating system. See the relevant [R packages book section](https://r-pkgs.org/setup.html#setup-tools) for more details.

Any other viable method to build and install an R package from source should also work here (e.g. R CMD build at the command line.)
for detailed instructions on building and installing R packages see [Karl Broman's guide](https://kbroman.org/pkg_primer/pages/build.html) or [Hadley Wickham and Jenny Bryan's book](https://rpkgs.org). 

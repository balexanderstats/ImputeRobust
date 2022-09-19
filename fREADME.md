# ImputeRobust
Modification of ImputeRobust R package to remove the dependency on extremevalues. The extremevalue package depends on the gWidgets2tcltk R package
which is not compatiable with headless RStudio servers (https://community.rstudio.com/t/problem-installing-gwidgets2tcltk/137978.  
The necessary functions from extremevalues are now included in the package.

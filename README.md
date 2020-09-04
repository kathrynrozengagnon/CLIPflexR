Readme
------

**CLIPflexR** is a R wrapper for the [**CLIP Tool Kit
(CTK)**](https://zhanglab.c2b2.columbia.edu/index.php/CTK_Documentation)
and additional functions to call other external libraries into an R
environment.

**CLIPflexR** makes use of the **CondaSysReqs** library to install conda
dependencies to your machine within a conda environment using the
reticulate libraries.

The [**CTK
toolkit**](https://zhanglab.c2b2.columbia.edu/index.php/CTK_Documentation)
is installed as part of the **CLIPflexR** **install\_ctk()** function.

### How to install CLIPflexR

    install.packages("devtools")
    devtools::install_github("kathrynrozengagnon/CLIPflexR")
    devtools::install_github("RockefellerUniversity/CondaSysReqs")
    CondaSysReqs::install_CondaSysReqs("CLIPflexR")
    CLIPflexR::install_ctk()

Path to Conda tools and Ctk/czplib
----------------------------------

    library(CLIPflexR)

    ## CLIPflexR_0.1.18 conda env found at /Users/kathryn/Library/r-miniconda/envs/CLIPflexR_0.1.18

    ## ctk found  at /Users/kathryn/Library/r-miniconda/envs/CLIPflexR_0.1.18/bin/ctk

    ## czplib found  at /Users/kathryn/Library/r-miniconda/envs/CLIPflexR_0.1.18/lib/czplib

    getOption("CLIPflexR.condaEnv")

    ## [1] "/Users/kathryn/Library/r-miniconda/envs/CLIPflexR_0.1.18"

    getOption("CLIPflexR.ctk")

    ## [1] "/Users/kathryn/Library/r-miniconda/envs/CLIPflexR_0.1.18/bin/ctk"

    getOption("CLIPflexR.czplib")

    ## [1] "/Users/kathryn/Library/r-miniconda/envs/CLIPflexR_0.1.18/lib/czplib"

See our vignettes for [installing CLIPflexR](installCliPR.html) processing CLIP data using [CTK](StandardandBrdU_Processing_CTK.html) or [CLIPflexR](Processing_to_matrix.html), or mix and match!
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Please bear with us while we update examples and vignettes, coming soon....
---------------------------------------------------------------------------

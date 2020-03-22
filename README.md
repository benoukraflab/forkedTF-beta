# forkedTF (ver. 0.1.0)
forkedTF builds a graph model of multiple PWMs (or forked-logos) to have a better representation of the different motif a TF uses when binding in complex with other transcription factors.
### Installation
Type the following in R console:
```r
install.packages("devtools")
library("devtools")
devtools::install_github("https://github.com/benoukraflab/forkedTF")
```
#### Dependencies
   - [TFregulomeR](https://github.com/benoukraflab/TFregulomeR) (>= 1.2)
   - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) (>= 3.2.1)
   - [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) (>= 2.3)
   - [gridGraphics](https://cran.r-project.org/web/packages/gridGraphics/index.html) (>= 0.4.1)
   - [ggseqlogo](https://cran.r-project.org/web/packages/ggseqlogo/index.html) (>= 0.1)
   - [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html) (>= 1.0.0)
   - [grid](https://www.rdocumentation.org/packages/grid/versions/3.6.1) (>= 3.6.1)

## Documentation
For an extended view of the options please visit the documentation in our [vignette](http://benoukraf-lab.com/vignettes/vignette-FPWM.html)

## Quick usage
#### MiniCofactor Report
This function helps in the exploration of binding partners in a cell line. Input the `mainTF`, `cell` line of interest and the `cobinding_threshold` to generate a PDF report of the binding partners.

```r
library("forkedTF")
miniCofactorReport( TF = "CEBPB", cell = "K562",cobinding_threshold=.06)
```
<div align="center">
<a name="miniCofactorReport"/>
<img src="./inst/MM1_HSA_K562_CEBPB_cofactor_minireport.png" alt="miniCofactorReport" width="490" height="630" ></img>
</a>
</div>

#### FPWM creation and plot
Use the `createFPWM` function to extract the motif, from empirical datasets, that a TF uses when binding with a partner TF. `plotFPWM` helps in visualizing the FPWM.
```r
fpwm <- createFPWM(mainTF ="CEBPB",
                        partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                        cell = "K562", 
                        forkPosition = 5,
                        flipMatrix = FALSE)

plotFPWM(fpwm,pdfName="fpwm_plot.pdf")
```
<div align="center">
<a name="fpwm_plot"/>
<img src="./inst/fpwm_plot.png" alt="fpwm_plot" width="340" height="630" ></img>
</a>
</div>

#### Writing FPWM
Save the FPWM to a local file can be used in matrix scanning or matrix clustering in `transfact` format or `FPWMtransfact` format. `Transfact` format will have a matrix for each interacting partner in the FPWM, while `FPWMtransfact` will output a single matrix.
```r
write.FPWM(FPWM = fpwm, format = "transfac", fileName = "FPWM.transfact" )
write.FPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "FPWM.FPWMtransfac" )
```

# Empirical analysis - in server

This part of the analysis was performed in Denmark Statistics' secure servers, using confidential Danish register data. As such, this part of the code is not immediately replicable. However, our access is not unique and others can gain similar access by following a procedure described by Statistics Denmark at [http://www.dst.dk/en/TilSalg/Forskningsservice.aspx](http://www.dst.dk/en/TilSalg/Forskningsservice.aspx). For future reference, the number of our project, approved by Denmark STatistics, was 705109.

The do files perform three tasks.

1. **Data build:** (`deathbuild.do` and `rawbuild.do`). `rawbuild.do` combines the many different raw administrative data files into a tabular dataset with a format suitable for analysis. It outputs  `workdata.dta`, the dataset used for our empirical results. `deathbuild.do` constructs an aggregated dataset on the whole population we use to calculate death probabilities and other population  moments necessary for the calibration of the structural models.

2. **Analysis:** (`main_analysis.do`, `descriptives.do`, and `moments.do`). The main file is `main_analysis.do`. This file performs all regressions for the empirical results of the paper and robustness checks, and outputs the relative `.ster` (stata results) files in server folder. We could not get the permission of downloading all `.ster` files from the server (by the end of the revision process, there were hundreds of them). We only downloaded the ones relative to our main results and moments that we use to calibrate the structural model. These are in the `./estimates/` folder. `descriptives.do` simply exports a few descriptives that we report in the paper. `moments.do` exports the other (descriptive) moments used in the calibration of the model. All these files are originally save in specific folders in the remote server. In the current project, we have saved them under `model/data/*.csv` for convenience.

3. **Output:** (`tables.do`, `graphs.do`, `moments.do`). These files take the inputs provided by steps 1 and 2 and outputs results in the form that either directly appear in the paper (and can be found in this repository in the form of graphs or tex files) or are used as estimation moments in the remainder of the paper.

The whole suite could be executed by the `master.do` file, were the raw administrative datasets available.

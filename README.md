# isolate-chem-mixtures

Data and code relating to our paper *High-throughput characterization of bacterial responses to complex mixtures of chemical pollutants* in Nature Microbiology: https://doi.org/10.1038/s41564-024-01626-9 

## Background

We created all possible mixtures of 8 chemicals (256 mixtures, including a control with no chemicals) and tested the growth of bacterial isolates in those chemical mixtures.

The chemical stressors for this initial work were designed to represent different types of chemical stress entering water ways in urban and rural environments. For each environment there is an antibiotic, a fungicide, a herbicide and an insecticide.

| Chemical      | Environment | Target  |
| ----------- | ----------- | --------- |
| Amoxicillin | Urban | Antibiotic |
| Chlorothalonil | Rural | Fungicide |
| Diflufenican | Rural | Herbicide |
| Glyphosate | Urban | Herbicide |
| Imidacloprid | Urban | Insecticide/Nematicide |
| Metaldehyde | Rural | Insecticide/Molluscicide |
| Oxytetracycline | Rural | Antibiotic |
| Tebuconazole | Urban | Fungicide |

The growth of bacteria was assayed by optical density in the presence of these chemicals. The bacterial isolates were from pristine stream systems in Iceland – i.e. bacteria not expected to have any history of being exposed to chemical stress before. Two "lab" bacteria were also tested for comparison - E. coli and Aliivibrio fischeri. A. fischeri is a luminescent marine bacterium for which a luminescence toxicology assay has been developed.

12 Strains of bacteria used, plus a mixture of isolates:

* 74	- Neobacillus soli (Firmicutes)
* 100	- Pseudomonas baetica (Gammaproteobacteria)
* 302	- Flavobacterium glaciei (Bacteroidetes)
* 306	- Arthrobacter humicola (Actinobacteria)
* 331	- Pseudomonas baetica (Gammaproteobacteria)
* 371	- Rhizobium herbae (Alphaproteobacteria)
* 419	- Sphingomonas faeni (Alphaproteobacteria)
* 448	- Carnobacterium gallinarum (Firmicutes)
* 487	- Aeromonas popoffii (Gammaproteobacteria)
* 527	- Arthrobacter humicola (Actinobacteria)
* E. coli 	(Gammaproteobacteria)
* A. fischeri	(Gammaproteobacteria)
* Mixture - mixture of of all pristine strains (i.e. excluding E. coli and A. fischeri)

Bacteria were diluted 1:100 from carrying capacity in LB media supplemented with chemical stressor additions at 0.1mg/L. Absorbance at 600nm (OD600) was then recorded for 72hrs at 1hr interval timepoints to capture the entire growth curves for these bacteria. Each bacteriaXchemical treatment has 4 replicates.

## /data/

* The raw data consists of 20 text files for each bacterial strain tested. 256 chemical combinations were assayed, using the internal 60 wells of 96-well plates (i.e. excluding Rows A and H and columns 1 and 12). This requires 5 plates, which were replicated 4 times, for a total of 20 plates - one file for each plate. The files record the OD600 in each well at each timepoint. These are named in the following system: isolate-\<number/name>-plate-\<number>.txt, e.g. isolate-74-plate-1.txt.

* chemical-matrix.csv - matrix of chemicals designating what is in each well.

* tidy-data.csv - single file combining the raw data files after some data cleaning operations (performed by `code/01-clean-data.R`)

## /code/

Data manipulation and analysis scripts:

* 01-clean-data.R - combine raw data files and tidy them up for future use.
* 02-growth-curves.R - calculate some parameters from the growth curves and output these as well as some summary statistics. Currently this is simply fitting a spline curve and integrating to find the area under the curve (AUC) as a measure of total "growth" over the 72hrs.
* 03-bootstrap-interactions.R - use the bootstrapping methodology described in the paper to determine whether chemical responses display interactive effects. Ouputs `results/bootstrapped-emergent-interactions.csv` and `results/bootstrapped-net-interactions.csv` *Note: this is a bootstrapping approach; re-running the code will result in slightly different numbers output.*

Post-analysis and results visualisation scripts:

* 04-plotting-heatmap.R - Create the aligned heatmaps in Figure 1A.
* 05-basic-visualisations.R - Create most of the remaining analysis figures, and produce statistics for the main text.
* 06-interaction-networks.R - Produce the visualisations of the interactions as networks, for figure 5 and supplementary figure S6.
* 07-phylogenetic-tests.R - Perform phylogenetic analyses including mantel test and lambda and K tests for phylogenetic signal. 


## /results/

* spline_fits.csv - AUCs from fitting spline curves to each bacteria/stressor combination.
* model_summaries.csv - AUCs from replicate wells summarised into a mean and SD.
* bootstrapped-emergent-interactions.csv - results of bootstrapping tests for emergent interactions.
* bootstrapped-net-interactions.csv - results of bootstrapping tests for net interactions.
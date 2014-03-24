# Modular color evolution facilitated by a complex nanostructure

Eliason, C. M., Maia, R. M., Shawkey, M. D.

## Introduction

<!-- Complex traits --> Evolution, either by selection or drift, requires heritable phenotypic variation {Futuyma}. Complex morphological traits like the jaws of cichlids {Liem?} or mandibles of mice {ref} (maybe beetle horns?) have numerous parameters (e.g., structural, morphological) and therefore greater potential for phenotypic variation than simple traits with fewer parameters {Vermeij 1973}. <!-- Target of selection -->In nature, selection often acts on functional properties of complex traits rather than directly on morphology. For example, sexual selected song characteristics in reed warblers are functional properties of the vocal apparatus {Hasselquist et al. 1996}. By contrast, sexual selection on tail length in widowbirds acts directly on morphology {Andersson 1992}). With functional traits, the way that form maps to function can be important in understanding how a trait responds to selection (i.e. its evolvability; {Wainwright 2007}).

<!-- Comparing different form-function maps - In simple traits (1 parameter), tradeoffs are present so less evolvable. More complex traits often have complex form-function relationships that can enhance evolvability. You can have a complex morphology, but all produce same function (many-to-one), or all independent functions (many-to-many). Thus, it is important to understand the form-function map to understand how functional properties can evolve. - Proximate studies of the form-function map can help explain differences in functional and morphological diversity among animal groups (e.g., Hulsey & Wainwright {2002} found that morphological and mechanical diversity were decoupled in XX).  -->

<!-- Bird color mechanisms --> Color is a well-studied functional property of the metazoan integument {ref}. Bird are well-known for their diverse plumage colors, and such colors are thought to evolve primarily through sexual selection {Andersson 1994}. <!-- Most studies focus on external forces (selection) in driving divergence. However, internal features of organisms (e.g., developmental constraints, form-function relationships) may also influence phenotypic divergence. Evolution needs variation, thus traits that produce variation should promote evolutionary divergence. BUT mechanism may also be  mportant. -->Feather colors can be produced either by selective light absorption by pigments that birds deposit in their feathers (pigment-based colors) or light scattering from nanoscale variation in refractive index (structural colors). Because color is a multidimensional trait described by three optical parameters (hue, saturation and brightness), an increase in the number of tunable morphological parameters in color-producing nanostructures may increase the number free optical parameters. Structural colors are thus an excellent system for studying how the functional architecture of complex traits influences phenotypic evolvability.

<!-- _We need to answer these questions: Why not just study color? and Who cares about morphology?_ -->

<!-- Hypothesis --> We described a photonic crystal (PC) structure in duck flight feathers that is formed by the arrangement of small (100-200 nm) melanin-containing organelles (melanosomes) into a hexagonal pattern {Eliason 2012}. PCs are therefore complex structures (made up of repeating subunits) that cause fluctuations in refractive index along certain directions. This periodic variation in refractive index, along with the small size of subunits (near the wavelengths of visible light, 400-700 nm), causes certain wavelength of light to be reflected more strongly than others (through constructive interference), producing visible color.<!-- This type of structure is flexible, producing diverse iridescent colors across the bird-visible spectrum through only minor changes in the size or spacing of melanosomes {Eliason 2012}. --> Compared to simple structures like thin-films responsible for the rainbow colors of soap bubbles, <!-- in which only the thickness of the film can vary to vary color {optics ref},  -->PCs in ducks are complex with numerous parameters controlling form (e.g., the size and spacing of melanosomes) {Eliason 2012, Joannopoulos 2008}. <!-- Add in that previous results suggest independent contributions of spacing and density to color properties hue and saturation/brightness? -->Here we ask how this increased morphological complexity affects the independent variation of color traits (evolvability). We hypothesized that color evolvability is enhanced by the way that morphology maps to color. This hypothesis makes the following predictions: i) color variation will be explained by variation in morphology, ii) "modular" (functionally independent) color traits will evolve independently and at different rates, and iii) rates of color evolution parallel those in underlying functional morphology.


## Methods

### Species sampling

Of the 54 recognized species in the sub-family Anatini {IOC ref}, 44 have available DNA sequences. We sampled nearly all of these species (42/44, 95%) at either the Field Museum of Natural History (FMNH) or the University of Michigan Museum of Zoology (UMMZ), with the exception of _Anas bernieri_ and _A. smithii_ (no suitable specimens available at either museum). Prior to analysis, we removed four species with non-iridescent wing patches (_Tachyeres pteneres_, _A. georgica_, _A. sibilatrix_, and _A. strepera_).

### Building a time-proportional phylogeny

To allow for the determination of rates of morphological and color evolution, we estimated a time- calibrated phylogeny for dabbling ducks using published mtDNA sequences (cytochrome b oxidase, NADH dehydrogenase; {refs}) and an uncorrelated relaxed clock algorithm implemented in the program BEAST 1.7.4 {Drummond}. We used available fossil data and exponential priors to place constraints on the root age and divergence times for the most recent common ancestor for six sets of taxa {REFS}. We used a starting tree estimated in MrBayes v. 3.1.2 {ref} with parameters specified in Gonzalez et al. {Gonzalez, 2009}. We ran five separate analyses in BEAST for 10 million generations each, sampling every 10 thousand generations, and checked that these independent analyses converged using Tracer 1.5 {ref}. In addition, we inspected effective sample sizes (>100) and posterior distributions for tree height (what else?). Finally, we computed the maximum clade credibility (MCC) tree, scaled all trees to a total depth of unity, and retained a posterior set of 500 trees (sampled according to their posterior probability) for further analyses.

### Measuring barbule nanomorphology (form)

In a photonic crystal, the features responsible for producing color in visible wavelengths are at the nanoscale, thus we used transmission electron microscopy (TEM) to view them (see ESM for details on standard lab protocol). We then used the program ImageJ {ref} to measure the following traits for 1-3 barbules per species: melanosome diameter (_d_), separation between melanosomes (_t_), thickness of the outer keratin cortex (_c_), and number of layers in a melanosome stack measured perpendicular to the barbule surface (_l_; see Fig. X for measurement schematic). We haphazardly chose 10 regions
within each barbule to measure to account for potential non-independence in the sizes of adjacent melanosomes. We natural log-transformed all morphological variables and computed species means. This transformation allowed us to compare evolutionary rates among traits with different length scales because it gives proportional change in units of _e_ (2.7) {ref}.

### Measuring spectral reflectance (optical function)

Any color can be characterized by three variables: hue (how red or blue a color is), saturation (how pure or vivid a color is), and brightness (how much light is reflected). To quantify color, we measured spectral reflectance of intact wing patches at normal incidence (with light perpendicular to the feather surface) from 300-700 nm using a spectrophotometer and attached Xenon light source (Avantes Inc., Broomfield, CO, USA). We took three measurements per bird from 1-5 individuals per species. Because brightness is highly sensitive to variation in alignment of the spectrophotometer (ref), we retained the brightest of the three spectra per bird for further analyses {Meadows?}. We then smoothed the spectra with local regression (LOESS) in `R` to minimize electrical noise from the spectrophotometer and thereby increase the accuracy of color variable determination. From these processed spectra, we determined hue as the wavelength of a reflectance peak, saturation as the half-width of a reflectance peak at 50% of peak reflectance (HWHM, see SI for details), and brightness as the height of the main reflectance peak. Saturation was calculated at the midpoint of the difference between the minimum and maximum reflectance values. We used half- width of the peak rather than full width because some of the peaks were near the limits of the wavelength range (i.e. 700 nm). These three color variables were chosen because they allow for direct comparison with optical model predictions in photonic crystals {Eliason 2012; Joannopoulos 2008?}. All color analyses were performed in the `R` package `pavo` {Maia 2013}.

### Linking form to function

To test our hypothesis that the form-function map biases the direction and rate of color evolution, we quantified the relationship between feather morphology and color using phylogenetic generalized least squares (PGLS) multiple regression. PGLS allowed us to account for phylogenetic effects on trait covariation (?) under an OU process using the `corMartins` function in the `geiger` package {Martins and Hansen 1997}. To account for phylogenetic uncertainty in parameter estimates, we ran the full additive PGLS model for each color variable on all 500 posterior trees and retained parameter estimates as well as their standard deviations. Uncertainty in parameter estimates is the combined effect of phylogenetic uncertainty ($\sigma^2_\text{phy}$) and uncertainty in parameter estimation ($\sigma^2_\text{est}$). We computed $\sigma^2_\text{phy}$ as the variance in parameter estimates across trees and $\sigma^2_\text{est}$ by averaging the variance in parameter estimates across all trees. Finally, we summed these two variance components {cite Revell? or you Raf?} and tested for significance of each slope by computing the _t_ statistic (_t_ = slope/se) using a _t_ distribution with 37 degrees of freedom (# of tips in the phylogeny - 1). To test whether the model fit the data, we calculated the difference in AIC values between the full additive model and an intercept-only model across the tree block. All PGLS analyses were done in `R` using the packages `nlme` {ref} and `ape` {Paradis 2008?}. We used multivariate Q-Q plots to assess normality and transformed variables when necessary (see ESM for details). Variables were centered and scaled by their standard deviations prior to analysis so that the standardized regression coefficients could be readily compared among traits {Schielzeth 2010}.

### Determining the mode of character evolution

To explore the evolutionary mode of color and morphology traits, we fit data for optical and morphological traits to a Brownian motion (BM) and Ornstein-Uhlenbeck (OU) model of character evolution using maximum likelihood (ML) in the `R` package `ouch` {Butler 2004}. To test the power to differentiate BM and OU models, as well as univariate and multivariate OU models, we used a phylogenetic Monte Carlo (pmc) approach {Boettiger 2012}. Briefly, following Boettiger {2012}, we did the following: i) estimated parameters with maximum likelihood (ML) for each posterior tree in `ouch`, ii) used these estimates to simulate evolution along a phylogeny under the null and test model for 1000 randomly chosen posterior trees, iii) re-estimated parameters with ML for each simulated dataset and under each model, and iv) computed the simulated
deviance as $\delta = 2*(\text{log}\mathcal{L}_\text{null} - \text{log}\mathcal{L}_\text{test})$. Comparing the overlap in the distribution of deviance gives the power to distinguish the two models, and comparing the observed deviance to the null distribution gives the probability of observing a value more extreme (i.e. the _p_-value) {Boettiger 2012}.

### Measuring rates of character evolution

Evolutionary rates describe how quickly variation accumulates over a given time period {Martins 1994}. For a simple Brownian motion model, variation increases continually with time. By contrast, in an Ornstein-Uhlenbeck model, variation increases at first, but over longer time periods variation reaches a plateau (i.e. the evolutionary rate slows down) as variance-generating mechanisms (drift, mutation) balance out variance-restraining forces (selection, developmental constraints) {Hansen 1997, Martins et al. 2002, Hunt 2012?}.<!-- For a single trait evolving by an OU process, the phenotypic variance at time _t_ is: $V_{OU} = \frac{\sigma^2}{2 \alpha} (1 - e^{-2 \alpha t})$ {Martins 1994}. For very small values of alpha (i.e. weak selection), the equation collapses into the familiar Brownian motion model: $V_{BM} = \sigma^2 t$. The rate of Brownian motion evolution is the slope of the graph of variance versus time.--> For this reason, it is impossible to calculate rates of change at specific times in the past given only data for extant taxa (ref?). However, it is possible to infer an average rate of evolution from the variance that has accumulated up to the present time. <!-- assuming a long enough time has passed, the rate of change will be zero and the stationary covariance matrix __V__  --> To do this, we used published equations {Bartoczek et al. 2012, Hansen and Martins 1996} to calculate the evolutionary variance-covariance matrix __V__ under an OU process (see equation B.8 in Bartoczek 2012). <!--  an estimate of how much variation has accumulated over the entire duration of the clade, from t = 0 to 1 (the scaled depth of the tree). -->We then used __V__ to compare the rates at which different traits have evolved (diagonals of __V__, rate metric $\omega$ in {Hunt 2012}) and test whether different pairs of traits evolve independently or in a correlated fashion (off-diagonals of __V__). We also simulated whether the estimated values of $\sigma^2$ and $\alpha$ would be expected to restrain further phenotypic evolution.

<!-- $\frac{d\textbf{V}}{dt} = -\textbf{A}\textbf{V}-\textbf{V}\textbf{A} + \textbf{S} = 0$ -->

<!-- $\textbf{V}(t) = \int_{t=0}^1 e^{-\textbf{A}t}\textbf{S}e^{-\textbf{A}^\text{T} t}dt$ -->

<!-- where __A__ and __S__ are the matrices of $\alpha$ and $\sigma^2$ values, respectively. Thus, as an estimate of how much variation evolved over the entire (time of the clade?), we computed __V__(_t_) after equation B.8 in Bartoczek et al. (2012).  -->

To calculate confidence intervals for the elements of __V__, we used the pmc approach described above to simulate evolution under the best fitting OU model. To test for correlated evolution, we evaluated whether the 95% confidence intervals of trait covariance estimates overlapped zero. To test whether traits evolve at significantly different rates, we calculated pairwise differences in evolutionary rate for all combinations of color and morphology traits. We then applied sequential Bonferroni correction to account for multiple tests {REF?} and assessed significance by examining whether the 95% confidence intervals overlapped zero. We determined the number of pmc simulations needed to capture variation in rate difference estimates using a rarefaction of variance analysis {see Claramunt 2010 for `R` code}. Briefly, we randomly select samples of a given size (ranging from 0-max MC sample size, in increments of 20), repeated this for 100 replicates for each sample size, calculated mean SD across replicates, and finally plotted mean SD versus sample size. A leveling off of the resulting rarefaction curve indicates increasing sample size will have little effect on the captured variation.

### Incorporating measurement error

Comparisons of evolutionary rates among different phenotypic traits are prone to measurement error (ME) {Ives, Revell, Adams}. Because TEM is both time and cost prohibitive, it is unrealistic to examine barbule morphology for many individuals in all species needed to compute ME. However, an estimate for intraspecific variation of color and morphological traits was available for a common species (mallard, _Anas platyrhynchos_; unpublished data). Thus, to account for ME we re-ran the above analyses incorporating these estimates of error (see {Harmon et al. 2010} for a similar approach). Although this assumes similar ME for all species (likely unreasonable), calculating rates with and without ME allowed us to evaluate the relative sensitivity of different traits to ME (i.e. slopes of evolutionary rate for different MEs), and therefore whether our results are robust.

<!-- We also ran analyses at 0.5X and 2X the observed ME to see how variation in ME would influence our results (this is what Luke Harmon et al. did in their 2010 paper). -->

### Outlier sensitivity analyses

Bivariate plots of morphological and optical variables revealed two potential outlier species: _A. aucklandica_ (very broad reflectance peak) and _A. rubripes_ (smallest melanosome diameter and spacing). To assess whether these values affected our analyses, we pruned them from the phylogenies and re-ran all analyses.


## Results

### Color variation is explained by variation in morphology

TEM analysis showed that nearly all dabbling duck species sampled (38/44, 86%) had a hexagonal PC within feather barbules. <!-- We removed the four non-iridescent species (_A. georgica_, _A. sibilatrix_, _A. strepera_, and _Tachyeres pteneres_) from our dataset before running comparative analyses. --> Melanosome diameter ranged from 90 - 180 nm, spacing between adjacent melanosomes varied from 18 - 87 nm, and the number of layers at the surface of barbules ranged from 4-8. Cortex thickness ranged from 110 - 478 nm. Peak location (hue) ranged from blue-red (463 - 647 nm), peak width (saturation) ranged from 32 - 69 nm, and peak reflectance (brightness) ranged from 5 - 70%. <!-- Unless otherwise indicated, all values are mean ± 1 S.E. of the mean. -->

Melanosome diameter had a significant positive effect on hue (phylogenetic generalized least squares regression, PGLS: $\beta$ = 0.36 ± 0.18, _t_ = 2.07, _p_ = 0.045) and a significant negative effect on saturation ($\beta = -0.50$ ± 0.18, _t_ = -2.77, _p_ = 0.0089) but did not significantly influence brightness (PGLS, $\beta$ = -0.35 ± 0.20, _t_ = -1.77, _p_ = 0.085; see Fig. X). Melanosome spacing had a significant positive effect on brightness (PGLS, $\beta$ = 0.43 ± 0.16, _t_ = 2.69, _p_ = 0.011) but did not significantly predict hue ($\beta$ = -0.064 ± 0.15, _t_ = -0.42, _p_ = 0.67) or saturation ($\beta$ = 0.22 ± -0.16, _t_ = 1.37, _p_ = 0.18; Fig. X). The number of melanosome layers did not significantly predict any color variables.

### "Modular" color traits evolve independently and at different rates

<!-- ### What is the mode of evolution in color and morphology? -->

<!-- Scatter plots showed weak relationships between morphological and optical traits (Fig. X). -->

The multivariate OU model was strongly preferred over the Brownian motion (BM) model for both morphological and optical traits (both _p_ < 0.001). `pmc` simulations showed that there was adequate power to tell these models apart, as none of the BM simulations had a log likelihood greater than the alternative OU simulations (see Fig. X). We further compared multivariate and univariate OU models (i.e. off-diagonal elements of $\sigma^2$ and $\alpha$ set equal to zero) and found strong support for the multivariate model for both suites of traits (both _p_ < 0.001; power = 99.3% and 98.0% for color and morphological traits, respectively). Figure X shows that ~5000(?) simulations were enough to capture the variation in estimates for both diagonal (Fig. Xa) and off-diagonal (Fig. Xb) components of __V__. Saturation was positively correlated with hue (_r_ = 0.43, CI: 0.08 - 0.70) and negatively correlated with brightness (_r_ = -0.43, CI: -0.78 - -0.11). Brightness and hue were not significantly correlated (_r_ = 0.18, 95% CI: -0.18 - 0.52). No morphological traits were significantly correlated, as all 95% confidence intervals overlapped zero (see Fig. x). Melanosome spacing evolved significantly faster (~5.5X) than the size of melanosomes or the number of layers, but the latter two traits did not evolve at significantly different rates (ratio = 1.9, see Fig. xx). Brightness evolved significantly faster than either hue (~73X) or saturation (~17X), and saturation evolved ~4X faster than hue (Fig. xx).

### Rates of evolution in color parallel those in morphology

Patterns of disparity among morphological traits paralleled those in color. Specifically, melanosome spacing (and consequently brightness) evolved faster than melanosome diameter (and consequently hue and saturation). Hue and saturation were positively correlated. This is likely because both of these variables stem from a common underlying morphological trait (the diameter of melanosomes), thus changes in diameter would produce correlated changes in either color trait. Hue and brightness evolved independently, as expected from their functional independence (see our PGLS results). Saturation and brightness showed a pattern of correlated evolution. Saturation decreased with melanosome diameter (reflectance peaks were broader for large melanosomes). Although the only significant predictor of brightness was the spacing between melanosomes ($\beta$ = 0.43 ± 0.16, _t_ = 2.69, _p_ = 0.01), the relationship between diameter and brightness was marginally significant ($\beta$ = -0.35 ± 0.20, _t_ = -1.77, _p_ = 0.08) and expected from theory (larger melanosomes should produce broad reflectance peaks and bright colors). Thus, the observed significant evolutionary correlation between these color traits might stem from a common morphological basis. Further work could explore this idea by experimentally manipulating the distance between melanosomes (during development?). Another possible explanation for the correlation could be disorder (more disordered melanosomes might cause duller colors and broader reflectance peaks). However, we did not find a significant relationship between the coefficient of variation in melanosome size or spacing (a measure of disorder) and any color variables (add results?).

### Effects of outliers and measurement error

When we removed _A. aucklandica_, the effect of diameter on hue was only marginally significant ($\beta$ = 0.34, _t_ = 1.93, _p_ = 0.06), the correlation between diameter and spacing was significantly positive (_r_ = Xxx), and saturation evolved slower (0.0085 _e_/time compared to 0.016 _e_/time). When we removed _A. rubripes_, the correlation between diameter and spacing was far from significance (_r_ = 0.13 [-0.23, 0.46] compared to 0.33 in the full analysis) and brightness evolved faster than spacing. In all cases, spacing evolved faster than diameter, diameter and saturation evolved at same rates, hue evolved slower than diameter, there was a fairly strong but non-signifant effect of diameter on brightness ($\beta_\text{auck}$ = -0.33, $\beta_\text{rubr}$ = -0.24), and no off-diagonal values for $\alpha$ or $\sigma$ were significantly different from zero.

As expected, incorporating ME decreased the estimated evolutionary rates for all traits. Morphological traits were more affected than optical traits because of their increased standard errors (see Table X). Indeed for one variable (_d_), the within-species estimate of SD was greater than the between-species variation in the trait. The ratio between intraspecific and interspecific variation determined the effect of ME on evolutionary rate {see Lynch 1990, eq. 3}. Including ME did not significantly alter the pattern of rate differences among morphological and optical traits. Melanosome spacing evolved significantly faster than diameter or the number of layers in both cases, and brightness evolved significantly faster than either hue or saturation. However, after accounting for ME the evolutionary rate for diameter was not significantly higher than hue.


## Discussion

### Overview

To our knowledge, this study represents the most dense TEM sampling of any
animal clade.

### Nanostructural development

Duck melanosomes are very small compared to most avian melanosomes {Li? dataset}, and these small sizes (and correspondingly high aspect ratios, {ref}) may allow them to be packed into hexagonal configurations. This could explain why hexagonal PCs of solid melanosomes are not found in many other taxa (any?) {refs saying this}. Non-close-packed crystals (having space between melanosomes) are hard to produce artificially because they are often unstable {refs}, yet nearly all Anatid species analyzed here (95%?) have this nanostructure. Thus, our results may provide inspiriation for engineering non-close-packed, and thus optimally bright, photonic crystals.

### Morphology predicts color

Patterns of covariation among traits can be due to correlated selection {Brodie 1992}, functional dependence (the __F__ matrix, e.g. {Holzman?}), or correlated response to selection (G matrix, e.g. xx). It is therefore important to study form-function relationship to draw conclusions about the potential roles of selection (agents of selection?) and drift (constraints?) in driving (producing?) patterns of phenotypic variation. Our PGLS results suggest that morphological traits melanosome diameter and spacing are, at least in part, functionally independent (contributing to separate color attributes). That the covariation between hue and brightness was near zero suggests that functional complexity eliminates tradeoffs in these optical traits {e.g., see Holzman 2011, Wainwright 2007}. By contrast, hue and saturation evolve in a correlated fashion and at similar rates because of their common morphological basis (diameter). Interestingly, the number of melanosome layers was not a significant predictor for any color variable. This could be due in part to the strong absorption by melanin in the upper few layers.

### Morphological complexity enhances color evolvability

<!-- ADD MORE...THIS IS THE COOLEST PART ABOUT THE PAPER! --> Thin-films are common in birds {Eliason 2010, Shawkey, other refs?} and have one parameter that can vary (cortex thickness). Thus, variation in thickness will cause correlated changes in different optical traits (hue, saturation and brightness). By contrast, PCs have multiple parameters that can vary to affect color (melanosome size and spacing), and are thus more optically complex than thin-films.

We found that this complexity enhances color evolvability.

That brightness evolved faster than melanosome spacing when incorporating ME could be explained either by overestimation of ME for spacing, or additional parameters important in brightness variation (e.g., barbule density, orientation, or within-barbule color variation). Interestingly, one species (_Anas querquedula_) showed cell-to-cell differences in pigmentation that would clearly reduce the correlation between nanostructure and brightness predicted from a single iridescent cell (see ESM Fig. X).


### Structural traits evolve faster than morphological traits

<!-- Why is phenotypic character A more variable than phenotypic character B?-->

The rate of evolutionary change depends on the strength of selection ($\alpha$) and amount of additive genetic variation ($\sigma^2$) {evo text ref}. We found that brightness evolves faster than hue because spacing evolves faster than diameter. Moreover, this pattern is likely due to greater genetic (or phenotypic) variation (i.e. the __G__ matrix) rather than relaxed selection or constraints (Figure S2). However, the question of why structural traits diverge faster than morphological traits (e.g., see {King & Jukes 1969}) remains unanswered. One possible explanation is that spacing may be determined by more underlying genes (polygenic) and therefore more sensitive to genetic variation {Mather 194x, Simpson Tempo & Mode}. <!-- Additionally, genotype networks make phenotypes more evolvable {Wagner 2012}. Complex gene-phenoype maps (epistasis? robustness?) can cause faster rates of phenotypic change compared to simple maps. --> Indeed, melanosome synthesis and deposition are complex processes involving a number of genes (reviewed in Marks 2001). <!-- For example, melanin synthesis (TRP1, TRP2, Pmel17), melanosome density or spacing (charge on melanosomes? rate of keratin polymerization/deposition-i.e. conversion from gene-protein; Maia 2011), and melanosome size (Pmel17 {ref}, OA; from Marks 2001). --> Moreover, melanosome spacing is likely sensitive to additional nongenetic factors (temperature, rate of melanosome deposition, keratin polymerization rates {Maia 2011, Prum?}) that could increase {==phenotypic variation==}. <!-- Experiments with artificial structures indicate that variation in the size of particles will automatically produce some correlated variation in spacing (Rengajaran 2005). Var(P) = Var(G) + Var(E) + Var(GxE) --> <!-- 3) differences in population variation? --> Complex developmental processes may evolve quickly due to more evolutionary degrees of freedom (steps in the process, genes involved in regulatory processes, etc.). Alternatively, if developmental plasticity in melanosome spacing brings some individuals closer to the 'preferred' color (brighter?), then selection can act on that variation and facilitate genetic change. <!-- Within species, the variation in spacing is greater than the variation in diameter (add figure?) and spacing is more variable than diameter within barbules (unpublished data?), suggesting that the developmental process causes more fluctuations in spacing. Greater population-level variation can cause greater rates of phenotypic evolution. Thus, the evolutionary significance of a given amount of phenotypic change becomes less significant with greater variation within a population {Lynch 1990}. (neutral genetic divergence expectation: Var(between)/t*Var(within) ~ 1) However, we found that differences in rates were still significant after standardizing by the amount of population-level variation in traits. Furthermore, all traits are linear and on the ratio scale, therefore rate differences for natural-log transformed traits are comparable {Gingerich 2009}. --> Finally, effective population size can vary across the genome <!-- (due to linkage dis?) -->, providing another possibility for why melanosome spacing is diverging faster than diameter. Further work should investigate the mechanistic bases for rate differences in complex phenotypic traits.

<!-- MISC: -->

<!-- Lenski 2001 (http://goo.gl/GaNE0L) - phenotypic evolution depends on beneficial mutations (b/c neutral DNA changes aren't seen in the phenotype and deleterious ones would be eliminated by natural selection) -->

<!-- Subramanian & Kumar 2004 (http://goo.gl/Fz3EJv) - showed that more highly expressed genes caused slower rates of protein evolution (I THINK) -->

An interesting question is how ducks get a wide variation in melanosome size: look up papers on organelle size control (fission-budding dynamics-that curr biology paper..)


### Implications

Phenotypic evolvability can be linked to diversification (Rabosky 2013, Maia 2013, others?). Many bird species are distinguishable only by plumage pattern {find ref}, implicating a role of such signals in reproductive isolation {ref}. More variable colors should therefore be correlated with rates of speciation {Maia et al. 2013}, possibly providing an explanation for the wide variation in diversity among Aves {Jetz 2012}. It would be interesting to look at patterns of diversification across Anseriformes and test whether lineages with morphological innovations (like PCs), through their demonstrated increase in color variability (this study), diversify faster than non-structurally (black, brown) colored lineages.

## Figure captions

Figure 1. PGLS results for mcc tree. Plots show conditional plots of color variables versus morphological variables. Relationship between y and x-axis calculated by holding all other variables in multiple regression at their held at median values. Lines are linear fits. Schematics along y-axes illustrate relevant variation in the shape of reflectance spectra.

![](figures/fig1.png?raw=true)

Figure 2. Evolutionary correlations for different pairs of morphological (a) and optical traits (b). Points show mean rate and line segments are 95% confidence intervals. Results are shown with (black) and without measurement error (grey).

![](figures/fig2.png?raw=true)

Figure 3. Phenotypic divergence rates for different morphological (a) and optical traits (b). Points show mean rate and line segments are 95% confidence intervals. Similar colors indicate functionally related traits from PGLS analysis. Results are shown with (bold colors) and without measurement error (faded colors). Different letters indicate significantly different rates (pairwise diffs, corrected for multiple comparisons). Note log-10 scale for y-axis. 'Traitgrams' in lower panels illustrate the relationship between time (y-axis) and evolutionary change in trait values (x-axis). Upper panels are histograms. Note natural-log scale. Scale range set to that of most variable trait to illustrate differences in evolutionary change.

![](figures/fig3.png?raw=true)


## ESM

### Defining evolvability

Phenotypic evolvability, as envisioned by Vermeij, is controlled by the _number_, _range_, and _independence_ of parameters. The phenotypic covariance matrix can provide insight into how independent different morphological traits are {Arnold?}. Similarly, the __F__ matrix (or form-function map) is a mathematical description of the morphological or physiological traits (rows of the matrix) and their relative contributions to different functional traits (columns of the matrix) {Ghalambor 2003, Walker 2007}. "Functional" evolvability can thus be inferred from the covariance matrix of functional traits (__P__), calculated as:

$\textbf{P}=\textbf{F}^\text{T}\textbf{C}\textbf{F}$,

where __F__ is the form-function map and __C__ is the phenotypic covariance matrix {Holzman 2011}. Considering the case where two traits (x and y) influence _only_ functions f1 and f2, respectively,

$\textbf{F}=\begin{bmatrix}
  1 & 0 \\
  0 & 1 \\
\end{bmatrix}$

and, given two morphological traits _x_ and _y_ that are independent (either genetically or developmentally), the off-diagonal elements of the performance covariance matrix (__P__) will be zero, indicating the two functions can evolve independently. Thus, we define evolvability here as the capacity for independent evolution of morphological or optical traits (i.e. off-diagonal elements of __C__ or __P__ equal to zero).

### Assessing multivariate normality

To assess variables for multivariate normality, we produced Q-Q plots of the quantiles of Mahalanobis distance (=d^2) versus those for a chi-squared probability distribution (?). Deviation from a 1:1 line indicates points that are potential multivariate outliers. Based on these plots, we used untransformed values for _t_ and _d_ in regressions involving brightness (B3). Q-Q plots involving saturation (hwhm) revealed a strongly right/positively skewed distribution, thus we transformed this variable as $Y^{-10}$ (following a boxcox transformation).


### Supplemental Tables

__Table S1__. Observed __F__ matrix (SE in parentheses). Matrix values are standardized partial regression coefficients for PGLS multiple regressions of different color traits (columns) on morphological traits (rows). Significant slopes are bolded.

| Trait     | Hue             | Saturation       | Brightness      |
| :-------- | -----------:    | -----------:     | -----------:    |
| Diameter  | __0.36 ± 0.18__ | __-0.50 ± 0.18__ | -0.35 ± 0.20    |
| Spacing   | -0.06 ± 0.15    | 0.22 ± 0.16      | __0.43 ± 0.16__ |
| Layers    | 0.09 ± 0.15     | -0.02 ± 0.17     | -0.21 ± 0.18    |

__Table S2__. Phylogenetic generalized least squares (PGLS) multiple regression model results for $\alpha$ parameter.

| Variable   | $\alpha$         |
| :---       | :---             |
| Hue        | 16.7 (3.1, 35.6) |
| Saturation | 10.2 (6.3, 16.4) |
| Brightness | 13.3 (2.4, 57.8) |


__Table S3__. Mean and confidence intervals for pairwise differences in evolutionary rates among morphological traits (_d_ = diameter, _t_ = spacing, _l_ = layer number). Alpha values adjusted for multiple comparisons.

| Trait pair | Mean                       | Alpha    |
| :------    | :--------                  | -------: |
| d-l        | -0.0154 (-0.0362, 0.0031)  | 0.0500   |
| t-l        | 0.0656 (0.0156, 0.1363)    | 0.0250   |
| d-t        | -0.0810 (-0.1581, -0.0318) | 0.0167   |

__Table S4__. Mean and confidence intervals for pairwise differences in evolutionary rates among color traits. Alpha values adjusted for multiple comparisons.

| Trait pair        | Mean                       | Alpha    |
| :---------------- | :--------                  | :------- |
| hue-sat           | -0.0123 (-0.0214, -0.0051) | 0.0500   |
| sat-bri           | -0.2595 (-0.4446, -0.1292) | 0.0250   |
| hue-bri           | -0.2718 (-0.4739, -0.1336) | 0.0167   |


__Table S5__. Standard deviations of natural-logged trait values measured for n=10 individual male mallards (_Anas platyrhynchos_).

| Trait    | SD(ln(x)) |
| :----    | :-------- |
| Hue      | 0.02422   |
| B3       | 0.34298   |
| hwhm_min | 0.09789   |
| d        | 0.05882   |
| t        | 0.10868   |
| layers   | 0.21903   |


### Supplemental figure captions

Figure S1. Selecting the best evolutionary model. Lines show deviance distributions (probability density) for 1,000 Phylogenetic Monte Carlo simulations under the null (grey shaded regions) and test models (blue shaded regions). (a,b) Brownian motion versus an Ornstein Uhlenbeck (OU) model. (c,d) Univariate OU versus multivariate OU model. Left panels (a,c) are results for morphological traits, right panels (b,d) are those for color traits. Solid vertical lines show observed deviance and dashed vertical lines 5% quantile of null distribution.

![](figures/figs1.png?raw=true)

Figure S2. Estimated values for $\alpha$ and $\sigma^2$. Boxplots show distribution of parameter values for $\alpha$ (red boxes) and $\sigma^2$ (green boxes) for functionally related morphology (a) and color variables (b). Note the different y-axes and log-10 scale.

![](figures/figs2.png?raw=true)

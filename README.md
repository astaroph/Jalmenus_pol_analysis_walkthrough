# Jalmenus_pol_analysis_walkthrough
A walkthrough of the stats and data analysis I did for the Rabideau Childers and Bernard et al. 2023 paper on polarization vision in Jalmenus butterflies. 
Richard A. Rabideau Childers, Gary D. Bernard, Heqing Huang, Cheng-Chia Tsai, Mary Caswell Stoddard, Benedict G. Hogan, Joel S. F. Greenwood, Edward R. Soucy, Mark Cornwall, Matthew Lek Min Lim, Marjorie A. Liénard, Nanfang Yu, Naomi E. Pierce; A hypothesis for robust polarization vision: an example from the Australian imperial blue butterfly, Jalmenus evagoras. J Exp Biol 1 April 2023; 226 (7): jeb244515. [https://doi.org/10.1242/jeb.244515]

![Fieldwork in Ebor, NSW](images/Australia_fieldwork_compilation.png)

# Analysis of Polarization Vision in Jalmenus Evagoras

## Introduction
This repository contains a guided walkthrough of the statistical analysis I conducted in our study on polarization vision in the Australian imperial blue butterfly, *Jalmenus evagoras* and the results we found. Our research investigates how these butterflies utilize polarization vision in mate selection and navigation.

## Datasets and Scripts
- **Datasets**: Available through Dryad [https://doi.org/10.5061/dryad.kprr4xh6t]
- **Analysis Scripts**: Published on Zenodo [https://doi.org/10.5281/zenodo.7139636]

## Statistical Analysis Walkthrough
The following is a step-by-step guide through the statistical analyses performed in our study, accompanied by visualizations of our findings.


### Analysis 1: Polarization-Discrimination Behavioral Experiment
- **Objective**: In this analysis I sought to determine if *Jalmenus evagoras* can discriminate between visual stimuli based on polarization content, focusing on how they respond to colored wing models that varied in polarization content in a mate-attraction context.
- **Methods**:
In my first real field experience in grad school, my PI, Naomi Pierce, and I travelled to Ebor, New South Wales, Australia, where *Jalmenus evagoras* is most common. Naomi had been working with Gary Bernard, who studied the vision of this species, for several years already, and his data seemed to suggest that this species should be able to make use of polarized light cues, based on their unique eye physiology (more on that below). Our task was to determine whether we could get them to respond differently to light that varied only in polarization content, and whether these preferences varied with color. As this species (and many lycaenids) don't behave well under lab conditions, we journeyed to the field. We had never tried any sort of rigorous behavioral experiment with free-flying individuals in the field before, though, and we had limited funding for extended field expeditions, so I had spent the preceding two weeks working with an engineer to come up with half a dozen different kinds of ways to display polarized and unpolarized butterfly wing models hoping that at least one them would get them to respond.

Fieldwork in Ebor was pretty relaxed: drive-by fieldwork I like to call it because you can cruise around in a car and see the host plants where the caterpillars and their attending ant colonies are pretty easily from the road. We identified a good 15 different focal sites that way and got to work trying our models.Fortunately for us, our wing models on our robotic wing beating devices attracted an incredible amount of attention, and we developed an assay on the fly where we would set out our paired robotic wing beaters, one polarized, one depolarized, and count the number of times individuals clearly approached one or the other. We also kept track of the overall level of activity and whether it was sunny (butterflies out) or shaded (no activity). Unfortunately for us, about as soon as we had gotten this assay developed, we were rained out of the second of our brief two weeks in the field and had to return home. Frustrated at being tantalizingly close to getting actual results from a behavioral experiment with Lycaenid butterflies in the field (no mean feat, given how intractable this whole family of butterflies often is!), we actually booked flights for Ebor only one week later, and managed to catch the last two weeks of the 1-2 month field season there, collecting data for a solid 10-12 sites for each of three mock wing colors (red, blue and yellow). This despite spraining my ankle halfway through and having to postpone our return flights after we both spent days of fever illness with Norovirus.

![Mock wing construction details](images/Figure_2_flapper_construction_details.png)

I considered using each response as the unit of replication, but because we had no way of knowing which of the dozen or so individuals present at each site would respond over a standard 20 minute choice experiment, and thus whether these were truly 'independent' events from a statistical standpoint, I treated each site as the unit of replication, summing the total interactions with polarized or depolarized at a site and then taking the difference between them, allowing for a simple one-sample two-sided t-test to determine whether there were significantly more interactions with polarized or depolarized wing models.
- **Code**: Australia_pol_trials_analysis.R
- **Results**:
We found that they did to our initial surprise, the butterflies exhibited subtle, but consistent and significant preferences for the *depolarized* blue wing models (Student’s t-test, t=−6.2185, **P=0.0008 for polarized−depolarized=0 at N=7 sites, 183 responses total, Figure 4), with more approaches to depolarized blue wing models than polarized at every single site we used, despite strong variability in their response to polarization for other colors at the same sites (so likely not just a function of site differences). However our subsequent analyses (below) later indicated that the depolarized wing models were more spectrally and polarizationally similar to female butterflies than male ones, which likely underpins at least some of this preference. In any case, we concluded that these butterflies could indeed unambiguously discriminate polarization cues, and that this did indeed seem to vary across different colors.

- **Visualization**:
Paper Figure 4:

![Behavioral results](images/Figure_4_behavioral_field_results.png)

### Analysis 2: Reflectance Spectrophotometry and Polarization Imaging

- **Objective**: To quantify the spectral and polarization differences in wing reflectance between male and female Jalmenus evagoras, and how these differences might be perceived by their visual system.

- **Methods**:
Having shown that they can perceive these polarization differences and that they exhibit these preferences in a mating context, we next wanted to characterize the spectral and polarizational properties of males and females themselves in order to determine whether males and females exhibited significant differences from each other, especially to their eyes, along with a similar analysis of our wing models to determine how closely we had been able to match the optical properties of this species and to see if they could help explain our earlier unexpected behavioral results. Our physics collaborators, Nanfang Yu and his team at Columbia University, profiled the spectrum from the UV (365nm) through the red (700nm) of the dominant "blue-green" wing patches in the dorsal side of a number of male and female pinned butterflies, and for polarization conducted a detailed set of Degree of Polarization (DOP) measurements at a range of viewing angles, along with a similar range of illumination angles from the UV or blue light sources used.

Male-female spectral/color differences:To determine how this species might perceive these reflectance properties, Cassie Stoddard and her Postdoc, Ben Hogan, used a tetrahedral color space analysis (R package pavo), informed by the light sensitivity spectra of their visual proteins that Gary had already determined, providing me with a range of Cartesian coordinates of the "stimulation" values of these measured spectra in tetracolor space. I then used the RRPP package in R, which allows for multivariate response variables, such as my cartesian coordinate stimulation values or raw reflectance spectra, in a traditional linear modeling context to test for the effect of Sex in specimen color and spectral values in order to determine whether males and females were significantly different from each other. We also quantified how similar our wing models were to male and female specimens, bringing their reflectance values into tetracolor space and computing the pairwise euclidean distances between each model's tetracolor values and each male or female specimen.

Male-female polarizational differences: Focusing on just the blue and UV wavelength bands where Gary's work on eye physiology had found the highest polarization sensitivity, I sought to make sense of the chaos of the DOP values taken at a range of different viewing and illumination values. I knew that polarization of butterfly wings is often is highly angle dependent, so I first sought to use viewing and illumination angles, as well as Sex, as fixed effects in a mixed effects linear model with degree of polarization in blue or UV light as the response and unique ID for each individual as a random effect to account for repeated measurements ('lmer' from the lme4 package in R). However, with only a dozen or so individual males and females analyzed, I also knew that three fixed effects and their interactions would likely overfit the model. However, in my graphical analysis, I had noticed that the illumination and viewing angles were not indpendendent of each other, but rather that DOP seemed to vary most strongly with the total difference in angle between the viewing and illumination angles, which I called the "Total Angular Contrast." Using this combined factor as my fixed effect drastically simplified my models while also providing a much better fit (via AIC and log likelihood as well as via predicted and adjusted R2 values provided by the MuMIn and semEff packages).
As the relationship between DOP and total angular contrast seemed very non linear, separately for UV and blue light responses, I first made a maximal model with Sex, and a third order polynomial tramsformation of Total Angular contrast, and their interactions. I then proceded via step-wise factor reduction, eliminating lower and lower order factors in sequence and determining the effect on model fit criteria. To test for the significance of individual effects (such as sex), was more complicated: unlike traditional linear models, the inclusion of random effects covariance structures in linear models complicates the traditional assumptions (specifically the null distributions) underlying significance testing of fixed effects. However, the application of brute computational force allowed me to sidestep this issue, employing large scale parametric bootstrapping to derive an empirical null distribution from a model that contained all the fixed effects except my effect of interest, against which a larger model that differed only by this effect could be compared to produce a simple test for whether the inclusion of this effect resulted in a significant difference from the empirical null (as implemented in the PBmodcomp function of the pbkrtest package in R), no assumptions needed! It turns out the maximal models were also the best fitting, and thus these were used to generate model fits and confidence bands for subsequent plotting, using the effects package in R.

- **Code**: Spectral/colorspace: Spectral_reflectance_analysis_RRPP_standardonly_v2.R ; DOP analysis: DOP_angular_contrast_script_cleaned.R

- **Results**:
With the high number of variables (wavelengths) in the raw data, there was insufficient power to resolve the differences between male and female spectra. However, when we looked at these differences 'through their eyes' as it were (Figure 5E), in the tetracolorspace values, we found a highly significant sex difference for both forewings and hindwings (F=12.51, adjusted P=0.002 for forewings, F=36.28, adjusted P=0.0004 for hindwings). Polarizationally, they showed the strongest differences in both blue and UV illumination at very oblique angles (Figure 5 C and D), but only for DoP in blue were the interactions of sex and polynomial contrast significant (sex×angular contrast first-order: t=11.45, P=0.0008; sex×angular contrast second-order: t=4.44, P=0.0372). This could explain why using Total Angular Contrast as a fixed effect provided such an improvement to my model fits, as it better reflects (if you will) the physical and biological reality of how polarization is produced and observed in this species, as males and females most critically evaluate each other facing head on, where the combined solar illumination and viewing angles are quite large. Overall, we could easily conclude that subtle, but highly significant differences in coloration and polarization underpinned the clearly visible sexual dimorphism of this species. We also observed that the blue wing models were the least dissimilar in tetracolorspace to the female and then male butterfly specimens (Figure 6C), and that the depolarized blue wing models were substantially closer in polarizational space to the actual specimens than the polarized blue wing models (Figure 6E and F). Taken together, this could potentially explain why these models had the strongest interest of the three colors we tested and the only significant polarizational preference we observed.

- **Visualization**:

Paper Figure 5

![Reflectance and tetracolor space](images/Figure5_reflectanc_pol_tetracolor.png)


Paper Figure 6:

![Specimen and wing model spectra and polarization](images/Figure6_pol_spectral_differences_and_wing_models.png)


### Analysis 3: Ommatidial Alignment Measurement

- **Objective**: To introduce a novel method for assessing the alignment of microvilli within the ommatidial arrays of Jalmenus evagoras and its implications for polarization vision.
- **Methods**: Gary's method (in a nutshell), relies on using a unique property of the eyes of some butterflies, their eyeshine or light reflected by their inner eyes when illuminated, to determine the orientation of their inner ommatidial waveguides ('rhabdoms'). Our overarching premise in this paper is that *Jalmenus evagoras* is a flying organism in a visually complex habitat that would greatly benefit from a means of polarization vision that does not rely on maintaining a fixed head angle with respect to their environment, unlike other species that use polarized light, such as crabs or shrimp. Our previous analysis had essentially established that *Jalmenus* has the 'motive and opportunity' for polarization vision, but Gary's method aimed to show that they also had the physiological 'means,' as well. He hypothesized that the substantial 'misalignment' of their arrays of rhabdoms in the individual ommatidia of their compound eyes could allow them to integrate information between these ommatidia with contrasting angles, and thus derive a means of polarization vision that was robust to changes in viewing angle. In order to see whether these misalignment patterns also varied across the eyes or between males and females, we sought to characterize a number of ommatidial eye 'patches' that he had measured, grouping them into subarrays of ommatidia that were misaligned and could thus allow for this novel polarization vision (PD for polarization detectors), and those that were more traditionally aligned and could thus be more useful for traditional detection of edges and boundaries (ED for edge detectors). I sought to characterize the frequency of these PD and ED ommatidia in relation to the vertical location of the patch on the eye (more dorsal regions often are involved in detecting polarization in butterflies), and critically, between malkes and females, and once again employed a similar mixed effects  linear modeling approach followed by parametric bootstrapping, as described above, to test for the effect of sex, eye patch elevation and their interaction on the percentage of PD or ED ommatidia.

- **Code**:Gary_PD_ED_linear_model_cleaned.R

- **Results**: 
We found that PD ommatidia did indeed increase significantly with increasing eye patch elevation (t=19.885, P=0.0001; Fig 7E below), as expected from the literature. We also found that females (yellow bars and points) tended to have more PD overall (t=7.977, P=0.012) , and especially at the lower eye elevations with which they would most likely observe male counterparts in mating interactions (Fig.7E below). 
- **Visualization**:

Paper Figure 7:

![Ommatidial Alignment](images/Figure7_EDPD_ommatidia_testing.png)


For further details on the methodology, datasets, and analyses, please refer to our published paper and the resources available on Dryad and Zenodo.


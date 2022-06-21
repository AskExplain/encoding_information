# Encoding large information structures in linear algebra

An experiment is shown with a linear mixed model. In genetics, a Genetic Relatedness Matrix is used as the covariance structure for the random effects. With encoding of the sample information, the sample structure of this covariance matrix can be reduced from (N samples) to (M samples) where (N > M).

This is for 100 permutation runs of a simulation with 1000 samples and 100 SNPs simulated according to:


![alt text](https://raw.githubusercontent.com/AskExplain/encoding_information/alpha_test_v2022.1/figures/encoded_vs_original_mixed_model.png)


Given information is being encoded, it is expected for there to be information loss leading to higher variability compared to the standard linear mixed model. However, due to reduced sample size via an encoding the runtime is faster - almost half the speed of the original mixed model according to the R package GMMAT. HE regression which is known to be much faster also has a greater bias in the heritability estimate.

Future experiments will involve mixture models with encoded features.

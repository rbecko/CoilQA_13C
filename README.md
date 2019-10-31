# CoilQA_13C

QA acquisition and reconstruction protocol for 13C RF receive coil SNR evaluation, to facilitate regular coil tests and comparisons locally, across sites, and in literature. Special care was taken to assure an unbiased SNR comparison both across vendors by accounting for scanner-specific filters, and between coils with different numbers of receive elements by performing sensitivity based coil combination.

Besides SNR estimation the source code also provides estimates of noise correlations, G-factors, and T2* values. The QA protocol is based on imaging with a non-enriched 13C ethylene glycol phantom.

The QA protocol follows the SNR estimation principles described by Kellman & Mcveigh.

The flowchart below provide an overview of the QA protocol. “s” refers to the extracted signal estimate, “σ^2” to the scaled noise variance, “Ψ” to the scaled noise covariance matrix, and “b” to the estimated coil sensitivity profiles. The factor of √2  in the final SNR estimation equations is included to account for the SNR definition in terms of the real channel noise component.

![alt text](https://github.com/rbecko/CoilQA_13C/blob/master/Fig1_ISMRM2020.png)

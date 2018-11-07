SeisFlows_SRVM
=================

=================update@Nov 7st,2018

3D local-scale SRVM numerical examples are oversized, but available upon request.

=================update@Nov 1st,2018 

with checkpointing function for Seisflows_SRVM

with simplest and robust fault tolerance solution

When seisflows_3D runs on supercomputer, the faults usually occur in the adjoint modelings in which massive I/Os are greatly needed. The forward modeling parts are basically safe throughout the inversion.

also, adjoint modeling is twice expensive as that of forward modeling, so we specify adaptive wall-time for them to make them more competitive in queuing on HPCs.

suggest using I/O optimisation softwares when run it on HPCs

=================update@Oct 31, 2018, 3D local-scale SRMV testings have been basically finished.

SeisFlows_SRVM is a special version of SeisFlows.

SeisFlows_SRVM contains SRVM, an additional optimization algorithm, in addition to SD, CG, and L-BFGS in the original Seisflows.

The full name of SRVM is square-root variable metric, which originates from DFP, the dual of BFGS.

SRVM based FWI allows for a posterior uncertainty estimation after the acoustic/elastic inversion.

SRVM works in a matrix-free vector version in elastic FWI, making it memory-affordable.

Documentation about the algorithm in Seisflows_SRVM is available online at readthedocs.org, except for an additional option: SRVM. 

One manuscript about how SRVM works during FWI is in preparation and will be submitted to GJI soon.

One EGU poster and one EAGE abstract can be found in this folder.

The codes for uncertainty estimation are written in Matlab. The posterior analysis part includes RSVD-SRVM, standard deviation map, eigenvectors and eigenvalues of the inverse Hessian, 2D Gaussian random samplers, 2D posterior sampling, marginal distributions, and null space. Please find the Matlab codes in this GitHub repository: Uncertainty_analysis_example-Marmousi 

Several SRMV-based manuscripts about the theories and methods of uncertainty estimation regarding acoustic/elastic FWI are in preparation, and will be submitted to GJIs and so on soon.

One manuscript about SRVM-based null-space shuttle in elastic FWI is in preparation.

SeisFlows is an open source seismic inversion package that

    Delivers a complete, customizable waveform inversion workflow

    Provides a framework for research in regional, global, and exploration seismology

Documentation is available online at readthedocs.org. A manuscript is currently under review in Computers and Geosciences.

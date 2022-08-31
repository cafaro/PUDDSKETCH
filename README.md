# PUDDSKETCH


A C++ implementation of a parallel version of the UDDSketch algorithm [\[1\]][1] [\[2\]][2] for quantile estimation.

We also provide:

- a parallel version of DDSketch (see the ParallelDDSketch directory);
- parallel version of both the KLL and REQ quantile algorithms, based on the Apache DataSketches library.

# Build
Tested on Linux and on MacOS; it shoukld work on other UNIX-like OSes.

If you use this software, please cite the following paper:

M. Cafaro, C. Melle, I. Epicoco, M. Pulimeno. Data stream fusion for accurate quantile tracking and analysis. Information Fusion, Elsevier, Volume 89, 2023, Pages 155-165, ISSN 1566-2535, DOI: 10.1016/j.inffus.2022.08.005

# References

\[1\] **I. Epicoco, C. Melle, M. Cafaro, M. Pulimeno and G. Morleo**. UDDSketch: Accurate Tracking of Quantiles in Data Streams. IEEE Access, vol. 8, pp. 147604-147617, 2020, ISSN: 2169-3536, DOI: 10.1109/ACCESS.2020.3015599

[2]: <https://github.com/cafaro/UDDSketch> "UDDSketch"

[1]: <https://www.sciencedirect.com/science/article/pii/S1566253522000975> "Data stream fusion for accurate quantile tracking and analysis"

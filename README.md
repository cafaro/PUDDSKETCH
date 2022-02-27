# PUDDSKETCH


A C++ implementation of a parallel version of the UDDSketch algorithm [\[1\]][2] for quantile estimation.

We also provide:

- a parallel version of DDSketch (see the ParallelDDSketch directory);
- parallel version of both the KLL and REQ quantile algorithms, based on the Apache DataSketches library.


# Build
Tested on Linux and on MacOS; it shoukld work on other UNIX-like OSes.



# References

\[1\] **I. Epicoco, C. Melle, M. Cafaro, M. Pulimeno and G. Morleo**. UDDSketch: Accurate Tracking of Quantiles in Data Streams. IEEE Access, vol. 8, pp. 147604-147617, 2020, ISSN: 2169-3536, DOI: 10.1109/ACCESS.2020.3015599

[2]: <https://github.com/cafaro/UDDSketch> "UDDSketch"

# Timing results
Note: 64 q-points, 64 thread blocks. My P2000 can run 168 at the same time.

Note: My GPU is bad at double precision.

## Silicon
Version | Old | New | New 2
--- | --- | --- | ---
Serial | 2.7 s | 2.6 s | 2.4 s
OpenMP | 0.13 s | 0.06 s | 0.05 s
CUDA | 0.5 s | 0.05 s | 0.05 s
CUDA V100 | 0.8 s | - | 0.017 s

## SiC
Version | Old | New | New 2
--- | --- | --- | ---
Serial | 62 s | 82 s | 70 s
OpenMP | 17 s | 20 s | 18 s
CUDA | 13 s | 13 s | 7.5 s
CUDA V100 | 2.7 s | - | 0.37 s

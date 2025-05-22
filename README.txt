Proof-of-concept SageMath implementation associated to the paper 
"Isogeny-based key exchange from orientations of large discriminant"

The main code can be found in the folder /scheme.
The 256-bit example in the paper can be evaluated using

sage benchmark_256bits.py

In the folder dCSIDH, we have implemented an unoptimized version
of the dCSIDH protocol as described in ia.cr/2023/793
for a prime of bitsize 4096.
This corresponds to the "conservative"
NIST level 1 parameter set.
It can be benchmarked using

sage benchmark.sage

The folder external_modules contains the packages

- deuring
by: Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni
from: ia.cr/2023/106
for: computing the effective Deuring correpondence

- KummerIsogeny
by: Giacomo Pope
from: github.com/GiacomoPope/KummerIsogeny
for: x-only arithmetic on Montgomery curves



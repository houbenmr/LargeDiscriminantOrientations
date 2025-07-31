Proof-of-concept SageMath implementation associated to the paper 
"Isogeny-based key exchange from orientations of large discriminant"

https://ia.cr/2025/1098

The main code can be found in the folder code/scheme/
The 256-bit example in the paper can be evaluated using

sage bench/benchmark_256bits.py

This does 25 full key exchanges (4 ideal actions each).
The average results are saved in bench/bench-256.out

An unoptimized version of the dCSIDH protocol from ia.cr/2023/793 can be evaluated using

sage bench/benchmark_dCSIDH4096.sage

This does 20 "half" key exchanges (2 ideal actions each).

The folder code/external_modules/ contains the packages

- deuring
by: Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni
from: ia.cr/2023/106
for: computing the effective Deuring correpondence

- KummerIsogeny
by: Giacomo Pope
from: github.com/GiacomoPope/KummerIsogeny
for: x-only arithmetic on Montgomery curves

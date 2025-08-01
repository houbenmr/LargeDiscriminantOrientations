# Bruno Sterner prime for lambda = 256, Apr 14, 2025
p = 0xaee61f0e7eb0b904142a14b1ba552ca011c6104b007e448fe6aed1f1af734f1f
r = 13
lmbda = 256

from time import time
from parameters import setup
params = setup(p, r, lmbda = lmbda)

import scheme as scheme
from scheme_cost import exps_to_B
scheme.params = params
B = exps_to_B(params.exps, lmbda)
print('using B =', B)

#benchmark

timings = []

for i in range(25):
    sk_alice, sk_bob = scheme.keygen(B), scheme.keygen(B)
    start = time()
    pk_alice = scheme.group_action(scheme.params.base_cycle, sk_alice, B)
    end = time()
    print('Alice public key generation took time', end-start)
    timings.append(end-start)
    start = time()
    pk_bob = scheme.group_action(scheme.params.base_cycle, sk_bob, B)
    end = time()
    print('Bob public key generation took time', end-start)
    timings.append(end-start)
    start = time()
    ss_alice = scheme.group_action(pk_bob, sk_alice, B)
    end = time()
    print('Alice shared secret generation took time', end-start)
    timings.append(end-start)
    start = time()
    ss_bob = scheme.group_action(pk_alice, sk_bob, B)
    end = time()
    print('Bob shared secret generation took time', end-start)
    timings.append(end-start)
    assert ss_alice == ss_bob

from numpy import median, std

with open('bench.out', 'w') as op:
    op.write('Timing results for a 256-bit group action' + '\n')
    op.write('The median is ' + str(median(timings)) + '\n')
    op.write('The standard deviation is ' + str(std(timings)))

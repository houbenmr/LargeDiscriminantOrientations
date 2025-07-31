# "more optimal" CSIDH-514.6 prime, May 12, 2025
# trace = 1569142699, 2119947157, 2322963581 for r = 8
# given Ms = (p+1) // 4 and Mt = 1

p = 0x62a78c7d5d89de462a94f1d2be9573eacd717abfca5e5e9993dc08861e702f7dc595cf96c3f2c8c488bc02e165a095ea79ebc348b269b96d12585a23966195b43
r = 8
lmbda = 256

import sys
import os
from pathlib import Path
current_file = Path(__file__).resolve()
current_folder = current_file.parent
sys.path.append(str(current_folder.parent / 'code/scheme'))

from time import time
from parameters import setup
params = setup(p, r, lmbda = lmbda, Ms = (p+1) // 4, Mt = 1)

import scheme_without_twist as scheme
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

with open(current_folder / 'bench-512.out', 'w') as op:
    op.write('Timing results for a 512-bit group action' + '\n')
    op.write('The median is ' + str(median(timings)) + '\n')
    op.write('The standard deviation is ' + str(std(timings)))

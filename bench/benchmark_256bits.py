# Bruno Sterner prime for lambda = 256, Apr 14, 2025
p = 0x53d1c0debdc0c1ba3edd760148d21a078008f642a9f26f65995c4d86fc6efff7
r = 13
lmbda = 256

import sys
import os
from pathlib import Path
current_file = Path(__file__).resolve()
current_folder = current_file.parent
sys.path.append(str(current_folder.parent / 'code/scheme'))

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

with open(current_folder / 'bench-256.out', 'w') as op:
    op.write('Timing results for a 256-bit group action' + '\n')
    op.write('The median is ' + str(median(timings)) + '\n')
    op.write('The standard deviation is ' + str(std(timings)))

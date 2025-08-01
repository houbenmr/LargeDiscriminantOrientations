# "more optimal" CSIDH-514.6 prime, May 12, 2025
# trace = 1569142699, 2119947157, 2322963581 for r = 8
# given Ms = (p+1) // 4 and Mt = 1

p = 0x62a78c7d5d89de462a94f1d2be9573eacd717abfca5e5e9993dc08861e702f7dc595cf96c3f2c8c488bc02e165a095ea79ebc348b269b96d12585a23966195b43

r = 8

lmbda = 256

from time import time
from parameters import setup

params = setup(p, r, generate = True, write = True, trace = 1569142699, lmbda = lmbda)

import scheme as scheme
from scheme_cost import exps_to_B

scheme.params = params

B = exps_to_B(params.exps, lmbda)

print('using B =', B)

sk_alice, sk_bob = scheme.keygen(B), scheme.keygen(B)

print(sk_alice, sk_bob)

start = time()
pk_alice = scheme.group_action(scheme.params.base_cycle, sk_alice, B)
end = time()
print('Alice public key generation took time', end-start)
start = time()
pk_bob = scheme.group_action(scheme.params.base_cycle, sk_bob, B)
end = time()
print('Bob public key generation took time', end-start)
start = time()
ss_alice = scheme.group_action(pk_bob, sk_alice, B)
end = time()
print('Alice shared secret generation took time', end-start)
start = time()
ss_bob = scheme.group_action(pk_alice, sk_bob, B)
end = time()
print('Bob shared secret generation took time', end-start)

assert ss_alice == ss_bob
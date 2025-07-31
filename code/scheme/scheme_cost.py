from sage.all import *
from math import log, prod

def ell_isogeny_cost(x):
    if x < 250: # Using Velu
        return 0.106 * x**1.028 + 1.705
    if x > 250: # Using sqrtVelu
        return x**0.650

def exps_to_B(exps, lmbda):
    
    # Input:
    # - exponents of prime factors of p^2-1
    # - security parameter lambda
    # Output:
    # - required number of iterations B to achieve a key space of size lambda

    B = 0
    
    while True:
        B += 1
        key_space_size = sum([log(B*e + 1., 2) for e in exps])
        if key_space_size > lmbda:
            break

    return B

def scheme_cost(ells, exps, lmbda):

    # Input:
    # - sorted list of prime factors of p^2-1
    # - their exponents
    # - security parameter lambda
    # Output:
    # - the cost of the scheme achieving security lambda

    assert len(ells) == len(exps)
    n = len(ells)

    B = exps_to_B(exps, lmbda)

    cost = B * sum([ell_isogeny_cost(ells[i])*exps[i] for i in range(n)])

    return cost

def minimize_scheme_cost(ells, exps, lmbda):

    # greedily remove exponents that decrease the cost the most

    assert len(ells) == len(exps)
    n = len(ells)
    
    while True:
        
        proceed = False

        best = scheme_cost(ells, exps, lmbda)
    
        for i in range(n):
            exps_try = [max(0, exps[j] - int(i == j)) for j in range(n)]
            cost_try = scheme_cost(ells, exps_try, lmbda)
            if cost_try < best:
                best = cost_try
                set_to_remove = i
                proceed = True

        if proceed:
            exps[set_to_remove] = exps[set_to_remove]-1
        else:
            break

    return best, exps, exps_to_B(exps, lmbda)

def prime_score(p, lmbda):

    p = ZZ(p)

    assert p%4 == 3
    
    tups = (p**2-1).odd_part().factor()
    exp_two = (p+1).valuation(2)
    ells = [2] + [ell for [ell,r] in tups]
    exps = [exp_two] + [r for [ell,r] in tups]

    score,_,_ = minimize_scheme_cost(ells, exps, lmbda)

    bitsize = log(p,2)

    return bitsize, score

def optimal_configuration(ells, bitlimit, lmbda):

    # Input:
    # - a list of ells
    # - a bit limit for prod(ell^exp)
    # - a security target lambda
    # Output:
    # - an exponent vector exps minimizing scheme_cost(ells, exps, lmbda)
    # Remark: this does not ensure any primality condition

    n = len(ells)
    exps = [1] * n
    
    m = 1

    while m:
        
        best = scheme_cost(ells, exps, lmbda)
        
        for b in reversed(range(n)):
            if exps[b] == m:
                break
                
        for a in reversed(range(n)):
            if exps[a] == m+1:
                a += 1
                break

        ell = ells[b]
        exps_temp = [exps[i] - int(i==b) for i in range(n)]
    
        for j in range(a,n):
            exps_temp[j] += 1
            M = prod([ells[k]**exps_temp[k] for k in range(n)])
            bitsize = log(M, 2)
            if bitsize > bitlimit:
                break
    
        new_score = scheme_cost(ells, exps_temp, lmbda)
        
        if new_score < best:
            best = new_score
            m += 1
        else:
            m -= 1

    return exps

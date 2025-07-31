import time

def cornacchiax(d,m): # finds the x-coordinate of a primitive solution to x^2 + dy^2 = m
    if not m > 0:
        return 0
    if not mod(-d,m).is_square():
        return 0
    r = int(mod(-d,m).sqrt())
    a = m
    b = r
    while  a^2 >= m:
        c = a%b
        a = b
        b = c
    s2 = (m-a^2)/d # if s2 = s^2 is a square, then a^2 + ds^2 = m is a primitive solution
    if s2.is_square():
        return a
    return 0

# Bruno's 256-bit twin smooth prime

l = 38295031104
p = 2*(l*(l+6)*(l+13)*(l+19) // 1080)^2-1

exp_two = (p+1).valuation(2)
power_of_two = 2^(exp_two)

tups = (p^2-1).odd_part().factor()

ells = [ell for [ell,r] in tups]
exps = [r for [ell,r] in tups]

M = prod([ell^r for [ell,r] in tups])

number_of_odd_primes = len(ells)
Fp = GF(p)

print('the value of p is',p)
print('the bitsize of p is',n(log(p,2)))
print('the bitsize of M is',n(log(M,2)))
print('the number of odd ell_i is',number_of_odd_primes)
print('the largest ell_i is',max(ells))

r = 11

target_norm = M^r

x = max(ells)

while True:
    while True:
        x = x.next_prime()
        factor_1 = 4*M^r-x^2
        if factor_1.is_pseudoprime():
            break
    
    y_squared_mod_p = Fp(4*target_norm-x^2)
    if y_squared_mod_p.is_square():
        y_mod_p = Fp(y_squared_mod_p).sqrt()
        y = ZZ(y_mod_p)
        if is_odd(y):
            y = p-y
        m = (4*target_norm - (x^2 + y^2)) // p
        if m%4 == 1 and m.is_pseudoprime():
            print('found candidate m', m)
            start = time.time()
            c = cornacchiax(1,m)
            end = time.time()
            print('cornacchias algorithm took time', end-start)
            if not c == 0:
                d = ZZ((m-c^2).sqrt())
                if is_odd(c):
                    c,d = d,c
                break
                
a = ZZ((x-d)/2)
b = ZZ((y-c)/2)

sigma = [a,b,c,d]
sigma_dual = [a+d,-b,-c,-d]

print('Orientation found.')

print('sigma =', sigma)

print('the trace of sigma is',x)

Disc_sigma = 4*M^r-x^2

print('Verifying that the discriminant is prime...')

start = time.time()

assert Disc_sigma.is_prime()

end = time.time()

print('Done. This took time', end - start)

print('The bitsize of the discriminant of sigma is',n(log(Disc_sigma,2)))

R.<X> = Fp[]
Fq.<i> = Fp.extension(X^2+1)
E_straight = EllipticCurve([Fq(1),0])
E_twist = E_straight.quadratic_twist()

ells_straight = [ell for [ell,r] in tups if ell.divides(p+1)]
ells_twist = [ell for [ell,r] in tups if ell.divides(p-1)]

ell_powers = [ell^r for [ell,r] in tups]
ell_powers_straight = [ell^r for [ell,r] in tups if ell.divides(p+1)]
ell_powers_twist = [ell^r for [ell,r] in tups if ell.divides(p-1)]

exps_straight = [r for [ell,r] in tups if ell.divides(p+1)]
exps_twist = [r for [ell,r] in tups if ell.divides(p-1)]

M_straight = prod(ell_powers_straight)
M_twist = prod(ell_powers_twist)

assert M == M_straight * M_twist

# Functions to evaluate endomorphisms from Z[1,i,pi,i*pi] on elements of E_straight and E_twist.

def aut_i(P):
    
    E = P.curve()
    a4 = E.a4()
    assert E == EllipticCurve([a4,0])
    
    if P == E(0):
        return P
    else:
        xP = P[0]
        yP = P[1]
        xQ = -xP
        yQ = i*yP
        return E(xQ,yQ)

def frob_p(P):
    
    E = P.curve()
    a4 = E.a4()
    assert E == EllipticCurve([a4,0])
    
    if P == E(0):
        return P
    else:
        # correction factor for a possible quadratic twist
        alpha = a4.sqrt() 
        alpha = sorted([alpha,-alpha])[0] # take a deterministic sqrt
        v = alpha^((1-p) // 2)
        
        xP = P[0]
        yP = P[1]
        xQ = v^2 * xP^p
        yQ = v^3 * yP^p
        return E(xQ,yQ)

def endo_basic(vec,P):
    
    # Input: vec = [a,b,c,d] in Z^4, P in E.
    # Output: the evaluation of a+b*i+c*pi+d*i*pi at P.
    
    [a,b,c,d] = vec
    iP = aut_i(P)
    jP = frob_p(P)
    kP = aut_i(jP)
    Q = a*P + b*iP + c*jP + d*kP
    return Q

# Evaluate endomorphisms from the maximal order Z[1,i,(i+pi)/2,(1+i*pi)/2].

def endo_max(vec,P):
    
    # Input: vec = [a,b,c,d] in Z^4, and P in E.
    # Output: the evaluation of a+b*i+c*(i+pi)/2+d*(1+i*pi)/2 in End(E) at P.
    
    E = P.curve()
    [a,b,c,d] = vec
    N = order(P)
    if is_odd(N): # In this case, dividing by two is easy.
        inverse_of_2_mod_N = mod(2,N)^(-1)
        P_divided_by_two = ZZ(inverse_of_2_mod_N) * P
    else:
        divide_by_two_field.<X> = P.division_points(2,poly_only=True).splitting_field()
        Eq2 = E.base_extend(divide_by_two_field)
        P_divided_by_two = (Eq2(P).division_points(2))[0]
    Q = endo_basic([2*a+d,2*b+c,c,d], P_divided_by_two)
    return E(Q)

# Find a basis P,Q of E[ell] such that P in ker(sigma) and Q in ker(sigma_dual)

def find_kernel_basis_max(E, sigma, ell):

    N = ZZ(sqrt(E.cardinality()))

    [a,b,c,d] = sigma
    
    sigma_dual = [a+d,-b,-c,-d]

    while True:
        R = E.random_point()
        cofactor = N // ell
        R *= cofactor
        P = endo_max(sigma_dual, R)
        if order(P) == ell:
            break
    
    while True:
        R = E.random_point()
        cofactor = N // ell
        R *= cofactor
        Q = endo_max(sigma, R)
        if order(Q) == ell:
            break

    return P,Q

def full_order_kernel_points(E, sigma, ell_powers):

    # We construct a basis P,Q for E[M] where M = prod ell_powers
    # such that P in ker(sigma) and Q in ker(sigma_dual).
    
    P = E(0)
    Q = E(0)
    
    for ell_power in ell_powers:
        P0,Q0 = find_kernel_basis_max(E, sigma, ell_power)
        P += P0
        Q += Q0

    assert P.order() == prod(ell_powers)
    assert Q.order() == prod(ell_powers)
    assert endo_max(sigma, P) == E(0)
    assert endo_max(sigma_dual, Q) == E(0)
    
    return P,Q

P0_straight, Q0_straight = full_order_kernel_points(E_straight, sigma, ell_powers_straight)

P0_twist, Q0_twist = full_order_kernel_points(E_twist, sigma, ell_powers_twist)

# Deuring

import sys

sys.path.append('deuring')

from deuring.broker import starting_curve
from deuring.randomideal import random_ideal
from deuring.correspondence import constructive_deuring, multiplyIdeals
from xonly import xPoint

E0, iota, O0 = starting_curve(Fq)

# Evaluate elements of O0 on points of E0

def endo_max_quat(alpha, P):
    [a,b,c,d] = alpha.coefficient_tuple()
    [a,b,c,d] = [a-d,b-c,2*c,2*d]
    return endo_max([a,b,c,d], P)

the_real_sigma = sum([sigma[k]*O0.basis()[k] for k in range(4)])
the_real_sigma_dual = the_real_sigma.conjugate()

positive_ideal_list = [O0.left_ideal([M^e, the_real_sigma]) for e in range(r)]
negative_ideal_list = [O0.left_ideal([M^e, the_real_sigma_dual]) for e in range(r)]

def push_orientation(E_base, negative_ideal, P0, M):

    # Input:
    # - a curve E_base, together with a negative_ideal connecting E0 to E_base
    # - a point P0 of order M on (possibly a quadratic twist of) E0
    # such that:
    # (*) the subgroups corresponding to negative_ideal and P0 intersect trivially.
    # Output: the x-coordinate of the image of P0 under the isogeny corresponding to negative_ideal
    # Note: By (*), the corresponding point still has order M

    E_home = P0.curve() # The curve where P0 lives

    assert E_home == EllipticCurve([E_home.a4(),0]) # Assert that it is isomorphic to E0

    banned_integer = ZZ(radical(M))
    
    print('e =', e, 'banned =', banned_integer.factor())
    start = time.time()
    E_temp, psi, beta, J, _ = constructive_deuring(negative_ideal, E0, iota, banned_integer)
    end = time.time()
    print('constructive deuring took time', end-start)
    
    assert E_temp.is_isomorphic(E_base)
    assert len(E_temp.isomorphisms(E_base)) == 2
    iso = E_temp.isomorphisms(E_base)[0]
    endo_ideal = multiplyIdeals(negative_ideal, J.conjugate(), beta=beta)
    assert endo_ideal.is_principal()
    assert endo_ideal == O0.left_ideal(beta)
    
    beta_P0 = endo_max_quat(beta, P0) # This lives on E_home, i.e. possibly a quadratic twist of E0
    
    a4 = E_home.a4()
    alpha = a4.sqrt() # the twisting parameter
    alpha = sorted([alpha,-alpha])[0] # replace by a deterministic sqrt

    xP = alpha^(-1) * beta_P0.x() # The x-coordinate of beta(P0) as a point on E0
    xP = xPoint(xP, E0) # Using the xPoint class to prepare evaluation under psi
    xP = xP.push(psi) # Now we are on E_temp
    xP = iso.x_rational_map()(xP.X) # Translate the point to E_base
    
    return xP

# Warning: no 2-isogenies allowed past this point!

# Using Kummer arithmetic

import sys

sys.path.append('KummerIsogeny')

from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny

Ms = M_straight
Mt = M_twist

a4 = E_twist.a4()
alpha = a4.sqrt() # the twisting parameter
alpha = sorted([alpha,-alpha])[0] # replace by a deterministic sqrt

def multiple_point_KummerLineIsogeny(K, kernel_list, degree_list, eval_list = []):
    while kernel_list:
        kernel_point = kernel_list.pop(0)
        degree = degree_list.pop(0)
        phi = KummerLineIsogeny(K, kernel_point, degree)
        K = phi.codomain()
        kernel_list = [phi(P) for P in kernel_list]
        eval_list = [phi(P) for P in eval_list]
    if eval_list:
        return K, eval_list
    else:
        return K

# Generate public parameters

data = []

K_home = KummerLine(E0)
Ps, Pt = K_home(P0_straight[0]), K_home(alpha^(-1) * P0_twist[0])
Qs, Qt = K_home(Q0_straight[0]), K_home(alpha^(-1) * Q0_twist[0])

assert Mt * Pt == K_home.zero()
assert Mt * Qt == K_home.zero()

A = K_home.a()
xPs, xPt = Ps.x(), Pt.x()
xQs, xQt = Qs.x(), Qt.x()

data.append([A, xPs, xPt, xQs, xQt])

for e in range(1,r):

    K_base, [Qs, Qt] = multiple_point_KummerLineIsogeny(K_home, [Ps, Pt], [Ms, Mt], [Qs, Qt])

    E_base = K_base.curve()

    A = K_base.a()
    
    xQs = Qs.x()
    xQt = Qt.x()
    
    negative_ideal = negative_ideal_list[-e]
    
    xPs = push_orientation(E_base, negative_ideal, P0_straight, M_straight)
    xPt = push_orientation(E_base, negative_ideal, P0_twist, M_twist)

    data.append([A, xPs, xPt, xQs, xQt])

    K_home = K_base

    Ps, Pt = K_home(xPs), K_home(xPt)

    assert Ms * Ps == K_home.zero()
    assert Mt * Pt == K_home.zero()

with open("order_11.txt", 'w') as op:
    op.write(str(data))

# test

start = time.time()
for e in range(r):
    A, Ps, Pt, Qs, Qt = data[e]
    K = KummerLine(Fq, [A,1])
    assert multiple_point_KummerLineIsogeny(K, [K(Ps), K(Pt)], [Ms, Mt]).j_invariant() == KummerLine(Fq, [data[(e+1)%r][0], 1]).j_invariant()
    assert multiple_point_KummerLineIsogeny(K, [K(Qs), K(Qt)], [Ms, Mt]).j_invariant() == KummerLine(Fq, [data[(e-1)%r][0], 1]).j_invariant()
end = time.time()

print('Estimate of the cost of two group actions:', end-start)

import time
import sys
import csv
import ast  # To safely evaluate string representations of lists

from sage.all import *

from pathlib import Path
current_file = Path(__file__).resolve()
current_folder = current_file.parent
sys.path.append(str(current_folder.parent / 'external_modules'))
sys.path.append(str(current_folder.parent / 'external_modules/KummerIsogeny'))

from deuring.broker import starting_curve
from deuring.randomideal import random_ideal
from deuring.correspondence import constructive_deuring, multiplyIdeals
from deuring.xonly import xPoint

from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny

from scheme_cost import minimize_scheme_cost

from math import prod, log

def cornacchiax(d, m): # finds the x-coordinate of a primitive solution to x^2 + dy^2 = m
    if not m > 0:
        return 0
    if not mod(-d,m).is_square():
        return 0
    r = int(mod(-d,m).sqrt())
    a = m
    b = r
    while  a**2 >= m:
        c = a%b
        a = b
        b = c
    s2 = (m-a**2)/d # if s2 = s^2 is a square, then a^2 + ds^2 = m is a primitive solution
    if s2.is_square():
        return a
    return 0

def find_sigma(p, r, M, trace = 0, num_cores = 1):
    
    if trace:
        print('trace prespecified; looking for an endomorphism with trace', trace)
        x = ZZ(trace - 1)
    else:
        print('trace not prespecified; computing orientation...')
        x = ZZ(1)
    
    target_norm = M**r

    if target_norm < p:
        raise ValueError('The desired norm is too small')
    
    Fp = GF(p)
    
    while True:
        while True:

            for _ in range(num_cores):
                x = x.next_prime() # silly method to parallelize

            candidate_discriminant = 4*M**r-x**2
            if r%2:
                if candidate_discriminant.is_pseudoprime():
                    break
            else:
                factor_1 = 2*M**(r//2)+x
                factor_2 = 2*M**(r//2)-x
                if factor_1.is_pseudoprime() and factor_2.is_pseudoprime():
                    break
        
        y_squared_mod_p = Fp(4*target_norm - x**2)
        if y_squared_mod_p.is_square():
            y_mod_p = Fp(y_squared_mod_p).sqrt()
            y = ZZ(y_mod_p)
            if is_odd(y):
                y = p-y
            m = (4*target_norm - (x**2 + y**2)) // p
            if m%4 == 1 and m.is_pseudoprime():
                c = cornacchiax(1,m)
                if not c == 0:
                    d = ZZ((m-c**2).sqrt())
                    if is_odd(c):
                        c,d = d,c
                    break
                    
    a = ZZ((x-d)/2)
    b = ZZ((y-c)/2)
    
    sigma = [a,b,c,d]
    sigma_dual = [a+d,-b,-c,-d]
    
    print('Found endomorphism of trace', x)
    
    if num_cores > 1:
        with open(current_folder / ('found_trace_started_from' + str(trace) + '.txt'), 'w') as op:
            op.write(str(x))
    
    Disc_sigma = 4*M**r-x**2

    start = time.time()
    if r%2:
#       Verifying that the discriminant is really prime
        if not Disc_sigma.is_prime():
            raise ValueError('the discriminant is pseudoprime but not prime', Disc_sigma) # This shouldn't happen
    else:
        factor_1 = 2*M**(r//2)+x
        factor_2 = 2*M**(r//2)-x
#       Verifying that the two factors of the discriminant are really prime
        if not factor_1.is_prime() or not factor_2.is_prime():
            raise ValueError('one of the factors of the discriminant is pseudoprime but not prime', Disc_sigma) # This shouldn't happen
    end = time.time()
    
#    print('Verification complete. This took time', end - start)
    print('The bitsize of the discriminant of sigma is', log(int(Disc_sigma), 2))

    return sigma, x
    
class SpecialExtremalCurve:
    def __init__(self, E, i):

        # i, a square root of -1 in the base field
        # E, a (possibly trivial) quadratic twist of y^2 = x^3 + x
        
        F = E.base_field()
        assert i**2 == F(-1)

        p = F.characteristic()
        assert p%4 == 3
        
        a4 = E.a4()
        assert E == EllipticCurve([a4, 0])
        
        self.curve = E
        self.base_field = F
        self.p = p
        self.i = i
        self.a4 = a4

        # parameters correcting the quadratic twist
        alpha = a4.sqrt()
        alpha = sorted([alpha,-alpha])[0] # take a deterministic sqrt
        v = alpha**((1-p) // 2)

        self.alpha = alpha
        self.v = v

    def __repr__(self):
        return f"Special extremal {self.curve}"

    def twist(self):
        return SpecialExtremalCurve(self.curve.quadratic_twist(), self.i)

    def aut_i(self, P):
        
        # Input: P in E, possibly living over an extension field.
        # Output: the evaluation of the automorphism i at P.
        
        E = P.curve() 
        F = E.base_field()
        assert self.base_field.is_subring(F)
        assert F(E.a4()) == F(self.a4)
        
        if P == E(0):
            return P
        else:
            xP = P[0]
            yP = P[1]
            xQ = -xP
            yQ = F(self.i)*yP
            return E(xQ,yQ)

    def frob_p(self, P):

        # Input: P in E, possibly living over an extension field.
        # Output: the evaluation of the p-Frobenius endormorphism pi at P.
        
        E = P.curve()
        F = E.base_field()
        assert self.base_field.is_subring(F)
        assert F(E.a4()) == F(self.a4)
        
        p = self.p
        v = F(self.v)
        
        if P == E(0):
            return P
        else:           
            xP = P[0]
            yP = P[1]
            xQ = v**2 * xP**p
            yQ = v**3 * yP**p
            return E(xQ,yQ)
            
    def endo_basic(self, sigma, P):
        
        # Input: P in E, and sigma = [a,b,c,d].
        # Output: the evaluation of a + b*i + c*pi + d*i*pi at P.
        
        [a,b,c,d] = sigma
        iP = self.aut_i(P)
        jP = self.frob_p(P)
        kP = self.aut_i(jP)
        Q = a*P + b*iP + c*jP + d*kP
        return Q

    def endo_max(self, sigma, P):

        # Input: P in E, and sigma = [a,b,c,d].
        # Output: the evaluation of a + b*i + c*(i+pi)/2 + d*(1+i*pi)/2 at P.

        [a,b,c,d] = sigma
        E = P.curve()
        N = order(P)
        if N%2: # In this case, dividing by two is easy.
            inverse_of_2_mod_N = mod(2,N)**(-1)
            P_divided_by_two = ZZ(inverse_of_2_mod_N) * P
        else:
            divide_by_two_field = P.division_points(2,poly_only=True).splitting_field('X')
            E_ext = E.base_extend(divide_by_two_field)
            P_divided_by_two = (E_ext(P).division_points(2))[0]
        Q = self.endo_basic([2*a+d,2*b+c,c,d], P_divided_by_two)
        return E(Q)

    def endo_max_quat(self, beta, P):

        # Input: beta in O0, P in E
        # Output: the evaluation of beta at P
        
        [a,b,c,d] = beta.coefficient_tuple() # this is with respect to the basis 1, i, pi, i*pi
        [a,b,c,d] = [a-d,b-c,2*c,2*d] # converted into the basis 1, i, i+pi, 1+i*pi
        return self.endo_max([a,b,c,d], P)
    
    def find_kernel_basis_max(self, sigma, ell):

        # Find a basis P,Q of E[ell] such that P in ker(sigma) and Q in ker(sigma_dual)

        E = self.curve

        N = ZZ(sqrt(E.cardinality()))

        assert N%ell == 0
    
        [a,b,c,d] = sigma
        sigma_dual = [a+d,-b,-c,-d]
    
        while True:
            R = E.random_point()
            cofactor = N // ell
            R *= cofactor
            P = self.endo_max(sigma_dual, R)
            if order(P) == ell:
                break
        
        while True:
            R = E.random_point()
            cofactor = N // ell
            R *= cofactor
            Q = self.endo_max(sigma, R)
            if order(Q) == ell:
                break
    
        return P, Q

    def full_order_kernel_points(self, sigma, M):
    
        # We construct a basis P,Q for E[M]
        # such that P in ker(sigma) and Q in ker(sigma_dual).

        E = self.curve
        
        if M == 1:
            return E(0), E(0)
            
        tups = ZZ(M).factor()
        ell_powers = [ell**r for ell,r in tups]

        [a,b,c,d] = sigma
        sigma_dual = [a+d,-b,-c,-d]
        
        P = E(0)
        Q = E(0)
        
        for ell_power in ell_powers:
            P0,Q0 = self.find_kernel_basis_max(sigma, ell_power)
            P += P0
            Q += Q0
    
        assert P.order() == M
        assert Q.order() == M
        assert self.endo_max(sigma, P) == E(0)
        assert self.endo_max(sigma_dual, Q) == E(0)
        
        return P, Q

def push_orientation(E_home, E_base, negative_ideal, P0, M):

    # Input:
    # - a special extremal curve E_home
    # - a curve E_base
    # - a negative_ideal connecting E_home to E_base
    # - a point P0 of order M on E_home
    # such that:
    # (*) the subgroups corresponding to negative_ideal and P0 intersect trivially.
    # Output: the x-coordinate of the image of P0 under the isogeny corresponding to negative_ideal
    # Note: By (*), the corresponding image point still has order M

    Fq = E_home.base_field
    E0, iota, O0 = starting_curve(Fq)
    banned_integer = ZZ(radical(M))
    E_temp, psi, beta, J, _ = constructive_deuring(negative_ideal, E0, iota, banned_integer)
    assert E_temp.is_isomorphic(E_base)
    assert len(E_temp.isomorphisms(E_base)) == 2
    iso = E_temp.isomorphisms(E_base)[0]
    endo_ideal = multiplyIdeals(negative_ideal, J.conjugate(), beta=beta)
    assert endo_ideal.is_principal()
    assert endo_ideal == O0.left_ideal(beta)
    
    beta_P0 = E_home.endo_max_quat(beta, P0) # This lives on E_home, i.e. possibly a quadratic twist of E0

    xP = E_home.alpha**(-1) * beta_P0.x() # The x-coordinate of beta(P0) as a point on E0
    xP = xPoint(xP, E0) # Using the xPoint class to prepare evaluation under psi
    xP = xP.push(psi) # Now we are on E_temp
    xP = iso.x_rational_map()(xP.X) # Translate the point to E_base
    
    return xP

def multiple_point_KummerLineIsogeny(K, kernel_list, degree_list, eval_list = []):
    while kernel_list:
        kernel_point = kernel_list.pop(0)
        degree = degree_list.pop(0)
        if degree == 1:
            continue
        phi = KummerLineIsogeny(K, kernel_point, ZZ(degree))
        K = phi.codomain()
        kernel_list = [phi(P) for P in kernel_list]
        eval_list = [phi(P) for P in eval_list]
    if eval_list:
        return K, eval_list
    else:
        return K

def generate_cycle(p, r, Ms, Mt, sigma, trace):

    Fp = GF(p)
    R = PolynomialRing(Fp, 'X')
    X = R.gen()
    Fq = Fp.extension(X**2+1, 'i')
    i = Fq.gen()
    
    E0, iota, O0 = starting_curve(Fq)
    embedded_sigma = sum([sigma[k]*O0.basis()[k] for k in range(4)]) # sigma as an element of O0
    embedded_sigma_dual = embedded_sigma.conjugate()
    E_straight = SpecialExtremalCurve(E0, i)
    E_twist = E_straight.twist()
    alpha = E_twist.alpha

    P_straight, Q_straight = E_straight.full_order_kernel_points(sigma, Ms)
    P_twist, Q_twist = E_twist.full_order_kernel_points(sigma, Mt) # These points live on E_twist
    
    K_home = KummerLine(E_straight.curve)
    A = Fq(K_home.a())

    # check temp.csv for a backup
    backup_found = 0
    with open(current_folder / 'orientation_data/temp.csv', mode='r') as file:
        reader = csv.DictReader(file, delimiter=';')
        data = []
        found = False
        for row in reader:
            if p == int(row['p']) and r == int(row['r']) and Ms == int(row['Ms']) and Mt == int(row['Mt']) and trace == int(row['trace']):
                raw_cycle = ast.literal_eval(row['base_cycle'])
                aff_cycle = [[x[0]+x[1]*i if not x[0] == None else None for x in tup] for tup in raw_cycle]
                backup_found = 1
                print('Found backup in temp.csv')
                print('aff_cycle =', aff_cycle)
                
    if not backup_found:
        
        # initialize base curve
        
        aff_cycle = []
        raw_cycle = []
        
        if Mt == 1:
    
            Ps = K_home(P_straight[0])
            Qs = K_home(Q_straight[0])
    
            Pt = K_home.zero()
            Qt = K_home.zero()
            
            xPs = Ps.x()
            xQs = Qs.x()
            
            xPt = None
            xQt = None
    
        else:
            
            Ps, Pt = K_home(P_straight[0]), K_home(alpha**(-1) * P_twist[0]) # These points live on K_home (the Kummer line of E_straight)
            Qs, Qt = K_home(Q_straight[0]), K_home(alpha**(-1) * Q_twist[0]) # These points live on K_home (the Kummer line of E_straight)
            
            assert Mt * Pt == K_home.zero()
            assert Mt * Qt == K_home.zero()
            
            xPs, xPt = Ps.x(), Pt.x()
            xQs, xQt = Qs.x(), Qt.x()
            
        tup = [A, xPs, xPt, xQs, xQt]
        raw_tup = []
        
        for x in tup:
            if x == None:
                raw_tup.append([None, None])
            else:
                raw_tup.append([(Fq(x).polynomial())[0], (Fq(x).polynomial())[1]])
                
        aff_cycle.append(tup)
        raw_cycle.append(raw_tup)
        
        new_row = [p, r, trace, Ms, Mt, str(raw_cycle).replace(" ","")]
        with open(current_folder / 'orientation_data/temp.csv', mode='w', newline="") as file:
            writer = csv.writer(file, delimiter=";")
            writer.writerow(['p','r','trace','Ms','Mt','base_cycle'])
            writer.writerow(new_row)

    starting_length = len(aff_cycle) # this equals 1 if we did not find a backup

    [A, xPs, xPt, xQs, xQt] = aff_cycle[-1]

    K_home = KummerLine(EllipticCurve([0,A,0,1,0]))

    if Mt == 1:
        Ps, Pt = K_home(xPs), K_home.zero()
        Qs, Qt = K_home(xQs), K_home.zero()
    else:
        Ps, Pt = K_home(xPs), K_home(xPt)
        Qs, Qt = K_home(xQs), K_home(xQt)
    
    for e in range(starting_length, r):
    
        K_base, [Qs, Qt] = multiple_point_KummerLineIsogeny(K_home, [Ps, Pt], [Ms, Mt], [Qs, Qt])
        E_base = K_base.curve()
        A = K_base.a()
        
        negative_ideal = O0.left_ideal([(Ms * Mt)**(r-e), embedded_sigma_dual])
        
        if Mt == 1:
            
            xQs = Qs.x()
            xQt = None

            # split up the torsion for more efficient KLPT

            tups = ZZ(Ms).factor()

            Ms0 = prod([x[0]**x[1] for x in tups[0::2]])
            Ms1 = prod([x[0]**x[1] for x in tups[1::2]])

            P_straight0 = Ms1*P_straight
            P_straight1 = Ms0*P_straight

            xPs0 = push_orientation(E_straight, E_base, negative_ideal, P_straight0, Ms0)
            xPs1 = push_orientation(E_straight, E_base, negative_ideal, P_straight1, Ms1)

            xPs = (E_base.lift_x(xPs0) + E_base.lift_x(xPs1))[0]            
            xPt = None

            K_home = K_base

            Ps = K_home(xPs)
            Pt = K_home.zero()

        else:
            
            xQs = Qs.x()
            xQt = Qt.x()
                        
            xPs = push_orientation(E_straight, E_base, negative_ideal, P_straight, Ms)
            xPt = push_orientation(E_twist, E_base, negative_ideal, P_twist, Mt)

            K_home = K_base

            Ps, Pt = K_home(xPs), K_home(xPt)

        assert Ms * Ps == K_home.zero()
        assert Mt * Pt == K_home.zero()

        tup = [A, xPs, xPt, xQs, xQt]
        raw_tup = []
        for x in tup:
            if x == None:
                raw_tup.append([None, None])
            else:
                raw_tup.append([(Fq(x).polynomial())[0], (Fq(x).polynomial())[1]])

        aff_cycle.append(tup)
        raw_cycle.append(raw_tup)
        
        new_row = [p, r, trace, Ms, Mt, str(raw_cycle).replace(" ","")]
        with open(current_folder / 'orientation_data/temp.csv', mode='w', newline="") as file:
            writer = csv.writer(file, delimiter=";")
            writer.writerow(['p','r','trace','Ms','Mt','base_cycle'])
            writer.writerow(new_row)

    return raw_cycle

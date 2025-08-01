import sys
import os

from pathlib import Path
current_file = Path(__file__).resolve()
current_folder = current_file.parent
sys.path.append(str(current_folder.parent / 'external_modules/KummerIsogeny'))

from random import randint # For secret key generation
from sage.all import EllipticCurve, ZZ
from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny
from parameters import Parameters

params = None

class OrientedKummerLine:
    r"""
    Class for x-only oriented Montgomery curves
    """
    def __init__(self, K, Ps, Qs):
        self.K = K
        self.Ps = Ps
        self.Qs = Qs

    def __repr__(self):
        return f"Oriented {self.K}"

    def IsogenyAndPush(self, kernel_point, degree, check = False, threshold = 250):
#        if kernel_point.is_zero() or kernel_point.curve_point().order() != ZZ(degree):
#            print('problem found with ', kernel_point)
#            print('the expected degree is ', degree)
#            raise ValueError("kernel_point is not of the correct order")
        phi = KummerLineIsogeny(self.K, kernel_point, ZZ(degree), check = check, threshold = threshold)
        return OrientedKummerLine(phi.codomain(), phi(self.Ps), phi(self.Qs))
            
    def conjugate(self):
        return OrientedKummerLine(self.K, self.Qs, self.Ps)

    def curve(self):
        return self.K.curve()

    def a(self):
        return self.K.a()

    def check_order(self): # for debugging
        if self.Qs.is_zero():
            ordQs = 1
        else:
            ordQs = (self.Qs).curve_point().order()
        return ordQs

    def j_invariant(self):
        return self.curve().j_invariant()

def OrientedKummer_to_AffineRepresentation(OK):
    return OK.K.a(), OK.Ps.x(), OK.Qs.x()

def AffineRepresentation_to_OrientedKummer(A, xPs, xQs):
    Fq = params.base_field
    K = KummerLine(Fq, [A, 1])
    return OrientedKummerLine(K, K(xPs), K(xQs))

class SecretKey:
    r"""
    Class for exponent vectors
    """
    def __init__(self, vec_straight):
        self.straight = vec_straight

    def __repr__(self):
        return f"Secret key with exponent vectors {self.straight}"

    def reduce(self, exps_straight, k):
        # if (s_i > k*e_i), reduce modulo e_i (with output in [1,...,e_i])
        # if (s_i <= k*e_i), set to 0
        vec_straight_red = [((x - 1) % m + 1) * int((x // m) > k) for x, m in zip(self.straight, exps_straight)]
        return SecretKey(vec_straight_red)

def keygen(B):
    vec_straight = [randint(0, B*exp) for exp in params.exps_straight]
    return SecretKey(vec_straight)

def IsogenyWalk(OK_right, OK_left, ells, exps, vec, cofactor):
    n = len(ells)
    for i in range(n):
        ell = ells[i]
        e = exps[i]
        s = vec[i]
        for j in range(e):
            cofactor //= ell
            if j < s:
                # Walk right
                kernel_point = cofactor * OK_right.Ps
#                OK_temp = OK_right.IsogenyAndPush(kernel_point, ell)
#                if OK_temp.check_order() != params.Ms:
#                    print(OK_right.Ps, OK_right.Qs, cofactor, ell, OK_right.K.j_invariant())
#                    raise ValueError("problem")
                OK_right = OK_right.IsogenyAndPush(kernel_point, ell)
                OK_left.Ps *= ell
            else:
                # Walk left
                kernel_point = cofactor * OK_left.Ps
                OK_left = OK_left.IsogenyAndPush(kernel_point, ell)
                OK_right.Ps *= ell

    return OK_right, OK_left


def group_action_home_base(aff_home, aff_base, vec_straight):

    OK_home = AffineRepresentation_to_OrientedKummer(*aff_home)
    OK_base = AffineRepresentation_to_OrientedKummer(*aff_base)
    
    OK_right, OK_left = OK_home, OK_base.conjugate()
    OK_right, OK_left = IsogenyWalk(OK_right, OK_left, params.ells_straight, params.exps_straight, vec_straight, params.Ms)
    
    A_right = params.base_field(OK_right.a())
    A_left = params.base_field(OK_left.a())
    
    E_right = EllipticCurve([0,A_right,0,1,0])
    E_left = EllipticCurve([0,A_left,0,1,0])
    
    isos = E_left.isomorphisms(E_right)
    assert len(isos) > 0
    extra_aut = False
    if len(isos) > 2:
        extra_aut = True
        print("Extra automorphisms detected; possible correction required.")
    if extra_aut:
        it = 0
        for iso in isos:
            it += 1
            iso_x = iso.x_rational_map()
            xPs = iso_x(OK_left.Qs.x())
            Ps = OK_right.K(xPs)
            cofactor = params.Ms // 3
            if not cofactor * Ps == cofactor * OK_right.Qs:
                print("Problem solved after",it,"iterations.")
                break
    else:
        iso_x = isos[0].x_rational_map()
        xPs = iso_x(OK_left.Qs.x())
        
    xQs = OK_right.Qs.x()
    
    aff_out = [A_right, xPs, xQs]

    return aff_out


def group_action_square_free(aff_cycle, sk):

    # Input:
    # - a list of affine representations of a cycle of oriented curves E_0, E_1, ... , E_{r-1} such that E_{j+1} = [e_1,...,e_n] * E_j;
    # - a secret key sk, representing an ideal class [a] = [s_1,...,s_n], split into the part on the curves and on their twists
    # where 0 <= s_i <= e_i for all i.
    # Output:
    # - a list of the affine representations of the cycle of oriented curves [a]E_0, [a]E_1, ... , [a]E_{r-1}.

    r = len(aff_cycle)
    
    return [group_action_home_base(aff_cycle[e], aff_cycle[(e+1)%r], sk.straight) for e in range(r)]


def group_action(aff_cycle, sk, B):
    
    # Input:
    # - a list of affine representations of a cycle of oriented curves E_0, E_1, ... , E_{r-1} such that E_{j+1} = [e_1,...,e_n] * E_j;
    # - a secret key sk, representing an ideal class [a] = [s_1,...,s_n], split into the part on the curves and on their twists
    # where s_i >= 0 for all i.
    # Output:
    # - a list of the affine representations of the cycle of oriented curves [a]E_0, [a]E_1, ... , [a]E_{r-1}.

    for k in range(B):
        sk_squarefree = sk.reduce(params.exps_straight, k)
        print(sk_squarefree.straight)
        aff_cycle = group_action_square_free(aff_cycle, sk_squarefree)

    return aff_cycle

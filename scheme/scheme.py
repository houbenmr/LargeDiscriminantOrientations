import sys
import os

sys.path.append(os.path.abspath('../external_modules/KummerIsogeny'))

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
    def __init__(self, K, Ps, Pt, Qs, Qt):
        self.K = K
        self.Ps = Ps
        self.Pt = Pt
        self.Qs = Qs
        self.Qt = Qt

    def __repr__(self):
        return f"Oriented {self.K}"

    def IsogenyAndPush(self, kernel_point, degree, check = False, threshold = 250):
        phi = KummerLineIsogeny(self.K, kernel_point, ZZ(degree), check = check, threshold = threshold)
        return OrientedKummerLine(phi.codomain(), phi(self.Ps), phi(self.Pt), phi(self.Qs), phi(self.Qt))

    def twist(self):
        return OrientedKummerLine(self.K, self.Pt, self.Ps, self.Qt, self.Qs)

    def conjugate(self):
        return OrientedKummerLine(self.K, self.Qs, self.Qt, self.Ps, self.Pt)

    def curve(self):
        return self.K.curve()

    def a(self):
        return self.K.a()

    def j_invariant(self):
        return self.curve().j_invariant()

def OrientedKummer_to_AffineRepresentation(OK):
    return OK.K.a(), OK.Ps.x(), OK.Pt.x(), OK.Qs.x(), OK.Qt.x()

def AffineRepresentation_to_OrientedKummer(A, xPs, xPt, xQs, xQt):
    Fq = params.base_field
    K = KummerLine(Fq, [A, 1])
    return OrientedKummerLine(K, K(xPs), K(xPt), K(xQs), K(xQt))

class SecretKey:
    r"""
    Class for exponent vectors
    """
    def __init__(self, vec_straight, vec_twist):
        self.straight = vec_straight
        self.twist = vec_twist

    def __repr__(self):
        return f"Secret key with exponent vectors {self.straight} and {self.twist}"

def keygen(B):
    vec_straight = [randint(0, B*exp) for exp in params.exps_straight]
    vec_twist = [randint(0, B*exp) for exp in params.exps_twist]
    return SecretKey(vec_straight, vec_twist)

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
                OK_right = OK_right.IsogenyAndPush(kernel_point, ell)
                OK_left.Ps *= ell
            else:
                # Walk left
                kernel_point = cofactor * OK_left.Ps
                OK_left = OK_left.IsogenyAndPush(kernel_point, ell)
                OK_right.Ps *= ell
    return OK_right, OK_left


def group_action_home_base(aff_home, aff_base, vec_straight, vec_twist):

    OK_home = AffineRepresentation_to_OrientedKummer(*aff_home)
    OK_base = AffineRepresentation_to_OrientedKummer(*aff_base)
    
    OK_right, OK_left = OK_home, OK_base.conjugate()
    OK_right, OK_left = IsogenyWalk(OK_right, OK_left, params.ells_straight, params.exps_straight, vec_straight, params.Ms)
    
    OK_right, OK_left = OK_right.twist(), OK_left.twist()
    OK_right, OK_left = IsogenyWalk(OK_right, OK_left, params.ells_twist, params.exps_twist, vec_twist, params.Mt)
    
    A_right = params.base_field(OK_right.a())
    A_left = params.base_field(OK_left.a())
    
    E_right = EllipticCurve([0,A_right,0,1,0])
    E_left = EllipticCurve([0,A_left,0,1,0])
    
    assert E_left.is_isomorphic(E_right)
    #assert len(E_left.isomorphisms(E_right)) == 2
    iso = E_left.isomorphisms(E_right)[0]
    iso_x = iso.x_rational_map()
    
    xPt = iso_x(OK_left.Qs.x())
    xPs = iso_x(OK_left.Qt.x())
    xQt = OK_right.Qs.x()
    xQs = OK_right.Qt.x()
    
    aff_out = [A_right, xPs, xPt, xQs, xQt]

    return aff_out


def group_action_square_free(aff_cycle, sk):

    # Input:
    # - a list of affine representations of a cycle of oriented curves E_0, E_1, ... , E_{r-1} such that E_{j+1} = [e_1,...,e_n] * E_j;
    # - a secret key sk, representing an ideal class [a] = [s_1,...,s_n], split into the part on the curves and on their twists
    # where 0 <= s_i <= e_i for all i.
    # Output:
    # - a list of the affine representations of the cycle of oriented curves [a]E_0, [a]E_1, ... , [a]E_{r-1}.

    r = len(aff_cycle)
    
    return [group_action_home_base(aff_cycle[e], aff_cycle[(e+1)%r], sk.straight, sk.twist) for e in range(r)]


def group_action(aff_cycle, sk, B):
    
    # Input:
    # - a list of affine representations of a cycle of oriented curves E_0, E_1, ... , E_{r-1} such that E_{j+1} = [e_1,...,e_n] * E_j;
    # - a secret key sk, representing an ideal class [a] = [s_1,...,s_n], split into the part on the curves and on their twists
    # where s_i >= 0 for all i.
    # Output:
    # - a list of the affine representations of the cycle of oriented curves [a]E_0, [a]E_1, ... , [a]E_{r-1}.

    for k in range(B):
        sk_squarefree = SecretKey([max(sk.straight[i] - k * params.exps_straight[i], 0) % params.exps_straight[i] for i in range(len(params.exps_straight))],\
                                  [max(sk.twist[i] - k * params.exps_twist[i], 0) % params.exps_twist[i] for i in range(len(params.exps_twist))])
        aff_cycle = group_action_square_free(aff_cycle, sk_squarefree)

    return aff_cycle

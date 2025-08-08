import sys

from pathlib import Path
current_file = Path(__file__).resolve()
current_folder = current_file.parent
sys.path.append(str(current_folder.parent / 'external_modules/KummerIsogeny'))

from random import randint # For secret key generation
from sage.all import EllipticCurve, ZZ
from kummer_line import KummerLine, KummerPoint
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
        if kernel_point.is_zero(): #or kernel_point.curve_point().order() != ZZ(degree):
            print('problem found with ', kernel_point)
            print('the expected degree is ', degree)
            raise ValueError("kernel_point is not of the correct order")
        phi = KummerLineIsogeny(self.K, kernel_point, ZZ(degree), check = check, threshold = threshold)        

        # Computing the multiplier (i.e. action on invariant differentials)
        # Currently quite costly. We should be able to obtain the multiplier
        # from the isogeny computation; this may be trickier in case of sqrtVelu

        # The multiplier is stored projectively 

        assert degree%2
        d = degree // 2

        pi_X = 1
        pi_Z = 1
        K_muls = kernel_point.multiples()
        for _ in range(d):
            Ki = next(K_muls)
            X, Z = Ki.XZ()
            pi_X *= X
            pi_Z *= Z

        # pi = prod([(i*kernel_point).x() for i in range(1, d+1)])

        # We invert pi to get the actual multiplier

        return OrientedKummerLine(phi.codomain(), phi(self.Ps), phi(self.Pt), phi(self.Qs), phi(self.Qt)), pi_Z, pi_X

    def twist(self):
        return OrientedKummerLine(self.K, self.Pt, self.Ps, self.Qt, self.Qs)

    def conjugate(self):
        return OrientedKummerLine(self.K, self.Qs, self.Qt, self.Ps, self.Pt)

    def curve(self):
        return self.K.curve()

    def AC(self):
        return self.K.extract_constants()

    def a(self):
        return self.K.a()

    def check_order(self): # for debugging; currently only supports straight points
        if self.Qs.is_zero():
            ordQs = 1
        else:
            ordQs = (self.Qs).curve_point().order()
        return ordQs == params.Ms

    def j_invariant(self):
        return self.curve().j_invariant()

def OrientedKummer_to_AffineRepresentation(OK):
    return OK.K.a(), OK.Ps.x(), OK.Pt.x(), OK.Qs.x(), OK.Qt.x()

def AffineRepresentation_to_OrientedKummer(A, xPs, xPt, xQs, xQt):
    Fq = params.base_field
    K = KummerLine(Fq, [A, 1])
    return OrientedKummerLine(K, K(xPs), K(xPt), K(xQs), K(xQt))

# Functions to convert between affine and projective representations

def OrientedKummer_to_AffineRepresentation(OK):
    return OK.K.a(), OK.Ps.x(), OK.Pt.x(), OK.Qs.x(), OK.Qt.x()

def AffineRepresentation_to_OrientedKummer(A, xPs, xQs):
    K = KummerLine(params.base_field, [A, 1])
    return OrientedKummerLine(K, K(xPs), K(xPt), K(xQs), K(xQt))

def OrientedKummer_to_ProjectiveRepresentation(OK):
    return *OK.K.AC(), *OK.Ps.XZ(), *OK.Pt.XZ(), *OK.Qs.XZ(), *OK.Qt.XZ()

def ProjectiveRepresentation_to_OrientedKummer(A_X, A_Z, Ps_X, Ps_Z, Pt_X, Pt_Z, Qs_X, Qs_Z, Qt_X, Qt_Z):
    K = KummerLine(params.base_field, [A_X, A_Z])
    Ps = KummerPoint(K, [Ps_X, Ps_Z])
    Pt = KummerPoint(K, [Pt_X, Pt_Z])
    Qs = KummerPoint(K, [Qs_X, Qs_Z])
    Qt = KummerPoint(K, [Qt_X, Qt_Z])
    return OrientedKummerLine(K, Ps, Pt, Qs, Qt)

def aff_cycle_to_proj_cycle(aff_cycle):
    r = len(aff_cycle)
    out_cycle = []
    for e in range(r):
        [A, xPs, xPt, xQs, xQt] = aff_cycle[e]
        out_cycle.append([A, 1, xPs, 1, xPt, 1, xQs, 1, xQt, 1])
    return out_cycle

def proj_cycle_to_aff_cycle(proj_cycle):
    r = len(proj_cycle)
    out_cycle = []
    for e in range(r):
        [A_X, A_Z, Ps_X, Ps_Z, Pt_X, Pt_Z, Qs_X, Qs_Z, Qt_X, Qt_Z] = proj_cycle[e]
        out_cycle.append([A_X/A_Z, Ps_X/Ps_Z, Pt_X/Pt_Z, Qs_X/Qs_Z, Qt_X/Qt_Z])
    return out_cycle

class SecretKey:
    r"""
    Class for exponent vectors
    """
    def __init__(self, vec_straight, vec_twist):
        self.straight = vec_straight
        self.twist = vec_twist

    def __repr__(self):
        return f"Secret key with exponent vectors {self.straight} and {self.twist}"

    def reduce(self, exps_straight, exps_twist, k):
        # if (s_i > k*e_i), reduce modulo e_i (with output in [1,...,e_i])
        # if (s_i <= k*e_i), set to 0
        vec_straight_red = [max(min(x - k*m, m), 0) for x, m in zip(self.straight, exps_straight)]
        vec_twist_red = [max(min(x - k*m, m), 0) for x, m in zip(self.twist, exps_twist)]
        return SecretKey(vec_straight_red, vec_twist_red)

    def min(self, other):
        vec_straight = [x-y for x, y in zip(self.straight, other.straight)]
        vec_twist = [x-y for x, y in zip(self.twist, other.twist)]
        return SecretKey(vec_straight, vec_twist)

def keygen(B):
    vec_straight = [randint(0, B*exp) for exp in params.exps_straight]
    vec_twist = [randint(0, B*exp) for exp in params.exps_twist]
    return SecretKey(vec_straight, vec_twist)

def IsogenyWalk(OK_right, OK_left, ells, exps, vec, cofactor):
    mult_sq_X = cofactor
    mult_sq_Z = 1
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
                OK_right, pi_X, pi_Z = OK_right.IsogenyAndPush(kernel_point, ell)
                mult_sq_X *= pi_X**2
                mult_sq_Z *= ell * pi_Z**2
                OK_left.Ps *= ell
            else:
                # Walk left
                kernel_point = cofactor * OK_left.Ps
                OK_left, pi_X, pi_Z = OK_left.IsogenyAndPush(kernel_point, ell)
                mult_sq_Z *= pi_X**2
                mult_sq_X *= ell * pi_Z**2
                OK_right.Ps *= ell

    return OK_right, OK_left, mult_sq_X, mult_sq_Z # mult_sq = (pi_right*M_left/pi_left)^2

def group_action_home_base(proj_home, proj_base, vec_straight, vec_twist):

    OK_home = ProjectiveRepresentation_to_OrientedKummer(*proj_home)
    OK_base = ProjectiveRepresentation_to_OrientedKummer(*proj_base)
   
    OK_right, OK_left = OK_home, OK_base.conjugate()
    OK_right, OK_left, mult_sq_X_s, mult_sq_Z_s = IsogenyWalk(OK_right, OK_left, params.ells_straight, params.exps_straight, vec_straight, params.Ms)
    
    OK_right, OK_left = OK_right.twist(), OK_left.twist()
    OK_right, OK_left, mult_sq_X_t, mult_sq_Z_t = IsogenyWalk(OK_right, OK_left, params.ells_twist, params.exps_twist, vec_twist, params.Mt)

    mult_sq_X = mult_sq_X_s * mult_sq_X_t
    mult_sq_Z = mult_sq_Z_s * mult_sq_Z_t    
 
    A_right_X, A_right_Z = OK_right.AC()
    A_left_X, A_left_Z = OK_left.AC()

    Ps_X, Ps_Z = OK_left.Qt.XZ() # Not corrected for isomorphism yet! 
    Pt_X, Pt_Z = OK_left.Qs.XZ() # "--------------------------------"

    Qs_X, Qs_Z = OK_right.Qt.XZ()
    Qt_X, Qt_Z = OK_right.Qs.XZ()
    
    proj_out = [A_right_X, A_right_Z, Ps_X, Ps_Z, Pt_X, Pt_Z, Qs_X, Qs_Z, Qt_X, Qt_Z]
   
    return proj_out, mult_sq_X, mult_sq_Z, A_left_X, A_left_Z


def group_action_square_free(proj_cycle, sk, alpha_sq):

    # Input:
    # - a list of affine representations of a cycle of oriented curves E_0, E_1, ... , E_{r-1} such that E_{j+1} = [e_1,...,e_n] * E_j;
    # - a secret key sk, representing an ideal class [a] = [s_1,...,s_n], split into the part on the curves and on their twists
    # where 0 <= s_i <= e_i for all i.
    # Output:
    # - a list of the affine representations of the cycle of oriented curves [a]E_0, [a]E_1, ... , [a]E_{r-1}.

    r = len(proj_cycle)
    mult_sq_total_X = 1
    mult_sq_total_Z = 1
    out_cycle = [[] for e in range(r)]
    
    for e in range(r):
        home = proj_cycle[e]
        base = proj_cycle[(e+1)%r]
        out_cycle[e], mult_sq_X, mult_sq_Z, A_left_X, A_left_Z = group_action_home_base(home, base, sk.straight, sk.twist)
        mult_sq_total_X *= mult_sq_X
        mult_sq_total_Z *= mult_sq_Z
            
    # Now we correct for an isomorphism
    [A_right_X, A_right_Z, Ps_X, Ps_Z, Pt_X, Pt_Z] = out_cycle[r-1][:6]
    
    u_sq_X = mult_sq_total_X
    u_sq_Z = alpha_sq * mult_sq_total_Z
    # u_sq = mult_sq / alpha_sq

    r_X = u_sq_X * A_left_Z * A_right_X - u_sq_Z * A_left_X * A_right_Z
    r_Z = 3 * u_sq_Z * A_left_Z * A_right_Z
    # r = (u_sq*A_right - A_left)/3

    # The Weierstrass isomorphism E_right -> E_left is x |-> u_sq*x + r
    
    Ps_X = u_sq_Z * (Ps_X * r_Z - r_X * Ps_Z)
    Ps_Z = u_sq_X * Ps_Z * r_Z
    # xPs = u_sq**(-1) * ((OK_left.Qs.x()) - r)

    Pt_X = u_sq_Z * (Pt_X * r_Z - r_X * Pt_Z)
    Pt_Z = u_sq_X * Pt_Z * r_Z
    # xPt = u_sq**(-1) * ((OK_left.Qt.x()) - r)

    out_cycle[r-1][:6] = [A_right_X, A_right_Z, Ps_X, Ps_Z, Pt_X, Pt_Z]
    
    return out_cycle

def group_action(aff_cycle, sk, B):
    
    # Input:
    # - a list of affine representations of a cycle of oriented curves E_0, E_1, ... , E_{r-1} such that E_{j+1} = [e_1,...,e_n] * E_j;
    # - a secret key sk, representing an ideal class [a] = [s_1,...,s_n], split into the part on the curves and on their twists
    # where s_i >= 0 for all i.
    # Output:
    # - a list of the affine representations of the cycle of oriented curves [a]E_0, [a]E_1, ... , [a]E_{r-1}.
    
#    skt = sk # for printing
    
    proj_cycle = aff_cycle_to_proj_cycle(aff_cycle)

    for k in range(B):
        sk_squarefree = sk.reduce(params.exps_straight, params.exps_twist, k)
#        skt = skt.min(sk_squarefree)
#        print(skt)
        proj_cycle = group_action_square_free(proj_cycle, sk_squarefree, params.alpha_sq)

    aff_cycle = proj_cycle_to_aff_cycle(proj_cycle)

    return aff_cycle

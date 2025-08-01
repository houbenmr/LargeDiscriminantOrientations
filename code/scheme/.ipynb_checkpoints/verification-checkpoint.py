import sys
import os
import time

from pathlib import Path
current_file = Path(__file__).resolve()
sys.path.append(str(current_file.parent.parent / 'external_modules/KummerIsogeny'))

from sage.all import ZZ
from math import prod
from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny

params = None

def multiple_point_KummerLineIsogeny(K, kernel_list, degree_list, eval_list = []):
    while kernel_list:
        kernel_point = kernel_list.pop(0)
        degree = degree_list.pop(0)
        phi = KummerLineIsogeny(K, kernel_point, ZZ(degree))
        K = phi.codomain()
        kernel_list = [phi(P) for P in kernel_list]
        eval_list = [phi(P) for P in eval_list]
    if eval_list:
        return K, eval_list
    else:
        return K

def verify_aff_repr(aff_repr):   
    [A, xPs, xPt, xQs, xQt] = aff_repr
    K = KummerLine(params.base_field, [A, 1])
    assert params.Ms * K(xPs) == K.zero()
    assert params.Mt * K(xPt) == K.zero()
    assert params.Ms * K(xQs) == K.zero()
    assert params.Mt * K(xQt) == K.zero()

def verify_cycle(aff_cycle):
    
    # checks that E_{j+1} = [e_1, ... ,e_n] E_j for all j

    print('verifying cycle of affine representations...')
    r = len(aff_cycle)
    
    for e in range(r):
        A, Ps, Pt, Qs, Qt = aff_cycle[e]
        K = KummerLine(params.base_field, [A,1])
        assert multiple_point_KummerLineIsogeny(K, [K(Ps), K(Pt)], [params.Ms, params.Mt]).j_invariant() == KummerLine(params.base_field, [aff_cycle[(e+1)%r][0], 1]).j_invariant()
        assert multiple_point_KummerLineIsogeny(K, [K(Qs), K(Qt)], [params.Ms, params.Mt]).j_invariant() == KummerLine(params.base_field, [aff_cycle[(e-1)%r][0], 1]).j_invariant()

def verify_M(ells_straight, ells_twist, exps_straight, exps_twist, Ms, Mt):
    assert Ms == prod([ells_straight[i]**exps_straight[i] for i in range(len(ells_straight))])
    assert Mt == prod([ells_twist[i]**exps_twist[i] for i in range(len(ells_twist))])
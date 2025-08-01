from sage.all import PolynomialRing, GF, ZZ, odd_part, factor # For Fq-generation and integer factorization
from math import prod
import csv
import ast  # To safely evaluate string representations of lists

from pathlib import Path
current_file = Path(__file__).resolve()
current_folder = current_file.parent

from scheme_cost import minimize_scheme_cost
from orientation_tools import find_sigma, generate_cycle

class Parameters:
    def __init__(self, p, r, Ms, Mt, trace, raw_cycle):

        M = Ms * Mt
        tups = ZZ(M).factor()
        ells = [ell for [ell,r] in tups]
        exps = [r for [ell,r] in tups]
        ells_straight = [ell for ell in ells if (p+1)%ell == 0]
        ells_twist = [ell for ell in ells if (p-1)%ell == 0]
        exps_straight = [exps[i] for i in range(len(exps)) if (p+1)%ells[i] == 0]
        exps_twist = [exps[i] for i in range(len(exps)) if (p-1)%ells[i] == 0]

        Fp = GF(p)
        R = PolynomialRing(Fp, 'X')
        X = R.gen()
        Fq = Fp.extension(X**2+1, 'i')
        i = Fq.gen()

        base_cycle = [[x[0]+x[1]*i for x in tup if not x[0] == None] for tup in raw_cycle]

        self.base_field = Fq
        self.p = p
        self.r = r
        self.Ms = Ms
        self.Mt = Mt
        self.M = M
        self.trace = trace
        self.ells = ells
        self.exps = exps
        self.ells_straight = ells_straight
        self.ells_twist = ells_twist
        self.exps_straight = exps_straight
        self.exps_twist = exps_twist
        self.raw_cycle = raw_cycle
        self.base_cycle = base_cycle
    
    def __repr__(self):
        return f"Scheme parameters for p = {self.p} and r = {self.r} using (p+1)-torsion {ZZ(self.Ms).factor()} and (p-1)-torsion {ZZ(self.Mt).factor()}"
    
def setup(p, r, Ms = 0, Mt = 0, trace = 0, search_database = True, generate = False, write = False, lmbda = 221):

    if not p%4 == 3:
        raise ValueError("p must be 3 mod 4")
    
    if Ms or Mt:
        if not (Ms and Mt):
            raise ValueError("Ms and Mt must either both be specified or both be unspecified")
            
        print('Ms and Mt are prespecified')
    
        if (p+1)%Ms:
            raise ValueError("Ms does not divide p+1")
    
        if (p-1)%Mt:
            raise ValueError("Mt does not divide p-1")

    else:
        print('Ms and Mt are not prespecified; computing optimal values for security level', lmbda)
        
        tups = ZZ(p**2-1).odd_part().factor()
        ells = [ell for [ell,r] in tups]
        exps = [r for [ell,r] in tups]
        best, exps, B = minimize_scheme_cost(ells, exps, lmbda)
        ells_straight = [ell for ell in ells if (p+1)%ell == 0]
        ells_twist = [ell for ell in ells if (p-1)%ell == 0]
        exps_straight = [exps[i] for i in range(len(exps)) if (p+1)%ells[i] == 0]
        exps_twist = [exps[i] for i in range(len(exps)) if (p-1)%ells[i] == 0]
        Ms = prod([ells_straight[i]**exps_straight[i] for i in range(len(ells_straight))])
        Mt = prod([ells_twist[i]**exps_twist[i] for i in range(len(ells_twist))])
    
    M = Ms * Mt
    
    if M%2 == 0:
        raise ValueError("can only handle odd-degree isogenies")

    print('Using Ms =', ZZ(Ms).factor(), 'and Mt =', ZZ(Mt).factor())
    
    if search_database:
        with open(current_folder / 'orientation_data/orientation_data.csv', mode='r') as file:
            reader = csv.DictReader(file, delimiter=';')
            data = []
            found = False
            for row in reader:
                if p == int(row['p']) and r == int(row['r']) and Ms == int(row['Ms']) and Mt == int(row['Mt']):
                    trace = int(row['trace'])
                    raw_cycle = ast.literal_eval(row['base_cycle'])
                    print('Found scheme parameters in database.')
                    return Parameters(p, r, Ms, Mt, trace, raw_cycle)
                    
        print('The data p, r, Ms, Mt does not appear in the database.')

    if generate:
        
        print('Generating scheme parameters...')

        sigma, trace = find_sigma(p, r, M, trace = trace)

        raw_cycle = generate_cycle(p, r, Ms, Mt, sigma, trace)

        print('found base cycle', raw_cycle)

        if write:
            new_row = [p, r, trace, Ms, Mt, str(raw_cycle).replace(" ","")]
            with open(current_folder / 'orientation_data/orientation_data.csv', mode='a', newline="") as file:
                writer = csv.writer(file, delimiter=";")
                writer.writerow(new_row)

        return Parameters(p, r, Ms, Mt, trace, raw_cycle)

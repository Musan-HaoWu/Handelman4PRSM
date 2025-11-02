from sympy import symbols, simplify, Abs
from itertools import product
from util.basic import *
from util.handelman import check_compatible, get_monoid
from gurobipy import Model, GRB


class Entailment:
    def __init__(self, vars, gs, f):
        self.K = gs
        self.f = f
        self.vars = vars

    def remove_absolute_values(self):
        """
        Takes an expression with possible absolute values 
        and removes them by splitting the cases.
        For example:
            Input: a*Abs(x-1) + b*Abs(x) + c, 
            Output: 
                [
                ([x-1, x], a*(x-1) + b*x + c),
                ([1-x, x], a*(1-x) + b*x + c),
                ([x-1, -x], a*(x-1) - b*x + c),
                ([1-x, -x], a*(1-x) - b*x + c)
                ]
            Nested example:
            Input: Abs(x-a*Abs(y))
            Output:
                [
                ([y, x-a*y], x-a*y),
                ([-y, x+a*y], x+a*y),
                ([y, a*y-x], a*y-x),
                ([-y, -x-a*y], -x-a*y)
                ]
        """
        cases = []

        # Find all absolute value terms in the expression
        # Note: the order f.atoms() returns is random!
        abs_terms = list(sorted(self.f.atoms(Abs), key=str, reverse=True) )  
        
        if not abs_terms:
            return [self]  # No absolute values, return the current entailment as is

        # Process the first absolute value term
        abs_term = abs_terms[0]
        inner_expr = abs_term.args[0]  # Extract the inner expression of Abs

        # Create two cases: one for the positive and one for the negative of the inner expression
        positive_case = self.f.subs(abs_term, inner_expr)
        negative_case = self.f.subs(abs_term, -inner_expr)

        # Create new entailments for each case
        cases.append(Entailment(self.vars, self.K + [inner_expr], positive_case))
        cases.append(Entailment(self.vars, self.K + [-inner_expr], negative_case))

        # Recursively handle remaining absolute values in each case
        result = []
        for case in cases:
            result.extend(case.remove_absolute_values())

        return result

    def remove_denominator(self):
        """
        If f is a rational expression, then remove its denominator.
        For example:
            Input: Abs(x)/(1+Abs(x))
            Output: Abs(x)

            Input: a*Abs(x)+1/(Abs(x)+1) * x/(2*x+1), 
            Output: a*Abs(x)*x+x
        """

        expr = self.f
        # Find all terms that are denominators in the expression
        numerator, denominator = expr.as_numer_denom()
        # Simplify the resulting expression
        expr = simplify(numerator)

        # Update the entailment's function
        self.f = expr

    def remove_duplicate(self):
        """
        Remove duplicate constraints in K
        """
        unique_constraints = list(set(self.K))
        self.K = unique_constraints

    def check_K_empty(self):
        """
            to check if  K is empty using Farkas' lemma
        """
        if len(self.K) == 0:
            return False
        
        monomial_list = generate_monomials(self.vars, degree=1)
        K_cmatrix = get_coeff_matrix(self.K, self.vars, monomial_list)
        not_empty = check_compatible(K_cmatrix)
        empty = not not_empty
        if self.f == 0:
            empty = True
        return empty


    def homogenize(self, x0):
        new_K = []
        for g in self.K:
            g_h = homogenization(g, self.vars, x0)
            new_K.append(g_h)
        
        f_h = homogenization(self.f, self.vars, x0)
        equality = -1 + x0
        for var in self.vars:
            if var in new_K:
                equality += var
            if -var in new_K:
                equality -= var
        
        self.K = new_K + [x0, equality, -equality]
        self.f = f_h
        self.vars = self.vars + [x0]
     

class Constrainst:
    def __init__(self, variables, parameters):
        self.constraints = []
        self.variables = variables # system variables
        self.parameters = parameters # unknown coeffcients in v
        self.m = Model("V_Constraints")
        self.__prepare_gurobi_variables()
        self.x0 = symbols('x0')  # homogenization variable

    def add_constraint(self, entailment):
        if isinstance(entailment, Entailment):
           self.constraints.append(entailment)
        else:
            raise TypeError("Only Entailment instances can be added as constraints.")

    def print_constraints(self):
        print("Number of constraints:", len(self.constraints))
        for i, constraint in enumerate(self.constraints):
            print(f"Constraint {i+1}: ", constraint.K, ">= 0 ==>", simplify(constraint.f), ">=0")

    def __remove_absolute_values_all(self):
        new_constraints = []
        for entailment in self.constraints:
            cases = entailment.remove_absolute_values()
            new_constraints.extend(cases)
        self.constraints = new_constraints

    def __remove_denominators_all(self):
        for entailment in self.constraints:
            entailment.remove_denominator()

    def __remove_duplicates_all(self):
        for entailment in self.constraints:
            entailment.remove_duplicate()
    
    def rewrite(self):
        self.__remove_denominators_all()
        self.__remove_absolute_values_all()
        self.__remove_duplicates_all()
        # self.__homogenize_all()

    def remove_empty_K(self):
        new_constraints = []
        for entailment in self.constraints:
            if not entailment.check_K_empty():
                new_constraints.append(entailment)
        self.constraints = new_constraints

    def __prepare_gurobi_variables(self):
        v_coeffs_gvars = []
        for i, coeff in enumerate(self.parameters):
            v_coeffs_gvars.append(self.m.addVar(name=f'vc{i}', lb=-GRB.INFINITY, ub=GRB.INFINITY))
        
        self.m.update() 
        self.gurobi_parameters = v_coeffs_gvars

    def sympy2gurobi(self, sympy_expr):
        expr_str = str(sympy_expr)
        for i in range(len(self.parameters)):
            expr_str = expr_str.replace(str(self.parameters[i]), f"self.gurobi_parameters[{i}]")
        expr = eval(expr_str)
        return expr
    
    def translate_handelman(self, D):
        """
        Adds Handelman constraints to the Gurobi model for polynomial p over the 
        constraints K up to the specified maximum degree.

        Args:
            D (int): Maximum degree for the Handelman monoid.
        """
        for e, entailment in enumerate(self.constraints):
            K = entailment.K
            f = entailment.f
            monomial_basis = generate_monomials(self.variables, degree=D)
            
            # Generate the Handelman monoid basis
            handelman_basis = get_monoid(K, D)
            # Get the coefficient matrix cmatrix for the Handelman basis
            cmatrix = get_coeff_matrix(handelman_basis, self.variables, monomial_basis)
            # Add constraints: cmatrix^T * lambda = p_coeff_vector
            cmatrix_T = cmatrix.T 
            # Define m_constraints as the number of rows in cmatrix
            monoid_rows = cmatrix.shape[0]

            # Get the coefficient vector of p with respect to the monomial basis
            f_coeffs = get_coeff_vector(f, self.variables, monomial_basis)

            lambdas = self.m.addVars(monoid_rows, name=f"lambda{e}", lb=0.0)
            for j in range(len(monomial_basis)):
                expr = sum(lambdas[i] * cmatrix_T[j, i] for i in range(monoid_rows))
                self.m.addConstr(expr == self.sympy2gurobi(f_coeffs[j]))

    def solve_for_v(self):
        self.m.setParam('OutputFlag', 0)
        self.m.setObjective(0, GRB.MINIMIZE)
        self.m.optimize()
        if self.m.status == GRB.OPTIMAL or self.m.status == GRB.SUBOPTIMAL:
            v_solution = [var.X for var in self.gurobi_parameters]
            return True, v_solution
        else:
            return False, None
        
    def __homogenize_all(self):
        self.variables = self.variables + [self.x0]
        for entailment in self.constraints:
            entailment.homogenize(self.x0)



if __name__ == "__main__":
    # tests
    x, y = symbols('x, y') 
    vars = [x, y]
    a = symbols('a')


    # K1 = [x+10, 10-x] 
    # f1 = 2*x
    # entailment1 = Entailment(vars, K1, f1)
    # constraints.add_constraint(entailment1)

    # K2 = [x-2, 3-x, y-2, 3-y]
    # f2 = Abs(x - Abs(1-y)) + 1
    # entailment2 = Entailment(vars, K2, f2)
    # constraints.add_constraint(entailment2)

    # K3 = [x+10, 10-x]
    # f3 = 2*Abs(3*x-1)/(1+Abs(x))*x/(2*x+1)
    # entailment3 = Entailment(vars, K3, f3)
    # constraints.add_constraint(entailment3)

    # K4 = [x+10, 10-x]
    # f4 = y+1/(1+Abs(a/(1+Abs(x))))
    # entailment4 = Entailment(vars, K4, f4)
    # constraints.add_constraint(entailment4)

    # constraints.print_constraints()

    # print("After removing denominators")
    # constraints.remove_denominators_all()
    # constraints.print_constraints()

    # print("After removing absolute values:")
    # constraints.remove_absolute_values_all()
    # constraints.print_constraints()

    v, v_coeffs, monomials = template_poly(vars, degree=2, name='v') 

    constraints = Constrainst(variables=vars, parameters=v_coeffs)

    K5 = [x, 1-x]
    entailment5 = Entailment(vars, K5, v)
    constraints.add_constraint(entailment5)

    constraints.print_constraints()
     
    x0 = symbols('x0') 
    constraints.homogenize_all()
    constraints.print_constraints()
import numpy as np
from util.basic import *
from gurobipy import Model, GRB, GurobiError
from itertools import combinations_with_replacement
from collections import Counter
# from lp import *

def check_empty(cmatrix):
    """
    Solves the LP to check if 1 can be represented as a non-negative linear 
    combination of the constraints whose coefficients are in cmatrix.
    
    The problem is formulated as finding lambda >= 0 such that: 
    cmatrix^T * lambda = b_target = [1, 0, 0, ...]
    """
    # The target vector b_target = [1, 0, 0, ...]
    # This means: Coeff of '1' must be 1.0, Coeff of 'x', 'y', etc. must be 0.
    
    if cmatrix.size == 0:
        return False
    
    m_constraints, n_monomials = cmatrix.shape

    b_target = np.zeros(n_monomials)
    b_target[0] = -1.0 

    if n_monomials != len(b_target):
        return False

    
    model = Model("FarkasCompatibilityLP")
    model.setParam('OutputFlag', 0) 
    # Define the lambda variables (non-negative multipliers)
    lambdas = model.addVars(m_constraints, name="lambda", lb=0.0)
    lambda_c = model.addVar(name="lambda0", lb=0.0)
    # Dummy objective (feasibility problem)
    model.setObjective(0, GRB.MINIMIZE)
    # Add constraints: cmatrix^T * lambda = b_target
    cmatrix_T = cmatrix.T 
    
    for j in range(n_monomials):
        # Sum_{i=1}^{m} lambda_i * cmatrix_T[j, i] == b_target[j]
        expr = sum(lambdas[i] * cmatrix_T[j, i] for i in range(m_constraints))
        if j==1:
            expr += lambda_c
        model.addConstr(expr == b_target[j])
    model.optimize()
    if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
        # multipliers = [lambdas[i].X for i in range(m_constraints)]
        # print(multipliers)
        empty = True
        return empty
    else:
        # print("Incompatible: No solution found.")
        empty = False
        return empty


def get_monoid(g_constraints, d):
    """
    Generates the elements of the Handelman Monoid H_d, which are products of 
    powers of the polynomials in g_constraints, such that the sum of the 
    exponents is less than or equal to d.

    Args:
        g_constraints (list): A list of SymPy polynomial expressions (g1, g2, ...).
        d (int): The maximum total degree of the exponent sum (e1 + e2 + ...).

    Returns:
        list: A list of SymPy expressions representing the monoid elements.
    """
    
    # Ensure the input is a list (or iterable)
    g_constraints = list(g_constraints)
    monoid_elements = []
    
    # Iterate through all possible sum of exponents (degree k) from 0 to d
    for k in range(d + 1):
        # combinations_with_replacement generates all ways to choose k items 
        # from g_constraints with replacement.
        # E.g., for (g1, g2) and k=2, it gives: (g1, g1), (g1, g2), (g2, g2)
        
        for combination in combinations_with_replacement(g_constraints, k):
            # The Counter counts the frequency of each constraint in the product.
            # E.g., for (g1, g1, g2), the counter is {g1: 2, g2: 1}
            exponent_map = Counter(combination)
            
            # Start the product element as 1 (for the degree 0 case)
            product_element = 1
            
            # Build the product expression: g1**e1 * g2**e2 * ...
            for g in g_constraints:
                # Get the exponent for the current constraint g. Defaults to 0.
                exponent = exponent_map.get(g, 0)
                product_element *= g**exponent
            
            # Use simplify() to ensure terms are expanded if necessary for clarity, 
            # though the expressions themselves (g1**2, g1*g2, etc.) are often preferred
            # in the context of Handelman's basis. We will keep it unexpanded 
            # as a product for a cleaner basis representation.
            monoid_elements.append(product_element)
            
    return monoid_elements

def add_handelman_constraints(model, K, f, variables, D):
    """
    Adds Handelman constraints to the Gurobi model for polynomial p over the 
    constraints K up to the specified maximum degree.

    Args:
        model (gurobipy.Model): The Gurobi optimization model.
        K (list): List of SymPy polynomial expressions representing constraints.
        p (sympy.Expr): The polynomial expression to be constrained.
        variables (list): List of SymPy symbols representing variables.
        D (int): Maximum degree for the Handelman monoid.
    """
    # Generate the Handelman monoid basis
    handelman_basis = get_monoid(K, D)

    # Get the coefficient vector of p with respect to the monomial basis
    monomial_basis = generate_monomials(variables, degree=D)
    f_coeffs = get_coeff_vector(f, variables, monomial_basis)

    # Get the coefficient matrix cmatrix for the Handelman basis
    cmatrix = get_coeff_matrix(handelman_basis, variables, monomial_basis)

    # Add constraints: cmatrix^T * lambda = p_coeff_vector
    cmatrix_T = cmatrix.T 

    # Define m_constraints as the number of rows in cmatrix
    monoid_rows = cmatrix.shape[0]

    lambdas = model.addVars(monoid_rows, name="lambda", lb=0.0)
    for j in range(len(monomial_basis)):
        expr = sum(lambdas[i] * cmatrix_T[j, i] for i in range(monoid_rows))
        model.addConstr(expr == f_coeffs[j])

    # Check compatibility using Farkas lemma
    
    # compatible, multipliers = check_compatible_farkas(cmatrix, monomial_basis)

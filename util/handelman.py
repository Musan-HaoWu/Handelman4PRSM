from util.basic import generate_monomials, get_coeff_vector, get_coeff_matrix
import numpy as np
from gurobipy import Model, GRB, GurobiError
from itertools import combinations_with_replacement
from collections import Counter


def check_compatible_farkas(A_LP, basis):
    """
    Solves the LP to check if 1 can be represented as a non-negative linear 
    combination of the constraints whose coefficients are in A_LP.
    
    The problem is formulated as finding lambda >= 0 such that: 
    A_LP^T * lambda = b_target = [1, 0, 0, ...]
    """
    # The target vector b_target = [1, 0, 0, ...]
    # This means: Coeff of '1' must be 1.0, Coeff of 'x', 'y', etc. must be 0.
    _, n = A_LP.shape
    b_target = np.zeros(n)
    b_target[0] = 1.0 
    return check_compatible(A_LP, basis, b_target)


def check_compatible(A_LP, basis, b_target):
    """
    Solves the LP to check if b_target can be represented as a non-negative linear 
    combination of the constraints whose coefficients are in A_LP.

    Args:
        A_LP (np.ndarray): Coefficient matrix (rows are constraints, columns are monomials).
        basis (list): The ordered basis of monomials.

    Returns:
        tuple: (bool, list): (is_compatible, non_negative_multipliers)
    """
    if A_LP.size == 0:
        return False, []
    
    m_constraints, n_monomials = A_LP.shape

    if n_monomials != len(b_target):
        return False, []

    try:
        model = Model("FarkasCompatibilityLP")
        model.setParam('OutputFlag', 0) 

        # Define the lambda variables (non-negative multipliers)
        lambdas = model.addVars(m_constraints, name="lambda", lb=0.0)

        # Dummy objective (feasibility problem)
        model.setObjective(0, GRB.MINIMIZE)

        # Add constraints: A_LP^T * lambda = b_target
        A_LP_T = A_LP.T 
        
        for j in range(n_monomials):
            # Sum_{i=1}^{m} lambda_i * A_LP_T[j, i] == b_target[j]
            expr = sum(lambdas[i] * A_LP_T[j, i] for i in range(m_constraints))
            
            # Naming the constraints based on the monomial basis
            # name = f"Coeff_{basis[j].subs({s:s.name for s in basis[j].free_symbols})}_is_{b_target[j]}"
            
            model.addConstr(expr == b_target[j])

        model.optimize()

        if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
            is_compatible = True
            multipliers = [lambdas[i].X for i in range(m_constraints)]
            return is_compatible, multipliers
        else:
            return False, []

    except GurobiError as e:
        print(f"Gurobi error: {e}")
        return False, []
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False, []

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

def add_handelman_constraints(model, K, p, variables, max_degree):
    """
    Adds Handelman constraints to the Gurobi model for polynomial p over the 
    constraints K up to the specified maximum degree.

    Args:
        model (gurobipy.Model): The Gurobi optimization model.
        K (list): List of SymPy polynomial expressions representing constraints.
        p (sympy.Expr): The polynomial expression to be constrained.
        variables (list): List of SymPy symbols representing variables.
        max_degree (int): Maximum degree for the Handelman monoid.
    """
    # Generate the Handelman monoid basis
    handelman_basis = get_monoid(K, max_degree)

    # Get the coefficient vector of p with respect to the monomial basis
    monomial_basis = generate_monomials(variables, degree=max_degree)
    p_coeff_vector = get_coefficient_vector(p, variables, monomial_basis)

    # Get the coefficient matrix A_LP for the Handelman basis
    A_LP = get_coefficient_matrix(handelman_basis, variables, monomial_basis)

    # Add constraints: A_LP^T * lambda = p_coeff_vector
    A_LP_T = A_LP.T 

    # Define m_constraints as the number of rows in A_LP
    monoid_rows = A_LP.shape[0]

    lambdas = model.addVars(monoid_rows, name="lambda", lb=0.0)
    for j in range(len(monomial_basis)):
        expr = sum(lambdas[i] * A_LP_T[j, i] for i in range(monoid_rows))
        model.addConstr(expr == p_coeff_vector[j])

    # Check compatibility using Farkas lemma
    
    # compatible, multipliers = check_compatible_farkas(A_LP, monomial_basis)

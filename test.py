from sympy import symbols, Poly, prod, Expr
import numpy as np
from itertools import combinations_with_replacement
from collections import Counter
from gurobipy import Model, GRB, GurobiError




# # --- Verification Example (using the previous setup) ---

# # Define variables (as per previous steps)
# x, y, z = symbols('x, y, z')
# variables = [x, y, z]
# d = 2

# # Monomial Basis (must be generated/ordered consistently)
# monomial_list = generate_monomials(variables, d)

# # Define the polynomial p(x,y,z)
# p = 1 + z + z**2 - 5*x*y + 5*x**2

# # 3. Get the coefficient vector (Now working!)
# c_vector = get_coefficient_vector(p, variables, monomial_list)

# # 4. Print results
# print("--- Debugged Function Output ---")
# print(f"Polynomial P: {p}")
# print(f"Monomial Basis Order (m_i):")
# print(monomial_list)
# print(f"Coefficient Vector (c_i):")
# print(c_vector)

def get_coefficient_matrix(g_constraints, variables, basis):
    """
    Extracts the coefficient matrix A_LP for the Farkas LP problem 
    from a list of linear constraints g_i(x) >= 0.

    Args:
        g_constraints (list): A list of SymPy expressions defining the constraints.
        variables (list) 
        basis (list): A list of monomials

    Returns:
        tuple: (A_LP) where 
               A_LP (np.ndarray): Matrix where A_LP[i, j] is the coefficient 
                                  of the j-th monomial in the i-th constraint g_i.
    """
    if not g_constraints:
        return np.array([]), []
        
    A_LP = []
    for g in g_constraints:
        # Each row is the coefficient vector of one constraint g_i(x)
        row_coeffs = get_coefficient_vector(g, variables, basis)
        A_LP.append(row_coeffs)
    
    return np.array(A_LP), basis

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
    
# 1. Define variables
# x, y = symbols('x, y')
# variables = [x, y]
# d = 2
# basis = generate_monomials(variables, d)
# 2. Define the linear constraints g_i(x) >= 0 (Example from before)
# g1: x - 1 >= 0   
# g2: y >= 0   
# g3: 3 - x - y >= 0 
# g_constraints = [x - 1, y, 3 - x - y]

# # 3. Step 1: Extract the Coefficient Matrix
# A_LP, basis = get_coefficient_matrix(g_constraints, variables, basis)

# print("--- Data Extraction ---")
# print(f"Monomial Basis: {basis}")
# print("Coefficient Matrix A_LP (Rows=g_i, Cols=Monomials):")
# print(A_LP)

# # 4. Step 2: Check Compatibility
# compatible, multipliers = check_compatible_farkas(A_LP, basis)

# print("\n--- Compatibility Check ---")
# print(f"System compatible (1 is a non-negative linear combination of constraints)? {compatible}")
# if compatible:
#     print(f"Found non-negative multipliers lambda_i: {multipliers}")
#     print("This confirms the existence of a positive linear combination equal to 1.")

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


# 1. Define variables and constraints
x, y = symbols('x, y')
g1 = x**2 - 1      # Constraint g1 >= 0
g2 = 2 - y         # Constraint g2 >= 0
g_constraints = [g1, g2]
d = 2

# 2. Generate the monoid elements
monoid_list = get_monoid(g_constraints, d)

# 3. Print the result
print(f"Constraints: (g1, g2) = ({g1}, {g2})")
print(f"Maximum Monoid Degree: d={d}")
print("Handelman Monoid Elements H_d:")
print(monoid_list)

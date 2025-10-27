import numpy as np
from itertools import combinations_with_replacement
from collections import Counter
from sympy import symbols, Poly, expand, Add, Mul
from gurobipy import Model


def purify(poly, variables):
        """
        Removes all monomials in the polynomial that contain any of the specified variables.

        Args:
            poly (sympy.Expr): The polynomial expression to be purified.
            variables (list): A list of sympy.Symbol objects representing the variables whose monomials should be removed.

        Returns:
            sympy.Expr: The modified polynomial with specified monomials removed.
        """

        # Filter out terms containing any of the variables
        filtered_terms = [term for term in poly.as_ordered_terms() if not any(var in term.free_symbols for var in variables)]
        
        # Reconstruct the polynomial
        modified_poly = Add(*filtered_terms)
        return modified_poly

def get_coeff_vector_unknown(p, variables, monomial_list):
    """
    Computes the coefficient vector of a polynomial with respect to a given list of monomials.
    Args:
        p (sympy.Expr): The polynomial expression to extract coefficients from.
        variables (list): A list of variables present in the polynomial.
        monomial_list (list): A list of monomials to match against the polynomial.
    Returns:
        list: A list of coefficients corresponding to the monomials in `monomial_list`.
              Each coefficient is purified with respect to the given variables.
    """
    p = expand(p)
    coeff_vector = []
    
    for monomial in monomial_list:
        coeff = p.coeff(monomial)
        coeff_vector.append(purify(coeff, variables)) 
        
    return coeff_vector


def generate_monomials(variables, degree):
    """
    Generates a list of monomials for a given list of variables with total degree <= d.

    Args:
        variables (list): A list of sympy.Symbol objects (e.g., [x, y, z]).
        degree (int): The maximum total degree.

    Returns:
        list: A list of sympy expressions representing the monomials.
    """
    num_vars = len(variables)
    monomials = []
    
    # Iterate through all possible total degrees k from 0 to d
    for k in range(degree + 1):
        for combination in combinations_with_replacement(variables, k):
            # The Counter counts the frequency of each variable in the combination.
            # E.g., for (x, x, y), the counter is {x: 2, y: 1}
            exponent_map = Counter(combination)
            
            # Start the monomial as 1 (for the degree 0 case)
            monomial = 1
            
            # Build the monomial expression: x1**e1 * x2**e2 * ...
            for var in variables:
                # Get the exponent for the current variable. Defaults to 0 if not present.
                exponent = exponent_map.get(var, 0)
                monomial *= var**exponent
            
            monomials.append(monomial)
            
    return monomials

def get_coeff_vector(p, variables, monomial_list):
    """
    Obtains the coefficient vector of a polynomial p according to a specific 
    list of monomials.

    Args:
        p (sympy.Expr): The polynomial expression.
        variables (list): The list of sympy.Symbol objects (e.g., [x, y, z]).
        d (int): The maximum total degree of the monomial basis.
        monomial_list (list): The ordered list of monomials (m_i) forming the basis.

    Returns:
        list: The coefficient vector (c_i) such that p = sum(c_i * m_i).
    """

    p = expand(p)
    coeff_vector = []
    
    for monomial in monomial_list:
        coeff = p.coeff(monomial)
        # coeff = p.coeff_monomial(monomial)
        coeff_vector.append(coeff) 
        
    return coeff_vector

def get_coeff_matrix(g_constraints, variables, basis):
    """
    Extracts the coefficient matrix A for the LP problem 
    from a list of linear constraints g_i(x) >= 0.

    Args:
        g_constraints (list): A list of SymPy expressions defining the constraints.
        variables (list) 
        basis (list): A list of monomials

    Returns:
        A (np.ndarray): Matrix where A[i, j] is the coefficient 
            of the j-th monomial in the i-th constraint g_i.
    """
    if not g_constraints:
        return np.array([]), []
        
    A = []
    for g in g_constraints:
        # Each row is the coefficient vector of one constraint g_i(x)
        row_coeffs = get_coeff_vector(g, variables, basis)
        A.append(row_coeffs)
    
    return np.array(A)
from itertools import combinations_with_replacement
from collections import Counter
from sympy import symbols, Add, Mul, simplify
from typing import List, Any
import numpy as np

def generate_monomials(variables: List[Any], degree: int) -> List[Any]:
    """
    Generates a list of monomials for a given list of variables with total degree <= d.

    Args:
        variables (list): A list of sympy.Symbol objects (e.g., [x, y, z]).
        degree (int): The maximum total degree.

    Returns:
        list: A list of sympy expressions representing the monomials.
    """
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

def template_poly(variables, degree, name: str):
    """
    Creates a template polynomial with symbolic coefficients for given variables and degree.

    Args:
        variables (list): List of sympy.Symbol objects.
        degree (int): Maximum total degree of the polynomial.

    Returns:
        sympy.Expr: The template polynomial expression.
    """

    monomials = generate_monomials(variables, degree)
    coeffs = symbols(f'{name}c0:{len(monomials)}')

    poly = Add(*[Mul(coeffs[i], monomials[i]) for i in range(len(monomials))])
    return poly, coeffs, monomials

def degree_in_variables(p, variables: List[Any]) -> int:
    """
    Computes the total degree of polynomial p with respect to the given variables.
    """
    max_degree = 0
    for monomial in p.as_ordered_terms():
        total_deg = sum(monomial.as_powers_dict().get(var, 0) for var in variables)
        if total_deg > max_degree:
            max_degree = total_deg
    return max_degree

def homogenization(p, variables: List[Any], h_var: Any) -> Any:
    """
    Homogenizes a polynomial p by adding a new variable h_var.

    Args:
        p: The polynomial to homogenize.
        variables (list): List of sympy.Symbol objects representing the variables in p.
        degree (int): The target total degree after homogenization.
        h_var: The new variable to add for homogenization.

    Returns:
        The homogenized polynomial.
    """
    degree = degree_in_variables(p, variables)
    # print("degree", degree)
    monomial_list = generate_monomials(variables, degree)
    
    coeffs = get_coeff_vector(p, variables, monomial_list)
    h_p = 0
    for c, m in zip(coeffs, monomial_list):
        if m == 1:
            total_deg = 0
        else:
            total_deg = sum(m.as_powers_dict().values())

        h_p += c * m * (h_var ** (degree - total_deg))
    
    return simplify(h_p)
    
    
    
    # homog_p = 0

    # for i, monomial in enumerate(monomial_list):
    #     if i==1:
    #         coeff = p.as_coeff_add()[0]
    #         total_deg = 0
    #     else:
    #         coeff = (p.coeff(monomial)).as_coeff_add()[0]
    #         total_deg = sum(monomial.as_powers_dict().values())
    #     h_exponent = degree - total_deg
    #     if h_exponent == 0:
    #         homog_monomial = monomial
    #     else:
    #         homog_monomial = monomial * (h_var ** h_exponent)
    #     print(coeff, monomial, total_deg, homog_monomial)
    #     homog_p += coeff * homog_monomial

    # return homog_p

def get_coeff_vector(p, variables: List[Any], monomial_list: List[Any]) -> List[Any]:
    """
    Obtains the coefficient vector of a polynomial p according to a specific 
    list of monomials.

    Example:
        input: 3*x**2 + 2*x*y + y + 5, [x,y], [1, x, y, x**2, x*y, y**2]
        output: [5, 0, 1, 3, 2, 0]
    """
    # Initialize the coefficient vector
    coeff_vector = []
    p = p.expand()

    # Iterate through the monomial list
    for monomial in monomial_list:
        # Extract the coefficient of the current monomial in the polynomial
        coeff = p.coeff(monomial)
        constant = coeff.as_independent(*variables, as_Add=True)[0]
        coeff_vector.append(constant)
    
    coeff_vector[0] = p.as_independent(*variables, as_Add=True)[0]
    return coeff_vector

def get_coeff_matrix(g_constraints: List[Any], variables: List[Any], monomial_list: List[Any]) -> List[List[Any]]:
    """
    Extracts the coefficient matrix A for the LP problem 
    from a list of linear constraints g_i(x) >= 0.
    """

    A = []
    for g in g_constraints:
        # Each row is the coefficient vector of one constraint g_i(x)
        row_coeffs = get_coeff_vector(g, variables, monomial_list)
        A.append(row_coeffs)
    
    return np.array(A)


if __name__ == "__main__":
    # tests
    # print("\nTest: get_coeff_vector")
    x, y = symbols('x, y') 
    variables = [x, y]
    monomial_list = generate_monomials(variables, degree=2)
    # print(monomial_list)
    # p1 = 3*x**2 + 2*x*y + y + 5
    # print("Polynomial p:", p1)
    # print("monomial_list:", monomial_list)
    # cv = get_coeff_vector(p1, variables, monomial_list)
    # print(cv, "==> Expect: [5, 0, 1, 3, 2, 0]")

    # print("\nTest: get_coeff_vector with parameters")
    # x, y = symbols('x, y') 
    # variables = [x, y]
    a, b = symbols('a, b')
    # monomial_list = generate_monomials(variables, degree=2)
    # print(monomial_list)
    # p2 = a*x**2 + b*x*y + y + 5
    # print("Polynomial p:", p2)
    # print("monomial_list:", monomial_list)
    # cv = get_coeff_vector(p2, variables, monomial_list)
    # print(cv, "==> Expect: [5, 0, 1, a, b, 0]")

    # print("\nTest: get_coeff_matrix")
    # cv = get_coeff_matrix([p1,p2], variables, monomial_list)
    # print(cv)

    print("\nTest: homogenization")
    p1 = x**2 + 2*x*y + 3*y + 4
    p2 = x**2 + a*x*x*y + b**2*y + 4
    x0 = symbols('x0')

    homog_p1 = homogenization(p1, variables, h_var=x0)
    print("Original polynomial p1:", p1)
    print("Homogenized polynomial homog_p1:", homog_p1)

    homog_p2 = homogenization(p2, variables, h_var=x0)
    print("Original polynomial p2:", p2)
    print("Homogenized polynomial homog_p2:", homog_p2)

    v, v_coeffs, monomial_list = template_poly(variables, degree=2, name='v')
    print("Template polynomial v:", v)
    print("homog v:", homogenization(v, variables, h_var=x0))
    
    # print("\nTest: template_poly", v)
    # print("v_coeffs:", v_coeffs)

    # monomial_list = generate_monomials(variables, degree=4)
    # print("v**2", v**2)
    # c = get_coeff_vector(v**2, variables, monomial_list)
    # print("Coefficient vector of v:", c)
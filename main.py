from util.basic import generate_monomials, get_coeff_vector, get_coeff_vector_unknown, get_coeff_matrix
from util.handelman import check_compatible_farkas, check_compatible, add_handelman_constraints
from util.entailment import Entailment
from sympy import symbols, simplify, Add, Mul, expand, Abs
from gurobipy import Model, GRB, GurobiError



def template_poly(variables, degree):
    """
    Creates a template polynomial with symbolic coefficients for given variables and degree.

    Args:
        variables (list): List of sympy.Symbol objects.
        degree (int): Maximum total degree of the polynomial.

    Returns:
        sympy.Expr: The template polynomial expression.
    """

    monomials = generate_monomials(variables, degree)
    coeffs = symbols(f'c0:{len(monomials)}')

    poly = Add(*[Mul(coeffs[i], monomials[i]) for i in range(len(monomials))])
    return poly, coeffs




def main():
    
    # 1. Define problem
    x, y = symbols('x, y') 
    variables = [x, y]


    # 2. Create templates
    p, p_coeffs = template_poly(variables, degree=2) 

    print("Polynomial p:", p)
    print("Polynomial coefficients:", p_coeffs)

    monomial_list = generate_monomials(variables, degree=2)
    cv = get_coeff_vector_unknown(p, variables, monomial_list)
    print(cv)

    # 3. Declare constraints in SymPy
    # cons = Entailment([x,1-x], p)

    # print(cons.premises)
    # print(cons.conclusion)

    print(remove_abolutee_values(2*abs(x-1) + 3*abs(x) + 4))


    
    # 4. Translate the constraints to Gurobi

    # model = Model("Handelman")
    # coeff_gurobi_vars = []
    # for i, coeff in enumerate(cv):
    #     coeff_gurobi_vars.append(model.addVar(name=f'c{i}', lb=-GRB.INFINITY, ub=GRB.INFINITY))
        
    # # for i in range(len(cv)):
    #     # cv[i].subs(dict(zip(p_coeffs, coeff_gurobi_vars)))


    # print(cv[0])
    # s = str(cv[0])
    # for i in range(len(p_coeffs)):
    #     s = s.replace(str(p_coeffs[i]), f'coeff_gurobi_vars[{i}]')
    # print("Substituted cv[0]:", eval(s))

    # q = p.subs(x, 2*x + y)
    # print("q", expand(q))
    # c = get_coeff_vector_unknown(q, variables, monomial_list)
    # print(c)

    # 1. Define problem
    # x, y = symbols('x, y')
    # variables = [x, y]
    # K1 = [x, 1 - x, y, 1 - y] 
    # K2 = [x-2, 3-x, y-2, 3-y]

    # p = template_poly(variables, degree=2) 
    # 

    # monomial_list = generate_monomials(variables, degree=2)
    # model = Model("Handelman")
    # # Declare n variables in the LP model
    # print(p.coeff(x))
    # print(get_coeff_vector(p,variables,monomial_list))
    
    
    # print(q.coeff(x))

    # v = get_coefficient_vector_unknown(model, p, variables, monomial_list)


    # model.setParam('OutputFlag', 0)
    # model.setObjective(0, GRB.MINIMIZE)

    # add_handelman_constraints(model, K1, p, variables, max_degree=2)
    # add_handelman_constraints(model, K2, -p, variables, max_degree=2)

    # model.optimize()
    # if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
    #     is_compatible = True
    #     multipliers = [lambdas[i].X for i in range(m_constraints)]
    #     return is_compatible, multipliers
    # else:
    #     return False, []

    # basis = generate_monomials(variables=[x,y], d=1)
    # # 2. Data Preparation
    # A_LP, basis = get_coefficient_matrix(K, [x,y], basis)

    # # 3. Solve
    # compatible, multipliers = check_compatible_farkas(A_LP, basis)

    # # 4. Report
    # print("--- Handelman/Farkas Compatibility Check ---")
    # print(f"Constraints: {K}")
    # print(f"Monomial Basis: {basis}")
    # print(f"Result: Compatible? {compatible}")
    # if compatible:
    #     print(f"Multipliers found: {multipliers}")

if __name__ == "__main__":
    main()
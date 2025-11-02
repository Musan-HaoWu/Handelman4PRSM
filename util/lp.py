from sympy import symbols
from typing import List, Any
from util.basic import *
from gurobipy import Model, GRB



# def prepare_gurobi_variables(model, sympy_coeffs: List[List[Any]]) -> List[List[Any]]:
#     v_coeffs_gvars = []
#     for i, coeff in enumerate(sympy_coeffs[0]):
#         v_coeffs_gvars.append(model.addVar(name=f'vc{i}', lb=-GRB.INFINITY, ub=GRB.INFINITY))
    
#     p_coeffs_gvars = []
#     for i, coeff in enumerate(sympy_coeffs[1]):
#         p_coeffs_gvars.append(model.addVar(name=f'pc{i}', lb=-GRB.INFINITY, ub=GRB.INFINITY))
    
#     d_coeffs_gvars = []
#     for i, coeff in enumerate(sympy_coeffs[2]):
#         d_coeffs_gvars.append(model.addVar(name=f'dc{i}', lb=-GRB.INFINITY, ub=GRB.INFINITY))

#     gurobi_coeffs = [v_coeffs_gvars, p_coeffs_gvars, d_coeffs_gvars]
#     model.update() 
#     return gurobi_coeffs

def sympy2gurobi(sympy_expr, sympy_coeffs: List[List[Any]], gurobi_coeffs: List[List[Any]]) -> Any:
    
    sympy_coeffs = [item for sublist in sympy_coeffs for item in sublist]
    gurobi_coeffs = [item for sublist in gurobi_coeffs for item in sublist]
    expr_str = str(sympy_expr)
    for i in range(len(sympy_coeffs)):
        expr_str = expr_str.replace(str(sympy_coeffs[i]), f"gurobi_coeffs[{i}]")
    expr = eval(expr_str)
    return expr


if __name__ == "__main__":
    # tests
    x, y = symbols('x, y') 
    variables = [x, y]
    
    model = Model("Test")
    v, v_coeffs = template_poly(variables, degree=2, name='v')
    p, p_coeffs = template_poly(variables, degree=1, name='p')
    d, d_coeffs = template_poly(variables, degree=3, name='d')

    sympy_coeffs = [v_coeffs, p_coeffs, d_coeffs]
    print("SymPy coefficient symbols:", sympy_coeffs)
    gurobi_coeffs = prepare_gurobi_variables(model, sympy_coeffs)

    # A = [v, p, v*p]

    monomial_list = generate_monomials(variables, degree=2)
    expr = get_coeff_vector(v+2, variables, monomial_list)[0]
    print("Sympy Expr:", expr, " of type: ", type(expr))

    # print(gurobi_coeffs)
    expr_g = sympy2gurobi(expr, sympy_coeffs, gurobi_coeffs)
    # print(type(g_c))
    print("Gurobi Expr:", expr_g, " of type: ", type(expr_g))
    
    # uu = model.addVar(name='u', lb=0.0, ub=GRB.INFINITY)
    # model.update()
    # print(uu>=0)
    
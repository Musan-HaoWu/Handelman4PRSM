from util.basic import generate_monomials, template_poly, get_coeff_vector, get_coeff_matrix, homogenization
from util.entailment import Entailment, Constrainst

from util.handelman import check_compatible_farkas, check_compatible, add_handelman_constraints
from sympy import symbols, simplify, Add, Mul, expand, Abs
from gurobipy import Model, GRB, GurobiError


def main():
    
    # 1. Define problem
    x = symbols('x') 
    variables = [x]

    x0 = symbols('x0') 
    
    # 2. Create templates
    v, v_coeffs = template_poly(variables, degree=1) 
    expect_v = x/(2*x+1)*v.subs({x: x - 1}) + (x+1)/(2*x+1)*v.subs({x: x + 1})
    # rational case
    v = v.subs({x: Abs(x)/(Abs(x)+1)})

    print("Template v:", v) 
    print("Expect v", expect_v)

    # 3. Declare constraints in SymPy
    prsm_constraints = Constrainst()
    prsm_constraints.add_constraint(Entailment([x, -x], -v))
    prsm_constraints.add_constraint(Entailment([x-1], v-1))
    prsm_constraints.add_constraint(Entailment([x-1], v-expect_v))
    
    prsm_constraints.print_constraints()
    prsm_constraints.rewrite()
    prsm_constraints.print_constraints()

    # 4. Translate the constraints to Gurobi

    D = 4

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
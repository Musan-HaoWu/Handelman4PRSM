from util.basic import generate_monomials, template_poly, get_coeff_vector, get_coeff_matrix, homogenization
from util.entailment import Entailment, Constrainst
from sympy import symbols, simplify, Add, Mul, expand, Abs
from gurobipy import Model, GRB, GurobiError
# from lp import prepare_gurobi_variables, sympy2gurobi

def main():
    
    # 1. Define problem
    x = symbols('x') 
    vars = [x]

    x0 = symbols('x0') 
    
    # 2. Create templates
    v, v_coeffs, monomials = template_poly(vars, degree=1, name='v') 
    # sympy_coeffs = [v_coeffs, [], []]
    # gurobi_coeffs = prepare_gurobi_variables(model, sympy_coeffs)
   
    
    # rational case
    if True:
        v = v.subs({x: x/(Abs(x)+1)})
    expect_v = x/(2*x+1)*v.subs({x: x - 1}) + (x+1)/(2*x+1)*v.subs({x: x + 1})
    

    print("Template v:", v) 
    print("Expect v", expect_v)

    # 3. Declare constraints in SymPy
    v_constraints = Constrainst(vars, v_coeffs)
    v_constraints.add_constraint(Entailment(vars, [x, -x], -v))
    v_constraints.add_constraint(Entailment(vars, [x-1], v-1))
    v_constraints.add_constraint(Entailment(vars, [x-1], v-expect_v))
    # v_constraints.add_constraint(Entailment(vars, [x-1, 2-x, x-3], v)) #redundant
    
    print(simplify(v - expect_v))

    # 4. Rewrite constraints into standard form (x in K => f(x) >= 0)
    v_constraints.print_constraints()
    v_constraints.rewrite()
    v_constraints.print_constraints()

    # 5. Remove redundant constraints by using Farkas' lemma
    v_constraints.remove_empty_K()
    print("remove empty K")
    v_constraints.print_constraints()

    # 6. Formulate D-th degree Handelman relaxation
    D = 6
    v_constraints.translate_handelman(D)

    # 7. Solve LP for v
    result, v_coeffs_sol = v_constraints.solve_for_v()
    if result:
        print("Solution: v = ", sum(x * y for x, y in zip(v_coeffs_sol, monomials)))
    else:
        print("No solution found.")

if __name__ == "__main__":
    main()
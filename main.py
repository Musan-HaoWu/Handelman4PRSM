from util.basic import generate_monomials, template_poly, get_coeff_vector, get_coeff_matrix, homogenization
from util.entailment import Entailment, Constrainst
from sympy import symbols, simplify, Add, Mul, expand, Abs
from gurobipy import Model, GRB, GurobiError
# from lp import prepare_gurobi_variables, sympy2gurobi

def main():

    # 1. Define problem
    x = symbols('x') 
    vars = [x]
    
    # # 1d_rw (v, p, d)
    # prob = [1/2, 1/2]
    # next_vars = [ x - 1, x + 1 ]
    # deg_v = 1
    # is_rational_v = False
    # deg_p = 1
    # deg_d = 1

    # # fair in the limit (v is ok, but cannot find p and d) 
    # prob = [x/(2*x+1), (x+1)/(2*x+1)]
    # next_vars = [ x - 1, x + 1 ]
    # deg_v = 1
    # is_rational_v = True
    # deg_p = 1
    # deg_d = 1

    # escaping (v, p, d)
    prob = [1/(x+1), x/(x+1)]
    next_vars = [ 0, x + 1 ]
    deg_v = 1
    is_rational_v = False
    deg_p = 1
    deg_d = 0



    # 2. Create templates
    v_poly, v_coeffs, monomials = template_poly(vars, degree=deg_v, name='v') 
    
    if is_rational_v:
        v = v_poly.subs({x: x/(Abs(x)+1)})
    else:
        v = v_poly
    
    expect_v = 0
    for p, nv in zip(prob, next_vars):
        expect_v += p * v.subs({x: nv})

    print("Template v:", v) 
    print("Expect v", expect_v)

    # 3. Declare constraints in SymPy
    v_constraints = Constrainst(v_coeffs)
    v_constraints.add_constraint(Entailment(vars, [x, -x], -v, is_bounded = True))  # v(x) >= 0 for all x
    v_constraints.add_constraint(Entailment(vars, [x-1], v-1))
    v_constraints.add_constraint(Entailment(vars, [x-1], v-expect_v))
    print("\n=== Initial Constraints for v ===")
    v_constraints.print_constraints()

    # 4. Rewrite and Prune
    v_constraints.rewrite()

    # 5. Solve for increasing degree D Handelman Relaxation
    min_D = v_constraints.max_degree_f()
    result = False
    v_sol = None
    for D in range(min_D, 10):
        print(f"\n=== Degree {D} Handelman Relaxation ===")    
        v_constraints.translate_handelman(D)
        result, v_coeffs_sol = v_constraints.solve()
        if result:
            v_sol = sum(x * y for x, y in zip(v_coeffs_sol, monomials))
            v_new = v_sol
            if is_rational_v:
                v_new = v_sol.subs({x: x/(Abs(x)+1)})
            print("Solution: v = ", v_new)
            break
        else:
            print("No solution found.")
            v_constraints.refresh_model()
        
    if not result:
        return
    
    # 6. Solving constraints for p and d, w.r.t. the found v
    print("\n=== Now solving for p and d ===")
    r = symbols('r')
    p_poly, p_coeffs, p_monomials = template_poly([r], degree=deg_p, name='p')
    d_poly, d_coeffs, d_monomials = template_poly([r], degree=deg_d, name='d')
    p = p_poly.subs({r: r/(Abs(r)+1)})
    d = d_poly.subs({r: r/(Abs(r)+1)})
    pd_constraints = Constrainst(p_coeffs + d_coeffs)
    pd_constraints.add_constraint(Entailment([r], [r, 1-r], -d_poly.diff(r), is_bounded=True))
    pd_constraints.add_constraint(Entailment([r], [r, -r], d_poly-0.01, is_bounded=True)) # d(r) > 0
    pd_constraints.add_constraint(Entailment([r], [r, 1-r], d_poly, is_bounded=True)) # d(r) > 0
    pd_constraints.add_constraint(Entailment([r], [r, 1-r], -p_poly.diff(r), is_bounded=True))
    pd_constraints.add_constraint(Entailment([r], [r, 1-r], p_poly, is_bounded=True)) # p(r) > v(r)
    pd_constraints.add_constraint(Entailment([r], [r, -r], p_poly-0.01, is_bounded=True)) # p(r) > v(r)
    pd_constraints.add_constraint(Entailment([r], [r, 1-r], 1 - p_poly, is_bounded=True)) # p(r) < 1
    
    pd_constraints.add_constraint(Entailment(vars, [x-1], v_new - v_new.subs({x: next_vars[0]}) - d.subs({r: v_new})))
    pd_constraints.add_constraint(Entailment(vars, [x-1], prob[0] - p.subs({r: v_new})))

    print("\n=== Initial Constraints for p and d ===")
    pd_constraints.print_constraints()
    pd_constraints.rewrite()
    min_D = pd_constraints.max_degree_f()
    result = False
    print("min_D: ",min_D)
    for D in range(min_D, 8):
        print(f"\n=== Degree {D} Handelman Relaxation for p,d ===")    
        pd_constraints.translate_handelman(D)
        result, pd_coeffs_sol = pd_constraints.solve()
        if result:
            p_sol = sum(x * y for x, y in zip(pd_coeffs_sol[:len(p_coeffs)], p_monomials))
            d_sol = sum(x * y for x, y in zip(pd_coeffs_sol[len(p_coeffs):], d_monomials))
            print("Solution: v = ", v_sol)
            print("Solution: p = ", p_sol)
            print("Solution: d = ", d_sol)
            break
        else:
            print("No solution found.")
            pd_constraints.refresh_model()

if __name__ == "__main__":
    main()
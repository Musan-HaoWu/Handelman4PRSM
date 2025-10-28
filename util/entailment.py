from sympy import symbols, simplify, Abs
from itertools import product

class Entailment:
    def __init__(self, gs, f):
        self.K = gs
        self.f = f

    def remove_absolute_values(self):
        """
        Takes an expression with possible absolute values 
        and removes them by splitting the cases.
        For example:
            Input: a*Abs(x-1) + b*Abs(x) + c, 
            Output: 
                [
                ([x-1, x], a*(x-1) + b*x + c),
                ([1-x, x], a*(1-x) + b*x + c),
                ([x-1, -x], a*(x-1) - b*x + c),
                ([1-x, -x], a*(1-x) - b*x + c)
                ]
            Nested example:
            Input: Abs(x-a*Abs(y))
            Output:
                [
                ([y, x-a*y], x-a*y),
                ([-y, x+a*y], x+a*y),
                ([y, a*y-x], a*y-x),
                ([-y, -x-a*y], -x-a*y)
                ]
        """
        cases = []

        # Find all absolute value terms in the expression
        # Note: the order f.atoms() returns is random!
        abs_terms = list(sorted(self.f.atoms(Abs), key=str, reverse=True) )  
        
        if not abs_terms:
            return [self]  # No absolute values, return the current entailment as is

        # Process the first absolute value term
        abs_term = abs_terms[0]
        inner_expr = abs_term.args[0]  # Extract the inner expression of Abs

        # Create two cases: one for the positive and one for the negative of the inner expression
        positive_case = self.f.subs(abs_term, inner_expr)
        negative_case = self.f.subs(abs_term, -inner_expr)

        # Create new entailments for each case
        cases.append(Entailment(self.K + [inner_expr], positive_case))
        cases.append(Entailment(self.K + [-inner_expr], negative_case))

        # Recursively handle remaining absolute values in each case
        result = []
        for case in cases:
            result.extend(case.remove_absolute_values())

        return result

    def remove_denominator(self):
        """
        If f is a rational expression, then remove its denominator.
        For example:
            Input: Abs(x)/(1+Abs(x))
            Output: Abs(x)

            Input: a*Abs(x)+1/(Abs(x)+1) * x/(2*x+1), 
            Output: a*Abs(x)*x+x
        """

        expr = self.f
        # Find all terms that are denominators in the expression
        numerator, denominator = expr.as_numer_denom()
        # Simplify the resulting expression
        expr = simplify(numerator)

        # Update the entailment's function
        self.f = expr

    def remove_duplicate(self):
        """
        Remove duplicate constraints in K
        """
        unique_constraints = list(set(self.K))
        self.K = unique_constraints

class Constrainst:
    def __init__(self):
        self.constraints = []

    def add_constraint(self, entailment):
        if isinstance(entailment, Entailment):
           self.constraints.append(entailment)
        else:
            raise TypeError("Only Entailment instances can be added as constraints.")

    def remove_absolute_values_all(self):
        new_constraints = []
        for entailment in self.constraints:
            cases = entailment.remove_absolute_values()
            new_constraints.extend(cases)
        self.constraints = new_constraints

    def remove_denominators_all(self):
        for entailment in self.constraints:
            entailment.remove_denominator()

    def print_constraints(self):
        print("Number of constraints:", len(self.constraints))
        for i, constraint in enumerate(self.constraints):
            print(f"Constraint {i+1}: ", constraint.K, ">= 0 ==>", simplify(constraint.f), ">=0")

    def remove_duplicates_all(self):
        for entailment in self.constraints:
            entailment.remove_duplicate()
    
    def rewrite(self):
        self.remove_denominators_all()
        self.remove_absolute_values_all()
        self.remove_duplicates_all()

if __name__ == "__main__":
    # tests
    x, y = symbols('x, y') 
    variables = [x, y]
    a = symbols('a')

    constraints = Constrainst()

    K1 = [x+10, 10-x] 
    f1 = 2*x
    entailment1 = Entailment(K1, f1)
    constraints.add_constraint(entailment1)

    K2 = [x-2, 3-x, y-2, 3-y]
    f2 = Abs(x - Abs(1-y)) + 1
    entailment2 = Entailment(K2, f2)
    constraints.add_constraint(entailment2)

    K3 = [x+10, 10-x]
    f3 = 2*Abs(3*x-1)/(1+Abs(x))*x/(2*x+1)
    entailment3 = Entailment(K3, f3)
    constraints.add_constraint(entailment3)

    K4 = [x+10, 10-x]
    f4 = y+1/(1+Abs(a/(1+Abs(x))))
    entailment4 = Entailment(K4, f4)
    constraints.add_constraint(entailment4)

    constraints.print_constraints()

    print("After removing denominators")
    constraints.remove_denominators_all()
    constraints.print_constraints()

    print("After removing absolute values:")
    constraints.remove_absolute_values_all()
    constraints.print_constraints()


from gurobipy import Model, GRB
from sympy import symbols, expand

class GurobiSympyPoly:
    def __init__(self, coeffs, sympy_vars):
        """
        Initialize the polynomial expression.
        :param coeffs: List of Gurobi variables [c0, c1, ..., cn].
        :param sympy_vars: List of SymPy variables [x, y, ...].
        """
        self.coeffs = coeffs
        self.sympy_vars = sympy_vars
        self.poly = sum(c * monomial for c, monomial in zip(coeffs, self._generate_monomials()))

    def _generate_monomials(self):
        """
        Generate monomials for the polynomial based on the SymPy variables.
        For example, for x, y: [1, x, y, x**2, x*y, y**2].
        """
        monomials = [1]  # Constant term
        monomials.extend(self.sympy_vars)  # Linear terms
        for i in range(len(self.sympy_vars)):
            for j in range(i, len(self.sympy_vars)):
                monomials.append(self.sympy_vars[i] * self.sympy_vars[j])  # Quadratic terms
        return monomials

    def __add__(self, other):
        """
        Overload addition for two GurobiSympyPoly objects.
        """
        if not isinstance(other, GurobiSympyPoly):
            raise TypeError("Can only add GurobiSympyPoly objects.")
        new_coeffs = [self.coeffs[i] + other.coeffs[i] for i in range(len(self.coeffs))]
        return GurobiSympyPoly(new_coeffs, self.sympy_vars)

    def __mul__(self, scalar):
        """
        Overload multiplication by a scalar.
        """
        new_coeffs = [scalar * c for c in self.coeffs]
        return GurobiSympyPoly(new_coeffs, self.sympy_vars)

    def extract_coefficients(self):
        """
        Extract the coefficient vector (list of Gurobi variables).
        """
        return self.coeffs

    def __repr__(self):
        """
        String representation of the polynomial.
        """
        return str(expand(self.poly))


# Example usage
if __name__ == "__main__":
    # Create a Gurobi model
    model = Model("poly_example")

    # Define Gurobi variables
    c0 = model.addVar(vtype=GRB.CONTINUOUS, name="c0")
    c1 = model.addVar(vtype=GRB.CONTINUOUS, name="c1")
    c2 = model.addVar(vtype=GRB.CONTINUOUS, name="c2")
    c3 = model.addVar(vtype=GRB.CONTINUOUS, name="c3")
    c4 = model.addVar(vtype=GRB.CONTINUOUS, name="c4")
    c5 = model.addVar(vtype=GRB.CONTINUOUS, name="c5")

    # Define SymPy variables
    x, y = symbols("x y")

    # Create a polynomial
    poly1 = GurobiSympyPoly([c0, c1, c2, c3, c4, c5], [x, y])
    print("Polynomial 1:", poly1)

    # Create another polynomial
    poly2 = GurobiSympyPoly([c0, c1, c2, c3, c4, c5], [x, y])
    print("Polynomial 2:", poly2)

    # Add two polynomials
    poly_sum = poly1 + poly2
    print("Sum of Polynomials:", poly_sum)

    # Multiply polynomial by a scalar
    poly_scaled = poly1 * 2
    print("Scaled Polynomial:", poly_scaled)

    # Extract coefficients
    coeffs = poly1.extract_coefficients()
    print("Coefficients:", coeffs)
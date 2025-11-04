from sympy import symbols

class Program:
    def __init__(self):
        pass



name = "Fair in the limit"

x = symbols('x')
vars = [x]

deg_v = 1
is_rational_v = True

deg_p = 1
deg_d = 2

prob = [x/(2*x+1), (x+1)/(2*x+1)]
next_vars = [ x - 1, x + 1 ]

print("include")
-- /* list3_5.lua */
   
require"Complex"; require"newton_1"
-- Cube root of 1X10^20
function f1(x) return x^3 - 1.E20 end 
print('f1',newton(f1,1.e10))
-- Real and complex roots of 4th degree polynomial   
function f2(x) return 5*x^4 - 4*x^2 + 2*x -3 end
print('f2',newton(f2,2))
print('f2',newton(f2,j))
print('f2',newton(f2,-j))
-- Function that varies rapidly with x
function f3(x) return math.exp(10*x)*(1-.5*x) -1 end
print('f3',newton(f3,4))
-- Rearranged f3 for more linear function
function f4(x) return (1-.5*x) - math.exp(-10*x) end
print('f4',newton(f4,4))
-- Linear equation
function f5(x) return 5*x - math.pi end
print('f5',newton(f5,0))
-- Complex roots of 2nd degree polynomial
function f6(x) return x^2 -4*x + 13 end
print('f6',newton(f6,j))
print('f6',newton(f6,-j))

-- /* File list4_20.lua */
-- Examples of Roots of Polynomials
require "Polynomial"
getfenv(nsolv).NMAX=200

p0 = Polynomial.new{1,1,4} -- 1 + x + 4x^2
rts = p0:roots() -- Or Polynomial.roots(p0)
print("Roots of p0 polynomial\n",p0)
table.foreach(rts,print)

p1 = Polynomial.new{1,2,3,4,5} -- 1 + 2x + 3x^2 + 4x^3 + 5x^4
rts = Polynomial.roots(p1)
print("\nRoots of p1 polynomial\n",p1)
table.foreach(rts,print)

p2 = Polynomial.new{-1,1}
p2 = p2*{-2,1}*{-3,1}*{-4,1} -- roots of 1,2,3,4
rts = p2:roots()
print("\nRoots of p2 polynomial\n",p2)
table.foreach(rts,print)

p3 = p1*p2
rts = Polynomial.roots(p3)
print("\nRoots of p1*p2 polynomial\n",p3)
table.foreach(rts,print)

p3 = Polynomial.new{-200,0,0,0,0,0,0,0,0,0,0,0,1} -- x^12 - 200 = 0
rts = p3:roots() -- 12 roots of x^12 = 2
print("\nRoots of p3 polynomial\n",p3)
table.foreach(rts,print)


p2 = Polynomial.scale(p1,100000)
rts = Polynomial.roots(p2)
print("\nPolynomial p1 scaled by factor of 1e5\n",p2)
table.foreach(rts,print)

p3 = p1:scale(0.001)
rts = p3:roots()
print("\nPolynomial p1 scaled by factor of 1e-3\n",p3)
table.foreach(rts,print)

p1 = p2*p3
rts = Polynomial.roots(p1)
print("\nProduct of scaled polynomials\n",p1)
table.foreach(rts, print)

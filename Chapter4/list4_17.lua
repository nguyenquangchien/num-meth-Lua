-- /* File list4_17.lua */ -- Examples of Polynomial operations 

require"Polynomial"

p1 = Polynomial.new(1,3,5) -- List of numbers
p2 = Polynomial.new{8,6,4,2} -- Table of numbers
p3 = Polynomial.new(4*p2) -- Another Polynomial object
p4 = p1^3 -- Powers of Polynomials
p5 = Polynomial.new(-2,0,0,0,0,0,1) -- x^6 - 2

print('p3 =',p3);
print('p2*p1 =',p2*p1)
print('p4/5 =',p4/5)
print('Value of p3(x=1.5) =',p3:value(1.5))
pq,pr = Polynomial.div(p4,p2)

print('Polynomial division of p4 by p2:');
print('Quotient =',pq); print('Remainder =',pr)

-- File list2.5 -- Complex numbers using Complex.lua
 
require"Complex" --Load Complex extensions
Cnew = Complex.new

c1 = Cnew(3,5) -- Define complex numbers
c2 = Cnew(2,8)

print(c1*c2); print(4*c1) -- Test math operations
print(c1/c2); print(c1+c2); print(-c1+c2); print(j*c1)
print(tostring(c1))
print('sin =',Complex.sin(c1))
print('sin =',c1:sin()); print('tan =',c1:tan())


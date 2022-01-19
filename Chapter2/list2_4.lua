-- File list2-4 -- Simple complex number operations

require"list2_3" -- Load defined complex extenstios
Cnew = Cmp.new -- Local define

c1 = Cnew(3,5) -- Define complex numbers
c2 = Cnew(2,8)
j = Cnew(0,1) -- pure imaginary of magnitude unity

print(c1*c2) -- Test math operations
print(c1/c2)
print(c1+c2)
print(-c1+c2)
print(j*c1)
print(tostring(c1))




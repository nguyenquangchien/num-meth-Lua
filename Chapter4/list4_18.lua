 --/* File list4_18.lua */ -- Example of divisor residuals
 
require"Polynomial"

p1 = Polynomial.new{1,-1}*{2,-1}*{3,-1}*{4,-1}
print(p1)

c1,c2 = -3,1
f1,f2 = p1:div({1,c1,c2})
print(f2)
i = 1;xi,v1,v2 = {},{},{}
for c2=1,3,.01 do
	f1,f2 = p1:div({c2,c1,1})
	xi[i],v1[i],v2[i] = c2,f2[1],f2[2]
	i = i+1
end
plot(xi,v1,v2)
write_data('list4_18.dat',xi,v1,v2)
	
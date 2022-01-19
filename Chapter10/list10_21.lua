--/* File list10_21.lua */
-- Programs to solve electric circuit equations using odeivse()
require"odeiv"

R1,R2,R3,R4,R5 = 4e3,2e3,4e3,2e3,2e3
C1,C2 = 1e-6,4e-6; L1 = 2.e-3; vs = 5

f = function(eqs,t,v,vp)
	eqs[1] = -v[4] + v[1]/R1 + (v[1]-v[2])/R2 + C1*(vp[1]-vp[3])
	eqs[2] = (v[2]-v[1])/R2 + v[2]/R3 + (v[2]-v[3])/R4
	eqs[3] = (v[3]-v[2])/R4 + v[3]/R5 + C1*(vp[3]-vp[1]) + C2*vp[3]
	eqs[4] = vs - v[1] - L1*vp[4]
end

s1,err1 = odeivse(f,{0,2e-2},{0,0,0,0})
print('#points=',#s1[1])

print(unpack(maxvalue(s1))); print(unpack(maxvalue(err1)))
plot(s1)

write_data("list10_21.dat",s1,err1)

-- /* File list10_30.lua */
-- Lorentz equations
require"odeiv"

a,b,c = 10,28,8/3
f = function(eqs,t,y,yp)
	eqs[1] = yp[1] + a*(y[1] - y[2])
	eqs[2] = yp[2] + y[2] + y[1]*y[3] - b*y[1]
	eqs[3] = yp[3] + c*y[3] - y[1]*y[2]
end

s1 = odeiv(f,{0,40,10000},{5,5,5})
s2 = odeiv(f,{0,40,10000},{5.0001,5,5})
s3 = odeiv(f,{0,40,10000},{5.000001,5,5})
plot(s1,s2)
write_data("list10_30.dat",s1,s2,s3)


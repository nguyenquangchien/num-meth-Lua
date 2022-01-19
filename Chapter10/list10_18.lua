-- /* File list10_18.lua */
-- Programs to integrate first order diff. equation using adaptive step size algorithm
require"odeiv"

mu = 20
f1 = function(eqs,t,y,yp) -- Stiff differential equation
	eqs[1] = yp[1] + 1001*y[1] - 999*y[2]
	eqs[2] = yp[2] - 999*y[1] + 1001*y[2]
end
f2 = function(eqs,t,y,yp) -- Sinusoidal differential equation
	eqs[1] = yp[1] -y[2]
	eqs[2] = yp[2] +y[1]
end
f3 = function(eqs,t,y,yp) -- Van der Pol equation 
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] - mu*(1 - y[1]^2)*y[2] + y[1]
end
--Now solve three equations
s1,err1 = odeivse(f1,{0,10},{0,2})
print(errstat(err1[2]));print(errstat(err1[3]))
s2,err2 = odeivse(f2,{0,100},{0,1})
print(errstat(err2[2]));print(errstat(err2[3]))
s3,err3 = odeivse(f3,{0,100},{1,0})
print(errstat(err3[2]));print(errstat(err3[3]))
write_data("list10_18.dat",s1,err1,s2,err2,s3,err3)
print(#s1[1],#s2[1],#s3[1])

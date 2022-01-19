-- /* File list10_8.lua */
-- Programs to integrate first order diff. equation using odeivqs()
require"odeiv"

f1 = function(eqs,t,y,yp) -- test of stiff differential eqns
	eqs[1] = yp[1] + 1001*y[1] - 999*y[2]
	eqs[2] = yp[2] - 999*y[1] + 1001*y[2]
end
f2 = function(eqs,t,y,yp) -- test of sinusoidal equation
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] + y[1]
end
mu = 20
f3 = function(eqs,t,y,yp) -- test of Van der Pol equation
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] - mu*(1 - y[1]^2)*y[2] + y[1]
end

s1 = odeivqs(f1,{0,{1.e-5,1e3}},{0,2}); plot(s1)
s2 = odeivqs(f2,{0,{1.e-5,1e3}},{1,0}); plot(s2)
s3 = odeivqs(f3,{0,{1.e-5,1e3}},{1,0}); plot(s3)

write_data("list10_8.dat",s1,s2,s3)

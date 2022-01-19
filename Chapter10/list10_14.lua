-- /* File list10_14.lua */
-- Programs to integrate first order diff. equation using odeive()
require"odeiv"

mu = 10
f1 = function(eqs,t,y,yp) -- Van der Pol equation 
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] - mu*(1 - y[1]^2)*y[2] + y[1]
end
-- TP solution with one step interval, also error evaluation
s1,err = odeive(f1,{0,{100},8000},{2,0})
print(errstat(err[2]));print(errstat(err[3]))

write_data("list10_14.dat",s1,err)
plot(s1[1],s1[2],err[2])

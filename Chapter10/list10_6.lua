-- /* File list10_6.lua */
-- Programs to integrate first order diff. equation using odebiv()

require"odeiv"

mu = 20
feqs = function(eqs,t,y,yp) -- test of stiff differential eqns
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] - mu*(1 - y[1]^2)*y[2] + y[1]
end

yin = {1,0} -- Initial values, 0 and 2
sol,nit = odebiv(feqs,{0,100,16000},yin)
print('Number of iterations =',nit); plot(sol[1],sol[2])
write_data("list10_6.dat",sol)

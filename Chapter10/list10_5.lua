-- /* File list10_5.lua */
-- Programs to integrate first order diff. equation using odebiv()

require"odeiv"

feqs = function(eqs,t,y,yp) -- test of stiff differential eqns
	eqs[1] = yp[1] + 1001*y[1] - 999*y[2]
	eqs[2] = yp[2] - 999*y[1] + 1001*y[2]
end

yin = {0,2} -- Initial values, 0 and 2
sol,nit = odebiv(feqs,{0,1,1000},yin)
print(nit); plot(sol)
write_data("list10_5.dat",sol)

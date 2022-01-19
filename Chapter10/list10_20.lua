---[[ /* File list10_20.lua */
-- Programs to integrate first or second order differential equations
require"odeiv"
odebiv = odeb12 -- Substitute 2nd order solver
-- mu = 20 Try different values
mu = 0 -- Gives pure sinusoidal solutions
-- Example of second order formulation
f1 = function(eqs,x,y,yp,ypp) -- Van der Pol equation 
	eqs[1] = ypp[1] -mu*(1-y[1]^2)*yp[1] + y[1]
end
-- Example of equivalent first order formulation
f2 = function(eqs,t,y,yp)
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] - mu*(1 - y[1]^2)*y[2] + y[1]
end
t1 = os.clock()
s1,err1 = odeive(f1,{0,100,8000},{1},{0})
print('first time =',os.clock()-t1); plot(s1)
t1 = os.clock()
s2,err2 = odeive(f2,{0,100,8000},{1,0})
print('second time =',os.clock()-t1); plot(s2[1],s2[2])

write_data("list10_20,dat",s1,s2,err1,err2)


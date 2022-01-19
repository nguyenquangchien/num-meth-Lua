-- /* list11_17.lua */
-- Solution of nonlinear BV problem with non-uniform spatial grid
require"odebvfd"

x = xlg(0,1,1.e-6,2000) -- Log distribution -- Try others below
--x = xgp(0,1,1.e-6,2000) -- Geometric distribution
--x = xll(0,1,1.e-6,2000) -- Log-log distribution

f = function(x,y,yp,ypp) -- Differntial equation
	return ypp + 4*yp^2
end

u,nn,err1,err2 = ode1fd(f,x,{0,2}) -- Solve equation

print(nn,err1,err2) -- # Newton iterations and errors
plot(u); write_data('list11_17.dat',u)

-- /* list11_18.lua */
-- Solution of nonlinear BV problem with non-uniform spatial grid
require"odebvfd"

f = function(x,y,yp,ypp) -- Differntial equation
	return ypp + 4*yp^2
end

x = xgp(0,1,1.e-6) -- Geometric distribution
s,err,nn = ode1fde(f,x,{0,2}) -- Solve equation
print('Number of Newton iterations = ',nn) -- # Newton iterations and errors
plot(s); write_data('list11_18.dat',s,err); plot(err)

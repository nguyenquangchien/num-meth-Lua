-- /* list11_19.lua */
-- Solution of nonlinear BV problem with non-uniform spatial grid
require"odebvfd"
require"intp"

f = function(x,y,yp,ypp) -- Differntial equation
	return ypp + 4*yp^2
end

xx = xgp(0,1,1.e-6,100) --Generate approximate solution, few grid points
s1,nn = ode1fd(f,xx,{0,2}) -- Solve equation
plot(s1);print('Number of Newton iterations = ',nn) -- # Newton iterations
x = {0} -- New spatial array, from approximate solution
for i=2,2001 do -- Expand to 2001 grid points
	y = 2*(i-1)/2000; x[i] = intp(s1[2],s1[1],y)
end
s,err,nn = ode1fde(f,x,{0,2}) -- More accurate solution
print('Number of Newton iterations = ',nn)
plot(s); write_data('list11_19.dat',s,err,s1); plot(err)

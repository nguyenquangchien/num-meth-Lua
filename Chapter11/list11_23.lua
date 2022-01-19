-- /* list11_23.lua */
-- Solution of nonlinear BV problem with coupled equation code

require"odefd"

f = function(eqs,x,y,yp,ypp) -- Differntial equation
	eqs[1] = ypp[1] + 4*yp[1]^2
end

u = {0,2}; x = {0,1,2000}
s,nn,err1 = odefd(f,x,u)
print(nn,err1[1])
plot(s[1],s[2]); write_data('list11_23.dat',s)


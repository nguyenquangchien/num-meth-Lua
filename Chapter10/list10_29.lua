-- /* File list10_29.lua */
--Shape of lamp reflector 

require"odeiv"

f = function(eqs,x,y,yp,ypp)
	eqs[1] = x*yp[1]^2 - 2*y[1]*yp[1] - x
end

s1 = odeivs(f,{-2,2},{0}) 
plot(s1); write_data("list10_29.dat",s1) 


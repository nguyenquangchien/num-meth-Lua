--- /* File list10_22.lua */
-- Programs to integrate second order differential equation with sliding friction

require"odeiv"
getfenv(odeiv).odebiv = odeb12

g = 983.21
w = 100
u = 0.02
eps = 1.e-2
ug = u*g

f1 = function(eqs,t,x,xp,xpp)
	eqs[1] = xpp[1] + ug*xp[1]/math.sqrt(xp[1]^2+eps) + w*x[1]
end

s1,err1 = odeivse(f1,{0,9},{10},{0})
plot(s1[1],s1[2],err1[2])
write_data("list10_22.dat",s1)
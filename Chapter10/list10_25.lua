--- /* File list10_25.lua */
-- Programs to integrate predator prey equations
require"odeiv"

alf1,alf2,del1,del2 = 2.0, 0.0002, 0.02, .8

f = function(eqs,t,p,pp)
	eqs[1] = pp[1] - alf1*p[1] + del1*p[1]*p[2]
	eqs[2] = pp[2] - alf2*p[1]*p[2] +del2*p[2]
end

s1 = odeivs(f,{0,30},{5000,100})
s2 = odeivs(f,{0,30},{5000,200})
s3 = odeivs(f,{0,30},{5000,300})
plot(s1);plot(s2);plot(s3)
write_data("list10_25.dat",s1,s2,s3)


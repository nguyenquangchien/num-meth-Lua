--- /* File list10_24.lua */
-- Programs to integrate second order F=ma equations
require"odeiv"
odebiv = odeb12

f = function(eqs,t,f,fp,fpp)
	eqs[1] = fpp[1] + f[1]/(f[1]^2 + f[2]^2)^1.5
	eqs[2] = fpp[2] + f[2]/(f[1]^2 + f[2]^2)^1.5
end

s1 = odeivs(f,{0,20},{0,1},{1,0})
s2 = odeivs(f,{0,20},{0,1},{.2,0})
s3 = odeivs(f,{0,20},{0,1},{1.2,0})
plot({s1[2],s1[3]},{s2[2],s2[3]},{s3[2],s3[3]})
write_data("list10_24.dat",s1,s2,s3)

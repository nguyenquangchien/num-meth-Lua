-- /* File list10_12.lua */
-- Programs to integrate first order diff. equation using odeive()
require"odeiv"

f1 = function(eqs,t,y,yp) -- test of sinusoidal differential eqns
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] + y[1]
end
-- TP solution with two step intervals, also error evaluation
s1,err = odeive(f1,{0,{.01,40},{100,2000}},{1,0})
print(errstat(err[2]));print(errstat(err[3]))
-- Calculate exact solution, corrected solution and actual errors obtained 
yexact,ycorrected,err1,err2,nd = {},{},{},{},#s1[1]
for i=1,nd do yexact[i] = math.cos(s1[1][i]) end
for i=1,nd do ycorrected[i] = s1[2][i] - err[2][i] end
for i=1,nd do 
	err1[i] = math.abs(yexact[i] - s1[2][i]) -- Error in TP solution
	err2[i] = math.abs(yexact[i] - ycorrected[i]) -- Error in corrected solution
end
print(errstat(err1));print(errstat(err2))
write_data("list10_12.dat",s1,ycorrected,err,err1,err2)
write_data(20,"list10_121.dat",s1,ycorrected,err,err1,err2)
plot(err)


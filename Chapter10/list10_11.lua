-- /* File list10_11.lua */
-- Programs to integrate first order diff. equation using odeiv()
require"odeiv"; exp = math.exp; abs = math.abs

f1 = function(eqs,t,y,yp) -- test of stiff differential eqns
	eqs[1] = yp[1] + 1001*y[1] - 999*y[2]
	eqs[2] = yp[2] - 999*y[1] + 1001*y[2]
end
-- TP solution with Quick Scan function, also error evaluation
s1,err = odeivqse(f1,{0,{1.e-6,10},400},{0,2})
print(errstat(err[2]));print(errstat(err[3]))
-- Calculate exact solution, corrected solution and actual errors obtained 
yexact,ycorrected,err1,err2,nd = {},{},{},{},#s1[1]
for i=1,nd do yexact[i] = exp(-2*s1[1][i])+exp(-2000*s1[1][i]) end
for i=1,nd do ycorrected[i] = s1[3][i] - err[3][i] end
for i=1,nd do 
	err1[i] = abs(yexact[i] - s1[3][i]) -- Error in TP solution
	err2[i] = abs(yexact[i] - ycorrected[i]) -- Error in corrected solution
end
print(errstat(err1));print(errstat(err2))
write_data("list10_11.dat",s1,err,err1,err2)
write_data(20,"list10_111.dat",s1,err,err1,err2)

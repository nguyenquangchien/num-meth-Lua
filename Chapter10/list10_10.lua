-- /* File list10._0.lua */
-- Programs to integrate first order diff. equation using odeiv()
require"odeiv"; exp = math.exp; abs = math.abs

f1 = function(eqs,t,y,yp) -- test of stiff differential eqns
	eqs[1] = yp[1] + 1001*y[1] - 999*y[2]
	eqs[2] = yp[2] - 999*y[1] + 1001*y[2]
end
-- TP Solution with step sizes of 2h and h
s1 = odeiv(f1,{0,{1.e-4,1.e-3,1.e-2,0.1,1,10},200},{0,2})
s2 = odeiv(f1,{0,{1.e-4,1.e-3,1.e-2,0.1,1,10},400},{0,2})

-- Evaluate estimated error given two solutions with h and 2h step sizes
erra = odeerror(s1,s2) -- Evaluation of estimated error
print(errstat(erra[2])) -- errstat() returns statistics of error estimate
print(errstat(erra[3])) -- returns RMS, maximum, point #, average for errors

-- Calculate exact solution, corrected solution and actual errors obtained 
yexact,ycorrected,err1,err2,nd = {},{},{},{},#s2[1]
for i=1,nd do yexact[i] = exp(-2*s2[1][i])+exp(-2000*s2[1][i]) end
for i=1,nd do ycorrected[i] = s2[3][i] - erra[3][i] end
for i=1,nd do 
	err1[i] = abs(yexact[i] - s2[3][i]) -- Error in TP solution
	err2[i] = abs(yexact[i] - ycorrected[i]) -- Error in corrected solution
end
print(errstat(err1))
print(errstat(err2))
plot(s1[1],err1)

-- Solve ODE by TP and return solution and error estimate -- odeive()
s3,err,n2 = odeive(f1,{0,{1.e-4,1.e-3,1.e-2,0.1,1,10},400},{0,2})
print(errstat(err[2]))
print(errstat(err[3]))
write_data("list10_10.dat",s2,s1,s3,erra,err,err1,err2)
write_data(20,"list10_101.dat",s2,s1,s3,erra,err,err1,err2)

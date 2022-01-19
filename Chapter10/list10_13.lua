-- /* File list10_13.lua */
-- Programs to integrate first order diff. equation using odeive()
require"odeiv"

--odebiv = odebrk
f1 = function(eqs,t,y,yp) -- test of nonlinear differential eqns
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2]  - math.exp(y[1])
end
sr2 =math.sqrt(2)
xmax = 2.221 -- must be less than 2.22144
-- TP solution with two step intervals, also error evaluation
s1,err1,n1 = odeive(f1,{0,{xmax*.99,xmax},1000},{0,0})
print(errstat(err1[2]));print(errstat(err1[3]));print('n = ',n1)
s2,err2,n2 = odeive(f1,{0,{xmax*.99,xmax},2000},{0,0})
print(errstat(err2[2]));print(errstat(err2[3]));print('n2 = ',n2)

y1,y2 = {},{} -- Evaluate exact solutions
for i=1,#s2[1] do 
	y1[i] = math.log(1/math.cos(s2[1][i]/sr2)^2) 
	y2[i] = sr2*math.tan(s2[1][i]/sr2)
end
write_data("list10_13.dat",s1,s2,y1,y2,err1,err2)
write_data(50,"list10_131.dat",s1,s2,y1,y2,err1,err2)
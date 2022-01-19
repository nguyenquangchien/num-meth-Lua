-- /* File list10_9.lua */
-- Programs to integrate first order diff. equation using odeiv()
require"odeiv"

f1 = function(eqs,t,y,yp) -- test of stiff differential eqns
	eqs[1] = yp[1] + 1001*y[1] - 999*y[2]
	eqs[2] = yp[2] - 999*y[1] + 1001*y[2]
end

s1 = odeiv(f1,{0,{1.e-4,1.e-3,1.e-2,0.1,1,10},400},{0,2}); plot(s1)
--s1 = odeivqs(f1,{0,{1.e-6,1.e3},400},{0,2}); plot(s1)
st1,st2 = {},{}
nd = #s1[1]
err,terr = 0,0
exp,abs = math.exp,math.abs
for i=1,nd do 
	time = s1[1][i]
	st1[i] = exp(-2*time)-exp(-2000*time)
	st2[i] = exp(-2*time)+exp(-2000*time)
	err1 = abs(st1[i] - s1[2][i])
	if err1 > err then err,terr = err1,time end
	err1 = abs(st2[i] - s1[3][i])
	if err1 > err then err,terr = err1,time end
end
print('Maximum error is ',err,' and occurs at time ',terr)
write_data("list10.9_dat",s1,st1,st2)

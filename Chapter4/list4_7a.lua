-- File list4_9a.lua Program for the voltages of a Full Wave Rectifier
-- Exponential functions present special convergence problems

require"nsolv"
exp = math.exp; pi2 = math.pi*2 -- Simpler to type

Is,Rl,Vt = 1.e-14,5000,.026 -- Basic parameters
Vm = 5 -- Magnitude of source voltage -- change as desired
getfenv(nsolv).NMAX = 200
NMAX = getfenv(nsolv).NMAX -- Get max number of iterations

fwr = function(I,v)
	I[1] = (v[1]-v[2])/Rl - Is*(exp((vs-v[1])/Vt)-1) - Is*(exp(-v[1]/Vt)-1) 
	I[2] = (v[2]-v[1])/Rl + Is*(exp((v[2]-vs)/Vt)-1) + Is*(exp(v[2]/Vt)-1) 
end

step = {-.8,-.8} -- Not too critical .2 to 1.2 seems to work fine
x = {0,0} -- Just need value to get started
xa,v1,v2,v3,vout = {},{},{},{},{} -- Arrays for values
for i=1,401 do
		xt = (i-1)*pi2/400; xa[i] = xt
		vs = Vm*math.sin(xt)
		nit,err = nsolv(fwr,x,step) -- Use previous solution as guess
		if nit==NMAX then print('Convergence error at i = ',i,err) end
		v3[i],v1[i],v2[i] = vs, x[1], x[2]
		vout[i] = x[1]-x[2] -- Output voltage
end
write_data("list4_9a.dat",xa,v1,v2,v3,vout)
plot(xa,v3,vout)
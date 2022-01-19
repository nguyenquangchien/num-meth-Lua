-- /* File list12_11.lua */
-- Program for transient solution of p-n junction

require"odefd"; require"pdeivbv"
getfenv(pdeivbv).nprint = 1
getfenv(ode2bvfd).umin = {5.e-4,5.e-4,5.e-4}
--getfenv(ode2bvfd).nprint = 1; getfenv(pdebivbv).nprint = 1

-- Model equations to be solved
L = 8.e-4; q = 1.6e-19; eps = 11.9*8.854e-14
Na = 1e16; Nd = 1.e19 -- P and N type doping
vth = .026; ni = 1.45e10
tno = 1.e-8; tpo = 2.e-8
unp = 1020 -- Minority carrier mobility in p region
unn = 100 -- Majority carrier mobility in n region
upn = 50 -- Minirity carrier mobility in n region
upp = 400 -- Majority carrier mobility in p region
vsat = 1.e7 -- Saturated velocity
qdep = q/eps
vj = vth*math.log(Na*Nd/ni^2)
va = 0

eq = function(fv,x,t,v,vp,vpp,vt) -- Equations with time and spatial derivatives
	if x<=0 then un,up,Nnet = unp,upp,-Na
	else un,up,Nnet = unn,upn,Nd end
	p = Na*math.exp((v[2]-v[1])/vth)
	n = (ni^2/Na)*math.exp((v[1]-v[3])/vth)
	U = (p*n-ni^2)/(tpo*(n+ni)+tno*(p+ni))
	fv[1] = vpp[1] + qdep*(p - n + Nnet) -- V term - Poisson's equation
	fv[2] = vt[1] - vt[2] + up*(vth*vpp[2] + vp[2]*(vp[2]-vp[1])) - vth*U/p-- Up term for holes
	fv[3] = vt[1] - vt[3] + un*(vth*vpp[3] + vp[3]*(vp[1]-vp[3])) + vth*U/n-- Un term for electrons
	end

efl = function(fv,v,vp) -- Using mixed boundary condition
	fv[1] = v[1] - va
	fv[2] = v[2] - va
	fv[3] = (1+unp/vsat*(vp[3]))-math.exp((v[3]-v[1])/vth)
end
efr = function(fv,v,vp) -- Using mixed boundary condition
	fv[1] = v[1] - vj 
	fv[2] = (1+upn/vsat*(vp[2]))-math.exp((v[1]-vj-v[2])/vth)
	fv[3] = v[3] 
end

ftzero = {0,0,0} -- Force zero time derivatives for Steady State
eqss = function(fv,x,v,vx,vxx) -- Steady State equations with zero time derivatives
	eq(fv,x,0,v,vx,vxx,ftzero)
end
x1 = {} -- Set up nonuniform spatial array
fact,step = 1.02,1.e-3; n1 =300; np1 = n1+1
x1[1] = 0
for i=2,np1 do -- Small to large steps
	x1[i] = x1[i-1] + step; step = step*fact
end
j = np1
for i=np1+1,2*n1+1 do -- Reverse array
	x1[i] = x1[i-1] + x1[j] - x1[j-1]; j = j-1
end
x = {}; j = 1
imax = #x1; xmax = x1[imax]
for i=1,imax do -- Scale to -2L/3 to L/3
	x[j] = (2*L/3)*(x1[j]/xmax - 1); j = j+1
end
for i=2,imax do
	x[j] = (L/3)*x1[i]/xmax; j = j+1
end
nx = #x; nxm1 = nx-1 -- End of array definition

v = {{},{},{}} -- Set initial voltage values
for i=1,nx do
	for j=1,3 do
		if j==1 then 
			if x[i]<0 then v[1][i] = 0 
			else v[1][i] = vj end
		else v[j][i] = 0 end
	end
end
va = 0
s= odefd({eqss,efl,efr},x,v) -- DC solution -- initial values

va = 1 -- Step input voltage of 1 volt
for i=1,nx do -- Add to existing voltages
	duu = va-(x[i]-x[1])/(x[nx]-x[1])*va
	for j=1,3 do v[j][i] = v[j][i] + duu end
end
--tvals = {0,{1.e-15,1.e-6},2,10} -- 20 steps/decade, save 2
tvals = {0,{1.e-15,1.e-6},2,5} -- 10 steps/decade, save 2
s = pdeivbvqs({eq,efl,efr},tvals,x,v)
write_data('list12_11.dat',s)	
	

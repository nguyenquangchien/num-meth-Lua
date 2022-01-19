-- /* File list11_28.lua */ -- Solution of semiconductor device equations

require"odefd"
getfenv(odefd).umin = {1.e-6, 1.e-6, 1.e-6}
getfenv(odefd).nprint=1

-- Material and device parameters
q = 1.6e-19; eps = 11.9*8.854e-14; L = 2.e-4
vt = .026; ni = 1.45e10
Na = 1e16; Nd = 1.e19 -- Doping densities
tno = 1.e-8; tpo = 2.e-8 -- Lifetimes
n1 = ni; p1 = ni
unp,unn = 1020, 100 -- Electron mobilities
upp,upn = 420, 50 -- Hole mobilities
vsat = 1.e7 -- Saturated drift velocity
qdep = q/eps; vj = vt*math.log(Na*Nd/ni^2); no = ni^2/Na
va = 0.0

eq = function(fv,x,v,vp,vpp)
	if x<=0 then Nnet, un, up = -Na, unp, upp
	else Nnet, un, up = Nd, unn, upn end
	p, n  = Na*math.exp((v[2]-v[1])/vt), no*math.exp((v[1]-v[3])/vt)
	U = (p*n-ni^2)/(tpo*(n+n1)+tno*(p+p1))
	fv[1] = vpp[1] + qdep*(p - n + Nnet) -- V term - Poisson's equation
	fv[2] = vpp[2] + vp[2]*(vp[2]-vp[1])/vt - U/(p*up)-- Up term for holes
	fv[3] = vpp[3] + vp[3]*(vp[1]-vp[3])/vt + U/(n*un)-- Un term for electrons
end

efl = function(fv,v,vp) -- Left boundary values
	fv[1], fv[2] = v[1]-va, v[2]-va
	fv[3] = math.exp((v[1]-v[3])/vt)*(1 + unp/vsat*vp[3]) - 1 
end
efr = function(fv,v,vp) -- Right boundary values
	fv[1], fv[3] = v[1] - vj, v[3]
	fv[2] = math.exp(v[2]/vt)*(1 + upn/vsat*vp[2]) - 1
end

x,xi, v, nx = {},{}, {}, 2001; dx = L/(nx-1)
for j=1,3 do v[j] = {} end
for i=1,nx do -- Set initial values -- step in v[1][]
	x[i] = (i-1)*dx - 3*L/4
	xi[i] = i-1
	if x[i]<0 then v[1][i] = 0 else v[1][i] = vj end -- Step in potential
	for j=2,3 do v[j][i] = 0 end -- Zero Quasi-Fermi levels
end
for i=0,11 do
	va = i/10
	s,nn = odefd({eq,efl,efr},x,v)
	print('Number of iterations is ',nn,' at va = ',va); io.flush()
	write_data('list11_28'..i..'.dat',s) -- Save solution
	plot(x,unpack(v))
	x,v = xad(x,v) 
end

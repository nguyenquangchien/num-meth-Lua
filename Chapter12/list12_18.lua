-- /* File list12_18.lua */
-- Program for transient solution of Schrodinger Wave Equation

require"Complex"; math.abs = Complex.abs -- Complex algebra
require"odefd"; require"pdeivbv"
require"intp"; require"intg"
getfenv(pdeivbv).nprint = 1 -- monitor progress
getfenv(pdebbviv).umin = {5.e-4,5.e-4,5.e-4}
getfenv(odefd).NMAX = 1

tomagsq = function(s) -- Convert to magnitude squared
	ns,nv = #s, #s[2]
	for i=2,ns do 
		for k=1,nv do s[i][k] = math.abs(s[i][k])^2 end
	end
	return s
end
sqarea = function(x,y) -- Calculate area under wavefunction
	local nx, sum = #x-1, 0
	for i=1,nx	do -- works for nonumiform spatial grid
		sum = sum + (x[i+1]-x[i])*(math.abs(y[i+1])^2+math.abs(y[i])^2)
	end
	return sum/2
end

-- Normalized Schrodinger equation to be solved
eq = function(fv,x,t,v,vp,vpp,vt) 
	fv[1] = j*vt[1] + vpp[1] -- V =0 term
end
efl = function(fv,v,vp) -- Zero boundary value
	fv[1] = v[1]
end
efr = function(fv,v,vp) -- Zero boundary value
	fv[1] = v[1]
end

x = {}; nx = 1000 -- Set initial parameters
L, alf, bet = 1, 75, 100
xo = .5*L
for i=1,nx+1 do x[i] = L*(-1 + (i-1)*2/nx) end

v = {{}} -- Set initial voltage values
vp,vpp = {{}}, {{}}
for i=1,nx+1 do
	xv = x[i]
	v[1][i] = math.exp(-alf*(xv+xo)^2)*Complex.exp(bet*j*xv)
end
v[1][1] = Complex.new(0,0); v[1][nx+1] = v[1][1]

sum =sqarea(x,v[1])^0.5; print('P = ',sum)
for i=1,nx+1 do v[1][i] = v[1][i]/sum end
print('Now P = ',sqarea(x,v[1]))

tvals = {0,{1.e-6,1.e-5},1,5} -- Starter solution, to conserve P
s = pdeivbvqs({eq,efl,efr},tvals,x,v,vp,vpp)
sum = sqarea(s[1],s[3]); print('After initial short interval, P = ',sum)

tvals = {0,.08,160,10,{-.75,-.5,-.25,0,.25,.5,.75}}
s,st = pdeivbvt({eq,efl,efr},tvals,x,v,vp,vpp)
sum = sqarea(s[1],s[#s]); print('After final solution, P = ',sum)

write_data('list12_18.dat',tomagsq(s)) -- Save data as probability density
write_data('list12_18t.dat',tomagsq(st))

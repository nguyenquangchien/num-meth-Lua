-- /* list11_20.lua */
-- Solution of semiconductor depletion layer as a BV problem
require"odebvfd"; exp = math.exp

L,Na,Nd = 2.e-5, 1e18, 1e18 -- Size and doping densities
q,eps,vt,ni = 1.6e-19, 11.9*8.854e-14, .026, 1.45e10
qdep = q/eps
v1,v2 = -vt*math.log(Na/ni), vt*math.log(Nd/ni) -- Boundary values

f = function(x,v,vp,vpp) -- Poisson's equation
	if x<=0 then Nnet = -Na else Nnet = Nd end
	return vpp + qdep*(ni*exp(-v/vt) - ni*exp(v/vt) + Nnet)
end

nx = 2001; dx = L/(nx-1); x,v = {},{}
for i=1,nx do -- Initial voltage approximation
	x[i] = (i-1)*dx - L/2
	if x[i]<0 then v[i] = v1 else v[i] = v2 end
end
s,err,nn = ode1fde(f,x,v) -- Solve equations, fixed boundary values
p,n = {},{} -- To calculate holes and electrons
for i=1,nx do -- Now calculate then
	p[i],n[i] = ni*exp(-s[2][i]/vt),ni*exp(s[2][i]/vt) end

print('Number of Newton iterations = ',nn)
plot(s); plot(s[1],n);plot(s[1],p)
plot(err); write_data('list11_20.dat',s,err,n,p)

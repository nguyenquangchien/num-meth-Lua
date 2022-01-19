-- File list12_34.lua --
-- Example of BV problem for p-n junction
require"pde2bv"
-- Material and device parameters
q = 1.6e-19; eps = 11.9*8.854e-14; L = 0.4e-4
vt = .026; ni = 1.45e10
Na = 1e19; Nd = 1.e17-- Doping densities
x1,y1,x2 = L/2, L/4, L*3/4; va = 0.0
qdep = q/eps; vj = vt*math.log(Na*Nd/ni^2); no = ni^2/Na

feqs = {
	function(x,y,uxx,ux,u,uy,uyy,i,j)
		if x>=x1 and y<=y1 then Nnet = -Na
		else Nnet = Nd end
		p, n  = Na*math.exp(-u/vt), no*math.exp(u/vt)
		return uxx + uyy + qdep*(Nnet + p - n)
	end,
	function(x,u,uy,ux,i) 
		if x>x2 then return u-va else return uy end
	end,
	function(x,u,uy,ux,i) return u-vj end,
	function(y,u,ux,uy,j) return ux end,
	function(y,u,ux,uy,j) return ux end
}	
	
Nx,Ny = 80,80; nx,ny = Nx+1,Ny+1
nyj = math.floor(Ny/4); nxj = math.floor(Nx/2)
x,y = setxy({0,L},{0,L},Nx,Ny)
u = Spmat.new(nx,-ny) 
for j = 1,ny do	-- Set initial values at 0 or Vj
	yv = y[j]
	u[j] = {}; for i = 1,nx do 
		xv = x[i]
		if xv>=x1 and yv<=y1 then u[j][i] = 0.0
		else u[j][i] = vj end
	end
end
getfenv(pde2bv).nprint = 1; --getfenv(pde2bvcoe).nprint=1
SPM,COE,SOR,ADI = 1, 2, 3, 4 -- 4 solution types 
t1 = os.clock()
u,errm = pde2bv(feqs,x,y,u,COE) -- Replace COE as desired
print('time =',os.clock()-t1)

pa,na = {},{}; ut = reversexy(u)
for i=1,nx do -- Calculate carrier densities
	paa,naa = {}, {}
	for j=1,ny do 
		ua = ut[i][j]
		paa[j] = math.log10(Na*math.exp(-ua/vt))
		naa[j] = math.log10(no*math.exp(ua/vt))
	end
	pa[i],na[i] = paa, naa
end
splot(ut); cplot(ut); splot('list12_34.emf',ut)
write_data('list12_34u.dat',ut)
write_data('list12_34p.dat',pa); write_data('list12_34n.dat',na)

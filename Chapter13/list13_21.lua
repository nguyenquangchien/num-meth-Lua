-- File list13_21.lua --
-- Example of BV problem for p-n junction by FE method
require"pde2fe"; exp = math.exp

-- Material and device parameters
q = 1.6e-19; eps = 11.9*8.854e-14; L = 0.4e-4
vt = .026; ni = 1.45e10
Na = 1e19; Nd = 1.e17-- Doping densities
x1,y1,x2 = L/2, L/4, L*3/4; va = 0.0
qdep = q/eps; vj = vt*math.log(Na*Nd/ni^2); no = ni^2/Na

nts = {readnts('list13_21',L)} -- Read spatial data
nds = nts[1] -- nodes data
npts,u = #nds, {} -- number of nodes
for i=1,npts do -- set initial value & scale dimensions
	xt,yt = nds[i][1],nds[i][2]
	if xt>=0.5*L and yt<=0.25*L then u[i] = 0.0
	else u[i] = vj end -- Set initial potential values
end

feqs = { -- Equation to be solved 
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		if mtn==2 then Nnet = -Na else Nnet = Nd end
		p, n  = Na*exp(-u/vt), no*exp(u/vt)
		return uxx + uyy + qdep*(Nnet + p - n)
	end,
	function(nd,u,un,nbs,kb) -- Boundary values
		if nbs==3 then return u-vj -- Top voltage contact
		elseif nbs==2 then return u-va -- Bottom contact
		else return un end -- Else zero normal derivative
	end
}	
	
getfenv(pde2fe).nprint = 1
SPM,COE,SOR = 1, 2, 3 -- 3 solution methods 

u,errm = pde2fe(feqs,nts,u,COE) -- Solve by FE method

x,y = {},{}; NT = 81 
for i=1,NT do -- Define uniform x-y grid for plotting
	x[i] = (i-1)*L/(NT-1); y[i] = x[i]
end
ut = toxysol(u,x,y,nts[1],nts[2])
pa,na = {},{} -- Calculate hole and electron densities
for i=1,#ut[3] do
	ua = ut[3][i]
	pa[i] = math.log10(Na*exp(-ua/vt))
	na[i] = math.log10(no*exp(ua/vt))
end
splot(ut); write_data('list13_21.dat',ut,pa,na)

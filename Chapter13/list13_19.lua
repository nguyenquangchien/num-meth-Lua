-- File list13_19.lua -- Time dependent nonlinear diffusion

require'pde2fe' -- Input FE code
L = 1.e-4
nts ={readnts('list13_9',L)} -- Get node, triangle, sides data
-- Model equations to be solved
D00,D10,E0,E1 = 0.05, 0.95, 3.5, 3.5 -- Diffusion coeff parameters
T = 1000+273 -- temperature
D0 = D00*math.exp(-E0/(0.026*T/300))
D1 = D10*math.exp(-E1/(0.026*T/300))
ni = 7.14e18; Un = 5e20
Um = 1.0 -- Use normalized value
xmax,ymax = 2*L, L
feqs = { -- Table of functions
	function(x,y,t,uxx,ux,u,uy,uyy,ut,utt) -- General point
		D = D0 + D1*Un*u/ni
		return D*(uxx + uyy) - ut 
	end,
	function(nd,t,u,un,ut,utt,nbs)  -- Boundary values
		if nbs==3 then return u - Um 
		else return un end
	end
}  -- End general point and boundary values
u = {}
for k=1,#nts[1] do 
	if nts[1][k][3] == 3 then u[k] = Um
	else u[k] = 0.0 end
end	

SPM,COE,SOR = 1, 2, 3 -- 3 solution methods 
t = os.time()
getfenv(pde1stp2fe1t).nprint=1
getfenv(pde2fe1tqs).nprint=1
getfenv(pde2fe).nprint=1
tvals = {0,{1,1000},2,20} 
u,errm = pde2fe1tqs(feqs,tvals,nts,u,SOR)
print('Time taken =',os.time()-t); io.flush()

x,y = {},{}; NT = 81 -- Change as desired
for j=1,NT do y[j] = ymax*(j-1)/(NT-1); x[j] = 2*y[j] end
nsol = #u
for i=1,nsol do -- Save 2D data files
	sfl = 'list13_19.'..i..'.dat'
	sol,solxy = toxysol(u[i],x,y,nts[1],nts[2])
	write_data(sfl,sol) 
	if i==nsol then 
		splot('list13_19a.emf',solxy); cplot('list13_19b.emf',solxy)
	end
end

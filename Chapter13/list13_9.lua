-- File list13_9.lua -- Example of time dependent BV problem in 2 dimensions 

require'pde2fe' -- Input FE code
nts ={readnts('list13_9')} -- Get node, triangle, sides data
print('Number of nodes, triangles =',#nts[1],#nts[2])
Um = 1.0 -- Define equations to be solved
D = 20 -- Constant diffusion coefficient
xmax,ymax = 2.0,1.0
feqs = { -- Table of functions
	function(x,y,t,uxx,ux,u,uy,uyy,ut,utt) -- General point
		return -ut + D*(uxx + uyy)
	end,
	function(nd,t,u,un,ut,utt,nbs)  -- Boundary values
		if nbs==3 then return u - Um 
		else return un end
	end
}  -- End general point and boundary values

u = {}; for k=1,#nts[1] do u[k] = 0.0 end
x,y = {},{}; NT = 41 -- Change as desired
for j=1,NT do y[j] = ymax*(j-1)/(NT-1); x[j] = 2*y[j] end

SPM,COE,SOR = 1, 2, 3 -- 3 solution methods 
t = os.time() -- Time results
getfenv(pde2fe1tqs).nprint=1 -- Observe progress in time?
--getfenv(pde2bvcoe).nprint=1 -- Observe COE convergence?
getfenv(pde2fe).nprint=1 -- Observe Newton results?
--getfenv(pde2fe1tqs).seeplot = {x,y} -- See popup plots?
tvals = {0,{1.e-4,1.e-1},2,10, -- Spatial points for saving on next line
	{{0,.5},{0,1},{2,0},{2,.75},{2,1},{1,.5},{1,.75}}}
u,uxyt,errm = pde2fe1tqs(feqs,tvals,nts,u,COE)

print('Time taken =',os.time()-t); io.flush()
nsol = #u
for i=1,nsol do -- Save 2D data files at saved times
	sfl = 'list13_9.'..i..'.dat'
	sol,solxy = toxysol(u[i],x,y,nts[1],nts[2]) -- to x-y arrays
	write_data(sfl,sol) -- Save x-y arrays
	if i==5 then -- Save popup plots at one time
		splot('list13_9a.emf',solxy); cplot('list13_9b.emf',solxy)
	end
end
write_data('list13_9a.dat',uxyt) -- Save time data

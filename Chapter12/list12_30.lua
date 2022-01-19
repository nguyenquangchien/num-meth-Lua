-- File list12_30.lua -- Example of time dependent BV problem in 2 dimensions 

require"pde2bv" 
Um = 1.0 -- Define equations to be solved
D = 20 -- Constant diffusion coefficient
xmax,ymax = 2.0,1.0
feqs = { -- Table of functions
	function(x,y,t,uxx,ux,u,uy,uyy,ut,utt,i,j) -- General point
		return -ut + D*(uxx + uyy)
	end,
	function(x,t,u,uy,ux,i) return uy end, -- Bottom boundary
	function(x,t,u,uy,ux,i) -- Top boundary
		if x>=xmax/2 then return u-Um else return uy end end,
	function(y,t,u,ux,uy,j) return ux end, -- Left boundary
	function(y,t,u,ux,uy,j) return ux end  -- Right boundary
}  -- End general point and boundary values

Nx,Ny = 100,50 -- Change as desired
nx,ny = Nx+1,Ny+1
x,y = setxy({0,xmax},{0,ymax},Nx,Ny)
u = {}; u[0] = {nx, ny, 1} -- Include size information 
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end
SPM,COE,SOR,ADI = 1, 2, 3, 4 -- possible solution methods 
t = os.time()
getfenv(pde2bv1tqs).nprint=1
getfenv(pde2bv1tqs).seeplot = 2 -- Uncomment for popup plots
tvals = {0,{1.e-4,1.e-1},2,10,{{0,.5},{0,1},{2,0},{2,.75},{2,1},{1,.5},{1,.75},}}
u,uxyt,errm = pde2bv1tqs(feqs,tvals,x,y,u,COE)
--tvals = {0,{.001,.01},5} -- Use these for linear time intervals
--u,errm = pde2bv1t(feqs,tvals,x,y,u,COE)

print('Time taken =',os.time()-t); io.flush()
nsol = #u -- Number of time solutions saved
for i=1,nsol do -- Save 2D data files
	sfl = 'list12_30.'..i..'.dat'
	write_data(sfl,reversexy(u[i])) -- note reversexy() before save
end
splot('list12_30a.emf',reversexy(u[5]))
cplot('list12_30b.emf',reversexy(u[5]))
write_data('list12_30a.dat',uxyt) -- Save t dependent data

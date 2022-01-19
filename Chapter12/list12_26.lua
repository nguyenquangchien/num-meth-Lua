-- File list12_26.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"
Vm = 5.0 -- Define equations to be solved
pi = math.pi;pi2 = 2*pi^2; sin=math.sin
feqs = { -- Table of functions
	function(x,y,uxx,ux,u,uy,uyy,i,j) -- General point
		return uxx + uyy + Vm*pi2*(sin(pi*x/xmax)*sin(pi*y/ymax))^20
	end,
	function(x,u,uy,ux,i) return u end, -- Bottom boundary -- Zero 
	function(x,u,uy,ux,i) return u end, -- Top boundary
	function(y,u,ux,uy,j) return u end, -- Right boundary
	function(y,u,ux,uy,j) return u end  -- Left boundary
}  -- End general point and boundary values

x,y,u = {},{},{}; xmax,ymax = 1.0,1.0 -- Try other values
Nx,Ny = 80,80
nx,ny = Nx+1,Ny+1; n = nx*ny
for i = 1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set x,y grid of initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end
a,b = setup2bveqs(feqs,x,y,u,1) -- Set up equations, x first numbering

ar,br = reorder(a,b) -- Reorder with y first numbering

wx,wy = 100,100 -- Set ADI factors -- Try other values
updatew(ar,wy); updatew(a,wx)

t = os.time()
ua = Spmat.new(nx,-ny) -- Array for solution
for j=1,ny do -- Linear array of solution values
	for i=1,nx do ua[i+(j-1)*nx] = u[j][i] end
end	
era,ier = {},1; nadi = 50
for k=1,nadi do -- Loop over ADI solutions
	ua,errm = trisolve(a,b,ua) -- X first labeling
	print('k, first max corr =',k,errm);io.flush()
	era[ier],ier = errm,ier+1
	ua = reorder(ua) -- X first to Y first 
	ua,errm = trisolve(ar,br,ua) --Y first labeling
	print('k, second max corr =',k,errm);io.flush()
	era[ier],ier = errm,ier+1
	ua = reorder(ua) -- Y first to X first
end
print('time = ',os.time()-t)
print(' At end maximum correction =',errm)
bb = to2darray(ua) -- Convert to 2D array for plots
splot(bb); cplot(bb) -- Surface and coutour plots
write_data('list12_26.dat',bb)
write_data('list12_26a.dat',era)
splot('list12_26.emf',bb)


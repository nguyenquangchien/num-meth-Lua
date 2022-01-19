-- File list12_26.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"; require"spgauss"
Vm = 5.0 -- Define equations to be solved
pi = math.pi;pi2 = 2*pi^2; sin=math.sin
feqs = { -- Table of functions
	function(x,y,uxx,ux,u,uy,uyy,i,j) -- General point
		return uxx + uyy + Vm*pi2*(sin(pi*x/xmax)*sin(pi*y/ymax))^20
	end,
	function(x,u,uy,ux,i) return u end, -- Bottom boundary
	function(x,u,uy,ux,i) return u end, -- Top boundary
	function(y,u,ux,uy,j) return u end, -- Right boundary
	function(y,u,ux,uy,j) return u end  -- Left boundary
}  -- End general point and boundary values

x,y,u = {},{},{}; xmax,ymax = 10,1
Nx,Ny = 100,100
nx,ny = Nx+1,Ny+1; n = nx*ny
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end
--a,b = setup2bveqs(feqs,x,y,u) -- Set up equations, x first numbering
--getfenv(spgauss).nprint=1
--ue = spgauss(a,b)

a,b = setup2bveqs(feqs,x,y,u,1) -- Set up equations, x first numbering

errmf =function(u)
	local erm = 0.0
	for m=1,n do erm = math.max(erm,math.abs(u[m]-ue[m])) end
	return erm
end

ar,br = reorder(a,b) -- Reorder with y first numbering

w = 1-- Set ADI factor
--wx = 100*w
wx = .83
wy = wx
updatew(ar,wy)
updatew(a,wx)

t = os.time()
ua = Spmat.new(nx,-ny) -- Array for solution
for j=1,ny do
	for i=1,nx do ua[i+(j-1)*nx] = u[j][i] end
end	
era,ier = {},1; nadi = 50
for k=1,nadi do
	ua,errm1 = trisolve(a,b,ua)
	--print('k, first max corr =',k,errm1);io.flush()
	--era[ier],ier = errmf(ua),ier+1
	era[ier],ier = errm1,ier+1
	ua = reorder(ua)
	ua,errm2 = trisolve(ar,br,ua)
	--print('k, second max corr =',k,errm);io.flush()
	era[ier],ier = errm2,ier+1
	ua = reorder(ua)
	--era[ier],ier = errmf(ua),ier+1
end
print('time = ',os.time()-t)
print(' At end maximum correction =',wx,errm2)
bb = to2darray(ua) -- Convert to 2D array
splot(bb); cplot(bb)
--write_data('list12_26.dat',bb)
--write_data('list12_26e.dat',era)


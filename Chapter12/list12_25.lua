-- File list12_25.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"; getfenv(pde2bvcoe).nprint=1; getfenv(pde2bvcoe).sverr=1
Vm = 1.0 -- Define equations to be solved
pi = math.pi;pi2 = 2*pi^2; sin=math.sin
feq = function(x,y,uxx,ux,u,uy,uyy,i,j)
	return uxx + uyy + Vm*pi2*sin(pi*x)*sin(pi*y)
end
fb = function(x,u,uy,ux,i) return u end
ft = function(x,u,uy,ux,i) return u end
fr = function(y,u,ux,uy,j) return u end
fl = function(y,u,ux,uy,j) return u end

x,y,u = {},{},{}; xmax,ymax = 1.0,1.0
Nx,Ny = 100,100
nx,ny = Nx+1,Ny+1; n = nx*ny
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end

a,b = setup2bveqs({feq,fb,ft,fl,fr},x,y,u)

t = os.time()
ua = {}
for j=1,ny do
	for i=1,nx do ua[i+(j-1)*nx] = u[j][i] end
end	

ua,errar = pde2bvcoe(a,b,ua)
print('time = ',os.time()-t)
bb = to2darray(ua,Nx,Ny); splot(bb)
write_data('list12_25.dat',bb)
write_data('list12_25_err.dat',errar)

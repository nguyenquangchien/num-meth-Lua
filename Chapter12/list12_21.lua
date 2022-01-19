-- File list12_21.lua --
-- Simple Example of BV problem in 2 dimensions -- rectangular grid 

require"pde2bv"; require"spgauss"
getfenv(spgauss).nprint=1; getfenv(spgauss).usage=2

Vm = 1.0 -- Define equations to be solved
feq = function(x,y,uxx,ux,u,uy,uyy,i,j)
	return uxx + uyy -- Poisson's equation
end
fb = function(x,u,uy,ux,i) return u end -- 0 at bottom
ft = function(x,u,uy,ux,i) return u - Vm end -- Vm at top
fr = function(y,u,ux,uy,j) return ux end -- zero derivatives
fl = function(y,u,ux,uy,j) return ux end

x,y,u = {},{},{}; xmax,ymax = 1.0,1.0
Nx,Ny = 40,40 -- 40 by 40 uniform grid
nx,ny = Nx+1,Ny+1
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}
	for i = 1,nx do u[j][i] = 0 end
end
-- Set up matrix equations
a,b = setup2bveqs({feq,fb,ft,fl,fr},x,y,u)

_,_,nel = spgauss(a,b) -- Collect usage statistics
print('max #matrix elements = ',nel[2][#nel[2]])
bb = to2darray(b) -- Convert to x,y array
write_data("list12_21.dat",bb)

udiff,mdiff,im,jm = {},0.0,0,0
udiff[0] = b[0]; Vmx = b[#b]
for j=1,ny do -- Calculate error in solution
	for i=1,nx do 
		m = i + (j-1)*nx
		udiff[m] = math.abs((Vmx*(j-1)/(ny-1)) - b[m])
		if udiff[m]>mdiff then mdiff,im,jm = udiff[m],i,j end
	end
end
print('max error = ',mdiff, '\nat j,i = ',jm,im)
sdsav = to2darray(udiff) -- 2D Array of errors
write_data('tmp12_21a.dat',sdsav) -- Save errors
write_data('tmp12_21b.dat',nel) -- Save fill statistics

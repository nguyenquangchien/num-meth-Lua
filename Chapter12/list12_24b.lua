-- File list12_25.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"
Vm = 1.0 -- Define equations to be solved

pi = math.pi;pi2 = 2*pi^2; sin=math.sin
feq = function(x,y,uxx,ux,u,uy,uyy,i,j)
	return uxx + uyy + Vm*pi2*sin(pi*x)*sin(pi*y)
	--return uxx + uyy 
end
--[[
fb = function(x,u,uy,ux,i) return u end
ft = function(x,u,uy,ux,i) return u - Vm end
fr = function(y,u,ux,uy,j) return ux end
fl = function(y,u,ux,uy,j) return ux end
--]]
fb = function(x,u,uy,ux,i) return u end
ft = function(x,u,uy,ux,i) return u end
fr = function(y,u,ux,uy,j) return u end
fl = function(y,u,ux,uy,j) return u end

x,y,u = {},{},{}; xmax,ymax = 1.0,1.0
Nx,Ny = 80,80
nx,ny = Nx+1,Ny+1; n = nx*ny
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end
function maxerr(ua)
	udiff,mdiff,im,jm = {},0.0,0,0
	Vmx = ua[#ua]
	for j=1,ny do
		yy = ymax*(j-1)/Ny
		for i=1,nx do 
			xx = xmax*(i-1)/Nx
			m = i + (j-1)*nx
			udiff[m] = math.abs(Vm*sin(pi*xx)*sin(pi*yy) - ua[m])
			if udiff[m]>mdiff then mdiff,im,jm = udiff[m],i,j end
		end
	end
	return udiff,mdiff
end

a,b = setup2bveqs({feq,fb,ft,fl,fr},x,y,u)

t = os.time()
ua = {}
for j=1,ny do
	for i=1,nx do ua[i+(j-1)*nx] = u[j][i] end
end	

rsp = ((x[2]-x[1])/(y[2]-y[1]))^2
rsp = ((math.cos(math.pi/nx)+rsp*math.cos(math.pi/ny))/(1+rsp))
--rsp = ((math.cos(2*math.pi/nx)+rsp*math.cos(2*math.pi/ny))/(1+rsp))
rsp = math.min(Nx,Ny)
rsp = math.cos(math.pi/rsp)
--rsp = rsp/1.002


print('rsp = ',rsp)
print('theoritical lan = ',2/(1+math.sqrt(1-rsp)))

ua,errar = pde2bvcoe(a,b,ua,rsp); ua[0] = {nx,ny,1}
--a,b = rev2deqs(a,b)
--ua,errar = pde2bvsor(a,b,ua,rsp)
--ua = rev2deqs(ua,ny,nx)
udiff,mdiff = maxerr(ua); udiff[0] = {nx,ny,1}
print('Final error = ',mdiff)

print('time = ',os.time()-t)
--bb = to2darray(x,y,ua)
bb = to2darray(ua)
write_data('list12.24a.dat',bb)
write_data('list12.24a_err.dat',errar)
--write_data('list12.24_60_err.dat',itarr,itarr1,errarr)
splot(bb)


sdsav = to2darray(udiff)
write_data('list12_24a_diff.dat',sdsav)

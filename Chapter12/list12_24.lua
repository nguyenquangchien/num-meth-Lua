-- File list12_24.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"
Vm = 1.0 -- Define equations to be solved
pi = math.pi;pi2 = 2*pi^2; sin=math.sin
feq = function(x,y,uxx,ux,u,uy,uyy,i,j)
	return uxx + uyy + Vm*pi2*sin(pi*x)*sin(pi*y)
end -- Now Zero boundary values
fb = function(x,u,uy,ux,i) return u end
ft = function(x,u,uy,ux,i) return u end
fr = function(y,u,ux,uy,j) return u end
fl = function(y,u,ux,uy,j) return u end

x,y,u = {},{},{}; xmax,ymax = 1.0,1.0
Nx,Ny = 50,50
nx,ny = Nx+1,Ny+1; n = nx*ny
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; 	for i = 1,nx do u[j][i] = 0 end
end

a,b = setup2bveqs({feq,fb,ft,fl,fr},x,y,u)

t = os.time()
ua = {}
for j=1,ny do
	for i=1,nx do ua[i+(j-1)*nx] = u[j][i] end
end	

rsp = ((x[2]-x[1])/(y[2]-y[1]))^2 -- Needed for lan
rsp = ((math.cos(math.pi/nx)+rsp*math.cos(math.pi/ny))/(1+rsp))
p4 = rsp^2/4; jprint=0; nprint = math.min(Nx,Ny)/4
print('theoritical lan = ',2/(1+math.sqrt(1-rsp^2)))

function sdloop(a,b,ua,itt)
	local uold,jpold,lan = {},0,1
	lan = 1.0
	for k=1,itt do
		ua,merr = sdsolve(a,b,ua,lan)
		jprint = math.floor(k/nprint)
		if jprint==jpold then
			jpold=jprint+1
			print("Completed iteration",k,"with correction",merr,lan); io.flush()
		end
		lan = 1/(1- lan*p4) -- Update lan
	end
	return ua
end	

ua = sdloop(a,b,ua,250)
print('time = ',os.time()-t)
bb = to2darray(ua,Nx,Ny)
splot(bb); write_data('list12_24.dat',bb)

-- File list12_23.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"
Vm = 1.0 -- Define equations to be solved
feq = function(x,y,uxx,ux,u,uy,uyy,i,j)
	return uxx + uyy
end
fb = function(x,u,uy,ux,i) return u end
ft = function(x,u,uy,ux,i) return u - Vm  end
fr = function(y,u,ux,uy,j) return ux end
fl = function(y,u,ux,uy,j) return ux end

x,y,u = {},{},{}; xmax,ymax = 1.0,1.0
Nx,Ny = 50,50
nx,ny = Nx+1,Ny+1
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}
	for i = 1,nx do u[j][i] = 0 end
end

a,b = setup2bveqs({feq,fb,ft,fl,fr},x,y,u) -- get matrix equations

ua,ue = {},{}
for j=1,ny do -- convert u to single array
	for i=1,nx do 
	m = i+(j-1)*nx
	ua[m],ue[m] = u[j][i],Vm*y[j] end
end	

function errins(u,ue)
	local err = 0
	for i=1,#u do err = math.max(math.abs(u[i]-ue[i]),err) end
	return err	
end

function sdloop(a,b,ua,itt) -- Iterative loop 
	local jpold = 1
	for k=1,itt do 
		ua,merr = sdsolve(a,b,ua) -- Single step
		jprint = math.floor(k/100)
		if jprint==jpold then
			jpold=jprint+1
			print("Completed iteration",k,"with correction",merr,errins(ua,ue));io.flush() 
		end
	end
	return ua
end	

ua = sdloop(a,b,ua,4000)
bb = to2darray(ua,Nx,Ny)
write_data('list12_23.dat',bb)
splot(bb)

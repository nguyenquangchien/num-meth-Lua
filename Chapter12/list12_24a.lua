-- File list12_24.lua --
-- Example of BV problem in 2 dimensions -- rectangular grid of points

require"pde2bv"
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
Nx,Ny = 80,80
nx,ny = Nx+1,Ny+1; n = nx*ny
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for i = 1,ny do y[i] = ymax*(i-1)/Ny end
for j = 1,ny do	-- Set zero initial values
	u[j] = {}
	for i = 1,nx do u[j][i] = 0 end
end

a,b = setup2bveqs({feq,fb,ft,fl,fr},x,y,u)

t = os.time()
ua,ub,ue = {}, {},{}
for j=1,ny do
	for i=1,nx do
		m = i+(j-1)*nx
		ua[m] = u[j][i]; ub[m] = ua[m]
		ue[m] = Vm*sin(pi*x[i])*sin(pi*y[j])
	end
end
n = #ua
function errins(u,ue)
	local err = 0
	for i=1,#u do err = math.max(math.abs(u[i]-ue[i]),err) end
	return err	
end

rsp = ((x[2]-x[1])/(y[2]-y[1]))^2 -- Needed for lan
--rsp = ((math.cos(math.pi/nx)+rsp*math.cos(math.pi/ny))/(1+rsp))^2
rsp = ((math.cos(math.pi/nx)+rsp*math.cos(math.pi/ny))/(1+rsp))
print('theoritical lan = ',2/(1+math.sqrt(1-rsp)))
lth = 2/(1+math.sqrt(1-rsp))
lth = 2/(1+math.pi/Nx)
nprint = math.min(Nx,Ny)/4;jprint=0

function sdloop(a,b,ua,itt)
	local uold,jpold,lan = {},0,1
	local p4 = rsp^2/4
	for i=1,#b do uold[i] = 0.0 end -- Array for k-1 values
	lan = 1.0
	for k=1,itt do
		ua,merr = sdsolve(a,b,ua)--,lan)
		for m=1,n do ua[m],ub[m] = lan*ua[m]+(1-lan)*ub[m],ua[m] end
		jprint = 0--math.floor(k/nprint)
		if jprint>=jpold then
			--jpold=jprint+1
			print("Completed iteration",k,"with correction",merr,lan,errins(ua,ue)); io.flush()
		end
		lan = 1/(1- lan*p4) -- Update lan
	end
	return ua
end	
function pde2bvcoe(a,b,u,rspin,umx) -- Odd/even with Chebyshev acceleration
	local rsp = rspin
	if rsp==nil then
		rsp = math.min(math.abs(a[0][1]),math.abs(a[0][2]))
		rsp = math.cos(math.pi/rsp)
	end
	local p4 = rsp^2/4
	local errst,uold = {{},{}}, {}
	local neq,vm,umx = #b, 0.0, umx or 0.0
	local lan1,lan2 = 1, 1/(1-2*p4) --math.min(1/(1-2*p4),1.5)
	local err = ERR/(2*math.sqrt(neq))
	iprint = math.floor(math.sqrt(neq)/4)
	for k=1,NMAX do
		u,mdiff,um = sdoesolve(a,b,u,lan1,lan2)
		if sverr==1 then errst[1][k],errst[2][k] = k, mdiff end
		if math.abs(mdiff)<err*math.max(umx,um) then 
			print('Exiting COE at iteration',k,'with correction',mdiff); io.flush(); break end 
		jprint = math.floor(k/iprint)+1
		lan1 = 1/(1-lan2*p4); lan2 = 1/(1-lan1*p4)
		--if jprint==nprint or iprint==0 then
			nprint=jprint+1
			print("Completed COE iteration ",k," with correction",mdiff,errins(u,ue)); io.flush() 
		--end
	end
	nprint = math.min(1,nprint)
	return u,errst
end
setfenv(pde2bvcoe,{ERR=1.e-14,NMAX=400,nerr=0,print=print,math=math,
table=table,io=io,sverr=1,nprint=0,sdoesolve=sdoesolve,errins=errins,ue=ue})
function pde2bvsor(a,b,u,rspin,umx) -- SOR relaxation
	local rsp = rspin
	if rsp==nil then 
		rsp = math.min(math.abs(a[0][1]),math.abs(a[0][2]))
		rsp = math.cos(math.pi/rsp)
	end
	local p4 = rsp^2/4
	local errst,uold = {{},{}}, {}
	local neq,lan,umx= #b, 1, umx or 0.0
	local err = ERR/(2*math.sqrt(neq))
	local iprint,jprint = math.floor(math.sqrt(neq)/4)
	for k=1,NMAX do
		u,mdiff,um = sdsolve(a,b,u,lan)
		if sverr==1 then errst[1][k],errst[2][k] = k, mdiff end
		if math.abs(mdiff)<err*math.max(umx,um) then 
			print('Exiting SOR at iteration',k,'with correction',mdiff); io.flush(); break end 
		jprint = math.floor(k/iprint)+1; lan = 1/(1-lan*p4)
		--if jprint==nprint then
			nprint=jprint+1
			print("Completed SOR iteration ",k," with correction",mdiff,errins(u,ue)); io.flush() 
		--end
	end
	nprint = math.min(1,nprint)
	return u,errst
end
setfenv(pde2bvsor,{ERR=1.e-14,NMAX=400,nerr=0,print=print,math=math,
table=table,io=io,sverr=1,nprint=0,sdsolve=sdsolve,errins=errins,ue=ue})
print('n =',n); getfenv(pde2bvsor).nprint=1
--ua = sdloop(a,b,ua,200)
ua = pde2bvcoe(a,b,ua,rsp)
print('time = ',os.time()-t)
bb = to2darray(ua,Nx,Ny)
splot(bb)
write_data('list12_24a.dat',bb)

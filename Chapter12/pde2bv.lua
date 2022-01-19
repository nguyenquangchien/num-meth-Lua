-- File pde2bv.lua --
-- Code for PDE BV problems in 2 dimensions -- rectangular grid of points
require'sdgauss'
setup2bveqs = function(eqs,x,y,u,ndg) -- Set up 2D equation set
	-- u is assumed to have the form u[j][i], j -> y grid and i -> x grid
	local u = u
	local uxx,ux,uij,uy,uyy,fxx,mxx,vx,xy
	local nx,ny,ndx,ndy = #x, #y, 0, 0
	local fctuxx,fctux,fctu,fctuy,fctuyy = FACT,FACT,FACT,FACT,FACT
	local a, b = Spmat.new(nx,ny), Spmat.new(nx,-ny)
	local sx,sy,ffxx,ffyy = {},{},{0,0,0}, {0,0,0}
	local alfxi,alfyi,fx1,fx2,fx3,fy1,fy2,fy3
	if #u~=nx*ny then -- for u in x,y format
		local ux = {} 
		for j=1,ny do for i=1,nx do m = i + (j-1)*nx; ux[m] = u[j][i] end end
		u = ux 
	end
	for i=1,nx-1 do sx[i] = x[i+1] - x[i] end
	for j=1,ny-1 do sy[j] = y[j+1] - y[j] end
	local function setx(fxc,ffxc) -- Store x oriented elements
		if fxc~=0.0 then for k=1,3 do
			mxx = mx+k-2
			if mxx==m and ndg~=nil then am[mxx][1] = am[mxx][1] + fxc*ffxc[k] 
			else am[mxx] = (am[mxx] or 0) + fxc*ffxc[k] end
		end end
	end
	local function sety(fyc,ffyc) -- Store y oriented elements
		if fyc~=0.0 then for k=1,3 do
			mxx = my +(k-2)*nx
			if mxx==m and ndg~=nil then am[mxx][2] = am[mxx][2] + fyc*ffyc[k]
			else am[mxx] = (am[mxx] or 0) + fyc*ffyc[k] end
		end end
	end
	local function bound(eq,x,uji,ux,uy,ij,nfeq)
		if nfeq<4 then 
			fv = eq(x,uji,uy,ux,ij)
			fu1 = (eq(x,uji,uy,ux+fctux,ij)-fv)/fctux
			fu2 = (eq(x,uji,uy+fctuy,ux,ij)-fv)/fctuy
			fu = (eq(x,uji+fctu,uy,ux,ij)-fv)/fctu
		else
			fv = eq(x,uji,ux,uy,ij)
			fu1 = (eq(x,uji,ux+fctux,uy,ij)-fv)/fctux
			fu2 = (eq(x,uji,ux,uy+fctuy,ij)-fv)/fctuy
			fu = (eq(x,uji+fctu,ux,uy,ij)-fv)/fctu
		end
		return fv, 0.0, fu1, fu, fu2, 0.0 
	end
	local function getfac(i,j,nfeq)
		eq,xi,yj = eqs[nfeq], x[i], y[j]; uji = u[i+(j-1)*nx] 
		if nfeq==1 then 
			fv = eq(xi,yj,uxx,ux,uji,uy,uyy,i,j)
			return fv, (eq(xi,yj,uxx+fctuxx,ux,uji,uy,uyy,i,j)-fv)/fctuxx,
			(eq(xi,yj,uxx,ux+fctux,uji,uy,uyy,i,j)-fv)/fctux, 
			(eq(xi,yj,uxx,ux,uji+fctu,uy,uyy,i,j)-fv)/fctu,
			(eq(xi,yj,uxx,ux,uji,uy+fctuy,uyy,i,j)-fv)/fctuy, 
			(eq(xi,yj,uxx,ux,uji,uy,uyy+fctuyy,i,j)-fv)/fctuyy
		else
			if nfeq<4 then return bound(eq,xi,uji,ux,uy,i,nfeq)
			else return bound(eq,yj,uji,ux,uy,j,nfeq) end
		end
	end
	mxx,vx,vy = 0.0, nx/abs(x[nx]-x[1]), ny/abs(y[ny]-y[1])
	for j=1,ny do for i=1,nx do -- Find maximum solution value
		mxx = max(mxx,abs(u[i+(j-1)*nx])) 
	end end
	if mxx~=0.0 then -- Scale probe factors to maximum solution value
		fctu = FACT*mxx
		fctuxx,fctux,fctuy,fctuyy = fctu*vx^2, fctu*vx, fctu*vy, fctu*vy^2
	end
	for j=1,ny do -- Loop over y-dimensions
		for i=1,nx do -- Loop over x-dimensions
			if j==1 then -- Bottom boundary row
				alfyi = sy[2]/sy[1]; fxx = 1/(alfyi*(sy[2]+sy[1]))
				fy1,fy2,fy3,ndy,nfeq = -alfyi*(2+alfyi)*fxx,(1+alfyi)^2*fxx,-fxx,1,2
			elseif j==ny then -- Top boundary
				alfyi = sy[ny-1]/sy[ny-2]; fxx = 1/(alfyi*(sy[ny-1]+sy[ny-2]))
				fy1,fy2,fy3,ndy,nfeq = alfyi^2*fxx,-(1+alfyi)^2*fxx,(1+2*alfyi)*fxx,-1,3
			else -- General interior point
				alfyi = sy[j]/sy[j-1]; fxx = 1/(alfyi*(sy[j]+sy[j-1]))
				fy1,fy2,fy3,ndy = -alfyi^2*fxx,(alfyi^2-1)*fxx,fxx,0
				fxx = 2/(sy[j]*(sy[j]+sy[j-1]))
				ffyy = {alfyi*fxx,-(alfyi+1)*fxx,fxx}
				if j>1 and j<ny then nfeq = 1 end
			end
			if i==1 then -- Left boundary
				alfxi = sx[2]/sx[1]; fxx = 1/(alfxi*(sx[2]+sx[1]))
				fx1,fx2,fx3,ndx = -alfxi*(2+alfxi)*fxx,(1+alfxi)^2*fxx,-fxx,1
				if j>1 and j<ny then nfeq = 4 end
			elseif i==nx then -- Right boundary
				alfxi = sx[nx-1]/sx[nx-2]; fxx = 1/(alfxi*(sx[nx-1]+sx[nx-2]))
				fx1,fx2,fx3,ndx = alfxi^2*fxx,-(1+alfxi)^2*fxx,(1+2*alfxi)*fxx,-1
				if j>1 and j<ny then nfeq = 5 end
			else -- General interior point
				alfxi = sx[i]/sx[i-1]; fxx = 1/(alfxi*(sx[i]+sx[i-1]))
				fx1,fx2,fx3,ndx = -alfxi^2*fxx,(alfxi^2-1)*fxx,fxx,0
				fxx = 2/(sx[i]*(sx[i]+sx[i-1]))
				ffxx = {alfxi*fxx,-(alfxi+1)*fxx,fxx}
				if j>1 and j<ny then nfeq = 1 end
			end
			-- Now evaluate derivatives
			if j==1 or j==ny then jj,ii = j+ndy, i
			elseif i==1 or i==nx then jj,ii = j, i+ndx
			else jj,ii = j, i end
			m = ii + (jj-1)*nx
			ujim,uji,ujip,ujmi,ujpi = u[m-1], u[m], u[m+1], u[m-nx], u[m+nx]
			ux,uxx = fx1*ujim + fx2*uji + fx3*ujip, ffxx[1]*ujim + ffxx[2]*uji + ffxx[3]*ujip
			uy,uyy = fy1*ujmi + fy2*uji + fy3*ujpi, ffyy[1]*ujmi + ffyy[2]*uji + ffyy[3]*ujpi 
			-- Now probe equations for dependences
			fv, fxx,fx,fu,fy,fyy = getfac(i,j,nfeq)
			m = i + (j-1)*nx; b[m] = -fv -- Equation number -- diagonal is at m,m
			am = a[m] -- work on mth row
			if ndg==nil then am[m] = 0 else am[m] = {0,0,0} end
			if fu~=0.0 then -- Now store diagonal elements
				if ndg==nil then am[m] = fu 
				else am[m][3] = fu end  
			end
			mx, my  = m+ndx, m+ndy*nx
			ffx, ffy = {fx1,fx2,fx3}, {fy1,fy2,fy3}
			setx(fxx,ffxx); setx(fx,ffx) -- Store x oriented elements
			sety(fy,ffy); sety(fyy,ffyy) -- Store y oriented elements
		end
	end
	return a,b -- Return matrix elements
end
setfenv(setup2bveqs,{table=table,Spmat=Spmat,abs=math.abs,max=math.max,FACT=1.e-4,
print=print})

function to2darray(u,Nx,Ny) -- Convert to s-y array from linear array
	local nx,ny,ixy
	if u[0]~=nil then nx,ny,ixy = u[0][1], u[0][2], u[0][3] 
	else nx,ny,ixy = Nx+1,Ny+1,1 end
	if ixy==2 then nx,ny=ny,nx end
	ut = {}
	for j=1,ny do
		ut[j] = {}
		for i=1,nx do m=i+(j-1)*nx; ut[j][i] = u[m] end
	end
	return ut
end
setfenv(to2darray,{})

function reversexy(u)
	local nx,ny,ixy = u[0][1],u[0][2],u[0][3]
	if ixy==2 then nx,ny = ny,nx end
	local ut = {}; ut[0] = {nx,ny,3-u[0][3]}
	for i=1,nx do
		ut[i] = {}; for j=1,ny do ut[i][j] = u[j][i] end
	end
	return ut
end
setfenv(reversexy,{})
	
updatew = function(a,w1,w2)
	local nx,ny,ixy = a[0][1],a[0][2],a[0][3]
	w2 = w2 or 0.0; w1 = w1-w2
	if ixy==2 then nx,ny = ny,nx end
	for j=1,ny do
		for i=1,nx do
			m = i + (j-1)*nx; v = a[m][m]
			v[1],v[2] = v[1]-w1,v[2]+w1
		end
	end
	return a
end
setfenv(updatew,{})

reorder = function(...) --function(a,b,nx,ny)
	local arg,af,nf,f,nyy,an = {...},{}
	local nm,lan = #arg, 0.0
	local v1,v2,nx,ny,ixy
	if type(arg[nm])=='number' then nm,lan = nm-1,arg[nm] end
	for nf=1,nm do
		local a = arg[nf]
		if type(a)=="SPMAT" then nx,ny,ixy = a[0][1],a[0][2],a[0][3] 
		else
			print("Size of martix not known in reorder"); return a
		end
		if type(a[1])=='table' then 
			na,ahere,nyy = #a, true, ny
			if na~=nx*ny then print("Incompatable sizes of u and nx,ny in reorder") end
		else ahere,nyy = false,-ny end
		an = Spmat.new(nx,nyy,3-ixy) -- New reversed table
		if ixy==2 then nx,ny = ny,nx end
		for j=1,ny do
			for i=1,nx do
				m = i + (j-1)*nx -- row location of a
				mn = j + (i-1)*ny -- row location in new matrix
				if ahere then 
					row = a[m]
					for k,v in pairs(row) do -- now get elements along row
						if k==m then -- Reverse order of elements
							an[mn][mn] = {}; rel = an[mn][mn]
							rel[ixy],rel[3-ixy],rel[3] = v[3-ixy],v[ixy],v[3]
						elseif abs(k-m)<nx then an[mn][mn+ny*(k-m)] = v
						else -- Need to insure integer values
							if k<m then sgn=-1 else sgn=1 end
							kk = (abs(k-m)+1)/nx
							kk = floor(kk)*sgn; an[mn][mn+kk] = v
						end
					end
				else an[mn] = a[m] end
			end
		end
		if ahere and lan~=0.0 then an = updatew(an,lan) end
		af[nf] = an
	end
	return unpack(af)
end
setfenv(reorder,{Spmat=Spmat,print=print,type=type,abs=math.abs,floor=math.floor,table=table,
unpack=unpack,updatew=updatew,pairs=pairs})

sdoesolve = function(a,b,u,lan1in,lan2in) -- Single step of Odd/even solving
	local n,epsm,an,row,eps,jm,im = #b,0.0
	local um = 0.0
	local lan1,lan2 = lan1in or 1,lan2in or lan1in or 1
	lan = {lan1,lan2}
	for ipass = 1,2 do -- Step over even and odd solutions
		for j=ipass,n,2  do -- Update odd or even values
			row,eps = a[j],-b[j]
			for i,v in pairs(row) do
				if i==j then an = v end
				eps = eps + v*u[i] 
			end
			eps = eps/an; u[j] = u[j] - lan[ipass]*eps -- lan1 or lan2 over-relaxation
			if abs(eps)>abs(epsm) then epsm,jm,im = eps,j,i end
			eps = abs(u[j]); if eps>um then um = eps end 
		end
	end
	return u, epsm, um, jm, im
end
setfenv(sdoesolve,{table=table,type=type,abs=math.abs,pairs=pairs})

sdsolve = function(a,b,u,lan) -- Single step of SOR
	local n,epsm,bb,row,eps,jm,im = #b,0.0
	local um = 0.0
	lan = lan or 1
	for j=1,n do -- Update each row in turn
		row,eps = a[j],-b[j]
		for i,v in pairs(row) do eps = eps + v*u[i] end
		eps = eps/a[j][j]; u[j] = u[j] - lan*eps
		if abs(eps)>abs(epsm) then epsm,jm,im = eps,j,i end
		um = max(um,abs(u[j]))
	end
	return u, epsm, um, jm, im
end
setfenv(sdsolve,{table=table,type=type,abs=math.abs,max=math.max,pairs=pairs})

function pde2bvsor(a,b,u,rspin,umx) -- SOR solution method
	local rsp = rspin
	if rsp==nil then 
		rsp = min(abs(a[0][1]),abs(a[0][2])) -- min of nx,ny
		rsp = cos(pi/rsp)
	end
	if nprint~=0 then print('SOR spectral radius =',rsp) end
	local p4 = rsp^2/4 -- Estimated best value
	local errst,uold = {{},{}}, {}
	local neq,lan,umx= #b, 1, umx or 0.0
	local err = ERR/(2*sqrt(neq))
	local iprint = floor((sqrt(neq)-1)/4)
	for k=1,NMAX do
		u,mdiff,um = sdsolve(a,b,u,lan) -- Single step of SOR solution
		if sverr==1 then errst[1][k],errst[2][k] = k, mdiff end
		if abs(mdiff)<err*max(umx,um) then 
			if nprint~=0 then print('Exiting SOR at iteration',k,'with correction',mdiff); io.flush() end
			break
		end 
		jprint = floor(k/iprint)+1; lan = 1/(1-lan*p4)
		if jprint==nprint then
			nprint=jprint+1
			print("Completed SOR iteration ",k," with correction",mdiff); io.flush() 
		end
	end
	nprint = min(1,nprint)
	return u,errst
end
setfenv(pde2bvsor,{ERR=1.e-4,NMAX=4000,nerr=0,print=print,abs=math.abs,min=math.min,pi=math.pi,
cos=math.cos,sqrt=math.sqrt,floor=math.floor,max=math.max,table=table,io=io,sverr=0,nprint=0,sdsolve=sdsolve})

function pde2bvcoe(a,b,u,rspin,umx) -- Odd/even with Chebyshev acceleration method
	local rsp = rspin
	if rsp==nil then
		rsp = min(abs(a[0][1]),abs(a[0][2]))
		rsp = cos(pi/rsp)
	end
	if nprint~=0 then print('COE spectral radius =',rsp) end
	local p4 = rsp^2/4
	local errst,uold = {{},{}}, {}
	local neq,vm,umx = #b, 0.0, umx or 0.0
	local lan1,lan2 = 1, 1/(1-2*p4) 
	local err = ERR/(2*sqrt(neq))
	iprint = floor(sqrt(neq)/4)
	for k=1,NMAX do
		u,mdiff,um = sdoesolve(a,b,u,lan1,lan2) -- One step in solution 
		if sverr==1 then errst[1][k],errst[2][k] = k, mdiff end
		if abs(mdiff)<err*max(umx,um) then 
			if nprint~=0 then print('Exiting COE at iteration',k,'with correction',mdiff); io.flush() end
			break
		end 
		jprint = floor(k/iprint)+1
		lan1 = 1/(1-lan2*p4); lan2 = 1/(1-lan1*p4)
		if jprint==nprint then
			nprint=jprint+1
			print("Completed COE iteration ",k," with correction",mdiff); io.flush() 
		end
	end
	nprint = min(1,nprint)
	return u,errst
end
setfenv(pde2bvcoe,{ERR=1.e-4,NMAX=4000,nerr=0,print=print,abs=math.abs,min=math.min,cos=math.cos,
pi=math.pi,sqrt=math.sqrt,floor=math.floor,max=math.max,table=table,io=io,sverr=0,nprint=0,sdoesolve=sdoesolve})

function trisolve(a,b,u,umx) -- Solve tri-diagonal equations
	local nx,ny,ixy = a[0][1],a[0][2],a[0][3]
	local a1,b1,c1,d1,c1p,a1p
	local a2,b2,c2,d2,k
	local um,errm,err,row = 0.0,0.0
	local e,f = {0},{0}
	local un = Spmat.new(nx,-ny,ixy) -- Storage for new solution
	if ixy==2 then nx,ny = ny,nx end
	for j=2,ny-1 do -- step over rows, skipping first and last
		for i=1,nx do -- Solve row by row
			k = i + (j-1)*nx -- Node number
			row = a[k] -- Row of matrix
			if i==1 then -- Left boundary nodes are special cases
				a1,b1,c1,c1p,d1 = 0.0,0.0,0.0,0.0,-b[k]
				for m,v in pairs(row) do
					if m==k then b1,d1 = v[1]+v[3], d1+v[2]*u[k]
					elseif m==k+1 then c1 = v 
					elseif m==k+2 then c1p = v
					else d1 = d1 + v*u[m] end 
				end
			elseif i==2 then -- Now get first g,h value
				a2,b2,c2,d2 = 0.0,0.0,0.0,-b[k]
				for m,v in pairs(row) do
					if m==k-1 then a2 = v
					elseif m==k then b2,d2 = v[1]+v[3],d2+v[2]*u[k]
					elseif m==k+1 then c2 = v
					else d2 = d2+v*u[m] end 
				end
				a1 = 1/(b1*b2 - a2*c1)
				e[3],f[3] = a1*(a2*c1p - b1*c2), a1*(a2*d1 - b1*d2)
				a1,a1p = 1/(b1+a2), (c2+c1p)/e[3]
				e[2],f[2] = -a1*(a1p+b2+c1), a1*(a1p*f[3]-d1-d2)
			elseif i==nx then -- Right boundary nodes are special cases
				a2p,a2,b2,c2,d2 = 0.0,0.0,0.0,0.0,-b[k]
				for m,v in pairs(row) do
					if m==k-2 then a2p = v
					elseif m==k-1 then a2 = v
					elseif m==k then b2,d2 = v[1]+v[3],d2+v[2]*u[k] 
					else d2 = d2 + v*u[m] end
				end
				c2 = -1/( b2+a2*e[nx]+a2p*e[nx]*e[nx-1])
				un[k] = c2*(d2+a2*f[nx]+a2p*(e[nx-1]*f[nx]+f[nx-1]))
			else -- Now have internal node along row
				a1,b1,c1,d1 = 0.0,0.0,0.0,-b[k]
				for m,v in pairs(row) do
					if m==k-1 then a1 = v
					elseif m==k then b1,d1 = v[1]+v[3], d1+v[2]*u[k]
					elseif m==k+1 then c1 = v
					else d1 = d1+ v*u[m] end
				end
				a2 = -1/(a1*e[i] +b1)
				e[i+1],f[i+1] = a2*c1, a2*(d1+a1*f[i])
			end
		end
		-- Now have last 2 node voltages, back substiture 
		for i=nx-1,1,-1 do
			k = i + (j-1)*nx; un[k] = e[i+1]*un[k+1] + f[i+1] 
		end
	end
	a2,b2,c2 = floor((nx+1)/2),nx,1 -- Middle outward
	for md=1,2 do
		for j=1,ny,ny-1 do -- Update lower and upper boundaries
			for i=a2,b2,c2 do -- now update row of boundary values
				k = i + (j-1)*nx
				row = a[k]; d1 = -b[k]
				for m,v in pairs(row) do 
					if m==k then b1 = v[1] + v[2] + v[3] 
					else d1 = d1 + v*un[m] end
				end
				un[k] = -d1/b1
			end
		end
		a2,b2,c2 = a2-1,1,-1-- Other end
	end
	for k=1,nx*ny do -- Evaluate max correction and use lan factor
		um = max(um,abs(un[k]))
		ep = un[k] - u[k]; err = abs(ep)
		if err>errm then errm = err end
		un[k] = u[k] + ep
	end
	if umx~=nil then if errm>umx then -- Limit change to umx, if specified
		umx,errm = umx/errm,umx
		for k=1,nx*ny do un[k] = u[k] + umx*(un[k]-u[k]) end
	end end
	return un,errm,um
end
setfenv(trisolve,{print=print,floor=math.floor,max=math.max,abs=math.abs,pairs=pairs,
Spmat=Spmat,type=type})

pde2bvadi = function(a,b,u,wxin,wyin,umx) -- Solution using ADI method
	local wx,wy,wxx,wyy,ixy = wxin,wyin or wxin,0,0
	local n,nx,ny,errst,erm1,erm2 = #b, a[0][1], a[0][2],{{},{}}
	local err,umx = ERR/10,umx or 0
	wx,wy = wx or WXX, wy or WXX -- WXX from experience -- set to 87 !
	if wx<0.0 then -- Attempt to scale for xmax, ymax
		wx = min(wx,wy) -- Use smallest interval
		wx = WXX/wx^2; wy = wx -- From experience
	end	
	for m=1,n do -- Check sign of diagonal matrix terms
		wyy = wyy + a[m][m][1]
	end
	if wyy>0 then wx,wy = -wx, -wy end -- Exchange if negative diagonal elements
	if nprint~=0 then print('ADI wx,wy =',wx,wy) end
	local iprint = floor(sqrt(n)/10)
	local ar,br = reorder(a,b) -- Reverse order of rows and columns
	updatew(a,wx); updatew(ar,wy)
	for k=1,NADI do
		u,erm1 = trisolve(a,b,u) -- Solve for x first ordering
		u = reorder(u)
		u,erm2,um = trisolve(ar,br,u) -- Solve for y first ordering
		u = reorder(u) -- One ADI cycle completed
		erm1 = max(erm1,erm2)
		if sverr==1 then errst[1][k],errst[2][k] = k, erm1 end
		if erm1<err*max(umx,um) then 
			if nprint~=0 then print('Exiting ADI at iteration',k,'with correction',erm1); io.flush() end
			break
		end 
		jprint = floor(k/iprint)+1
		if jprint==nprint then
			nprint=jprint+1
			print("Completed ADI iteration ",k," with correction",erm1); io.flush() 
		end
	end
	updatew(a,-wx); nprint = min(1,nprint) -- Restore a to original value
	return u,errst
end
setfenv(pde2bvadi,{print=print,min=math.min,sqrt=math.sqrt,floor=math.floor,max=math.max,
trisolve=trisolve,reorder=reorder,table=table,io=io,updatew=updatew,ERR=5.e-5,
NADI=500,WXX=250,sverr=1,nprint=0,to2darray=to2darray})
	
pde2bv = function(feqs,x,y,u,tpsola) -- 2 spatial dimensions PDE Solver
	if type(tpsola)~='table' then tpsola = {tpsola} end
	local tpsol,rx,ry = tpsola[1] or SPM, tpsola[2], tpsola[3] or tpsola[2]-- Default to sparse matrix solution
	local umx,errm,a,b,n,uxy,ndg = 0.0
	local nx,ny = #x, #y
	local uold,ua,n,itmx = {},{},nx*ny,0
	local ur = Spmat.new(nx,-ny)
	if linear==1 then nnmx=1 else nnmx=NMX end -- One Newton cycle for linear eqns
	if tpsol==ADI then ndg = 1 end
	uold[0] = {u[0][1],u[0][2],u[0][3]}
	ua[0] = {u[0][1],u[0][2],u[0][3]}
	if #u==n then uxy=false else uxy=true end
	for j=1,ny do
		for i=1,nx do
			k = i + (j-1)*nx
			if uxy then ua[k] = u[j][i] else ua[k] = u[k] end
		end
	end
	for k=1,n do uold[k],umx,ua[k] = ua[k], max(umx,ua[k]), 0.0 end
	for int=1,nnmx do -- Newton iterative loop
		a,b = setup2bveqs(feqs,x,y,uold,ndg)
		if tpsol==SPM then -- Solve with sdgauss, sparse matrix solver
			sdgauss(a,b); ua = b -- b is new solution
		elseif tpsol==COE then -- Solve with Chebychev odd/even SOR
			ua = pde2bvcoe(a,b,ua,rx,umx)
		elseif tpsol==SOR then -- Solve with SOR
			ua = pde2bvsor(a,b,ua,rx,umx)
		elseif tpsol==ADI then -- Solve with ADI
			if rx==nil then rx = -abs(x[nx]-x[1]) end
			if ry==nil then ry = -abs(y[ny]-y[1]) end
			ua = pde2bvadi(a,b,ua,rx,ry,umx)
		else print('Unknown type solution request:',tpsol,' in pde2bv')
		end
		errm,umx,itmx = 0.0,0.0,itmx+1
		for k=1,n do errm,uold[k] = max(errm,abs(ua[k])), uold[k]+ua[k] end
		for k=1,n do umx,ua[k] = max(umx,abs(uold[k])), 0.0 end
		if nprint~=0 then print('Completed Newton iteration',int,'with correction',errm); io.flush()
			if seeplot~=0 then if seeplot==1 then splot(to2darray(uold)) else cplot(to2darray(uold)) end end end
		if errm<ERR*umx then itmx = int; break end
	end
	if itmx==NMX then print('Maximum number of iterations exceeded in pde2bv!!')
		io.flush() end
	if uxy==false then return uold,errm,itmx
	else
		for j=1,ny do ur[j] = {}; for i=1,nx do k = i+(j-1)*nx; ur[j][i] = uold[k] end end
		return ur,errm,itmx
	end
end
setfenv(pde2bv,{sdgauss=sdgauss,max=math.max,abs=math.abs,table=table,
setup2bveqs=setup2bveqs,pde2bvcoe=pde2bvcoe,pde2bvsor=pde2bvsor,pde2bvadi=pde2bvadi,
io=io,getfenv=getfenv,linear=0,SPM=1,COE=2,SOR=3,ADI=4,NMX=50,ERR=1.e-5,print=print,
type=type,nprint=0,io=io,Spmat=Spmat,seeplot=0,cplot=cplot,splot=splot,to2darray=to2darray})

pde1stp2bv1t = function(feqs,tvals,x,y,ua,tpsola,ut,utt) -- 2 spatial dimensions plus time PDE Solver
	local j, neq, t, h, h2,h2sq,hs,hx,hy,hz -- Local variables for function
	local unn,un,jfirst = {},{},0
	local nit,nitt,errm = 0,0 -- Number of iterations
	local nx,ny = #x, #y
	local neq,uxy = nx*ny
	local tmin,tmax,ntval = tvals[1],tvals[2],tvals[3]
	local u,ur = Spmat.new(nx,-ny), Spmat.new(nx,-ny) -- Single array with nx,ny information	
	ut = ut or {}; utt = utt or {}
	if #ua==neq then uxy=false else uxy=true end
	-- Functions to add next time values and time derivatives 
	local fpde,fbb,fbt,fbr,fbl = feqs[1], feqs[2], feqs[3], feqs[4], feqs[5]
	feq = { -- Local functions to add time and time derivatives
		function(x,y,uxx,ux,u,uy,uyy,i,j) -- General spatial point
			local k = i + (j-1)*nx
			local ut,utt = (u - un[k])/h2, (u - unn[k])/h2sq
			return fpde(x,y,t,uxx,ux,u,uy,uyy,ut,utt,i,j)
		end,
		function(x,u,uy,ux,i) -- Bottom boundary
			local ut,utt = (u - un[i])/h2, (u - unn[i])/h2sq
			return fbb(x,t,u,uy,ux,ut,utt,i)
		end,
		function(x,u,uy,ux,i) -- Top boundary
			local k = i+neq-nx
			local ut,utt = (u - un[k])/h2, (u - unn[k])/h2sq
			return fbt(x,t,u,uy,ux,ut,utt,i)
		end,
		function(y,u,ux,uy,j) -- Left boundary
			local k = 1+(j-1)*nx
			local ut,utt = (u - un[k])/h2, (u - unn[k])/h2sq
			return fbr(y,t,u,ux,uy,ut,utt,j)
		end,
		function(y,u,ux,uy,j) -- Right boundary
			local k = j*nx
			local ut,utt = (u - un[k])/h2, (u - unn[k])/h2sq
			return fbl(y,t,u,ux,uy,ut,utt,j)
		end
	}
	for j=1,ny do -- Local array for storing solution values in linear array
		for i=1,nx do
			k = i + (j-1)*nx; if uxy then u[k] = ua[j][i] else u[k] = ua[k] end
		end
	end
	t = tmin -- Initial t value
	hs = (tmax - t)/ntval -- Equal increments in t used, no adjusting step size
	-- If initial derivative not specified, use Backwards differencing for first 4 points
	if #ut~=neq then for m=1,neq do ut[m] = ut[m] or 0 end end -- Complete initial derivatives 
	if #utt~=neq then for m=1,neq do utt[m] = 0 end
		jfirst,h = 0,0.25*hs; h2,h2sq,hy,hx,hz = h,h*h,h,0,0 -- Set to BD parameters
	else 
		if bd~=false then jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = h,h*h,h,0,0 -- BD parameters
		else jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = hs/2,h*h/4,h,h/2,h*h/4 end -- TP parameters
	end
	for k=1,ntval do -- Main loop for incrementing independent variable (t)
		repeat -- Use backwards differencing for first interval with 4 sub intervals of size h/4
			jfirst = jfirst+1
			-- Set up yn, and ynn arrays and get ready to solve equations
			for m=1,neq do
				un[m] = u[m] + hx*ut[m] -- hx = 0 or h/2
				unn[m] = u[m] + hy*ut[m] + hz*utt[m] -- hy=h, hz=0 or (h/2)^2
				u[m] = u[m] + h*ut[m] -- Predicted value of u array
			end
			t = t + h -- Now increment t to next t value
			-- Calculate new u values at next time step, new values are returned in u array
			u,errm,nitt = pde2bv(feq,x,y,u,tpsola) -- Solve PDE at time t
			if nitt>nit then nit = nitt end -- Monitor maximun number of iterations
			-- New derivative values, same function as in fnext 
			for m=1,neq do ut[m],utt[m] = (u[m] - un[m])/h2,(u[m] - unn[m])/h2sq end
		until jfirst>=4 -- End of first interval repeat using Backwards difference
		if k==1 then 
		if bd~=false then jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = h,h*h,h,0,0 -- BD parameters
			else jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = hs/2,h*h/4,h,h/2,h*h/4 end -- TP parameters
		end   
		if nprint~=0 then print('Completed time =',t,' with correction',errm); io.flush() end
	end -- End of main loop on t, now return solution array
	if uxy==false then return u,errm,nit
	else
		for j=1,ny do ur[j] = {}; for i=1,nx do k = i+(j-1)*nx; ur[j][i] = u[k] end end
		return ur,errm,nit
	end
end -- End of pde1stp2bv1t
setfenv(pde1stp2bv1t,{table=table,pde2bv=pde2bv,print=print,Spmat=Spmat,io=io,
bd=false,nprint=0})

local usave, utsave
pde2bv1t = function(feqs,tvals,x,y,uar,tpsola) -- Multi Time Step Solver
	local ua,ub,uc,tl,ns = {},{},{}
	local ntp,ni 
	local NMAX,ND = getfenv(pde2bv).NMX,10,0
	local nx,ny = #x, #y
	local neq,j = nx*ny
	local nit,tmin,ntkeep,ntval,dtkeep = 0
	local u,ut,utt
	if #uar<nx then u,ut,utt = uar[1], uar[2] or {}, uar[3] or {}
	else u,ut,utt = uar, {}, {} end
	if #u==neq then uxy=false else uxy=true end
	for j=1,ny do -- Save as linear array of soultion values
		for i=1,nx do
			k = i+(j-1)*nx; if uxy then ub[k] = u[j][i] else ub[k] = u[k] end
		end
	end
	ua = usave(ua,ub,nx,ny) -- Save initial values
	tvals[3],tvals[4] = tvals[3] or NKEEP, tvals[4] or NTVAL
	if type(tvals[2])=='number' then tvals[2] = {tvals[2]} end
	ntvs =  #tvals[2]
	if type(tvals[3])~='table' then 
		t3,tvals[3] = tvals[3], {}; for i=1,ntvs do tvals[3][i] = t3 end
	end
	if type(tvals[4])~='table' then 
		t4,tvals[4] = tvals[4], {}; for i=1,ntvs do tvals[4][i] = t4 end
	end
	tmin = tvals[1]
	if type(tvals[2])=='number' then
		nkeep,ntval = tvals[3][1], tvals[4][1]
		ntps,tmax = nkeep*ntval, tvals[2]
		tvals[2],tvals[3] = {},{}
		for i=1,ntps do tvals[2][i],tvals[3][i] = tmin+(i-1)*(tmax-tmin)/(ntps-1), ntval end
	else
		tmin = tvals[1]
		local t2,t3,k = {tmin}, {tvals[4][1]}, 1 
		for i=1,ntvs do
			nsav,tmax = tvals[3][i] or NKEEP, tvals[2][i]
			for j=1,nsav do
				k = k+1
				t2[k] = tmin + j*(tmax - tmin)/(nsav)
				t3[k] = tvals[4][i] or NTVAL 
			end
			tmin = tmax
		end
		tvals[2],tvals[3] = t2, t3
	end
	ntkeep = max(#tvals[3]-1,1); tmin = tvals[2][1]
	for i=1,ntkeep do -- Step over time increments
		ub,errm,ni = pde1stp2bv1t(feqs,{tmin,tvals[2][i+1],tvals[3][i]},x,y,ub,tpsola,ut,utt) 
		tmin = tvals[2][i+1]
		if ni==NMAX then 
			print("Error: Maximum number of iterations exceeded in pde2bv")
			print("Results may not be accurate! Time parameters =",tl,tvals[2][i],ns);io.flush()
		end
		nit = max(nit,ni); ua = usave(ua,ub,nx,ny)
		if seeplot~=0 then print('Plot at time =',tmin); io.flush()
			cplot(reversexy(ua[i+1])) end
	end
	return ua,errm,nit
end
setfenv(pde2bv1t,{table=table,type=type,print=print,io=io,splot=splot,seeplot=0,max=math.max,
NKEEP=10,NTVAL=10,pde1stp2bv1t=pde1stp2bv1t,getfenv=getfenv,reversexy=reversexy,
cplot=cplot})

pde2bv1tqs = function(feqs,tvals,x,y,uar,tpsola) -- Quick scan log time solver
	local nps,npts,tpts -- save per decade, additional steps per save
	local ttvals,nl,nu,fact = {}
	local nt,i,j,neq = #tvals,2,0
	local ua,ub,uta,nit,ni = {},{},{},0 
	local nx,ny,u,ut,utt = #x, #y
	local neq = nx*ny
	local nprint = nprint or getfenv(pde2bv).nprint
	local nts,tmp
	if tvals[5]==nil then tpts = false 
	else
		tpts = true
		if type(tvals[5])=='number' then
			nts = tvals[5]; tvals[5] = {}; tmp = x[nx] - x[1]
			for i=1,nts+1 do tvals[5][i] = (i-1)*tmp/nts + x[1] end
		end	
	end
	local tsave = tvals[5]
	if #uar<ny then u,ut,utt = uar[1], uar[2] or {}, uar[3] or {}
	else u,ut,utt = uar, {}, {} end
	if nt<2 then print('Error, must specify two times in pd2vbv1tqs') return end
	if #u==neq then uxy=false else uxy=true end
	nps,npts = floor(tvals[3] or NPS),floor(tvals[4] or NPTS)
	for j=1,ny do -- Save as linear array of soultion values
		for i=1,nx do
			k = i+(j-1)*nx; if uxy then ub[k] = u[j][i] else ub[k] = u[k] end
		end
	end
	ua = usave(ua,ub,nx,ny) -- Save initial values
	nts,nl = nps*npts,tvals[1]; tmp = (tvals[2][1]-tvals[1])/nts
	nu = nl + tmp
	if tpts then uta = utsave(uta,ub,x,y,nl,tsave) end -- save initial values
	for i=1,nts do -- loop over initial interval
		ub,errm,ni = pde1stp2bv1t(feqs,{nl,nu,1,1},x,y,ub,tpsola,ut,utt)
		nit = max(nit,ni) 
		if tpts then uta = utsave(uta,ub,x,y,nu,tsave) end -- save time values
		if nprint~=0 then print('In quick scan, initial time =',nu);io.flush() end
		nl,nu = nu,nu+tmp
	end
	ua = usave(ua,ub,nx,ny) -- Save values after initial interval
	if seeplot~=0 then print('Plot at time =',nl); io.flush()
	if seeplot==1 then splot(reversexy(ua[2])) 
		elseif seeplot==2 then cplot(reversexy(ua[2])) end
	end
	nps,npts = floor(tvals[3] or NPS),floor(tvals[4] or NPTS)
	fact = 10^(1/(nps*npts)) -- Factor between steps
	nl = 10^(floor(log10(tvals[2][1])))
	nu = 10^(ceil(log10(tvals[2][2])))*1.000001/fact
	while nl<=nu do -- Step over log time spacings
		ub,errm,ni = pde1stp2bv1t(feqs,{nl,nl*fact,1},x,y,ub,tpsola,ut,utt)
		nit,j,nl = max(nit,ni), j+1, nl*fact
		if tpts then uta = utsave(uta,ub,x,y,nl,tsave) end -- save if requested
		if nprint~=0 then print('In quick scan, t =',nl); io.flush() end
		if j==npts then -- Time to save array
			i,j = i+1,0 -- Track saved number
			ua = usave(ua,ub,nx,ny); if seeplot==1 then splot(reversexy(ua[i])) 
				elseif seeplot==2 then cplot(reversexy(ua[i])) end
			if nprint~=0 then 
				print('Saved at time = ',nl,' Number of iterations in quick scan = ',nit,'\n')
				io.flush(); nit = 0
			end
		end
	end
	if tpts then return ua,uta,errm,nit else return ua,errm,nit end
end
setfenv(pde2bv1tqs,{print=print,type=type,floor=math.floor,ceil=math.ceil,max=math.max,
log10=math.log10,nprint,io=io,NPS=1,NPTS=10,getfenv=getfenv,seeplot=0,usave=usave,
io=io,pde1stp2bv1t=pde1stp2bv1t,to2darray=to2darray,splot=splot,cplot=cplot,
nprint=0,reversexy=reversexy})

grad = function(x,y,u,ord)
	local ux,uy,row = {},{}
	local nx,ny,uxy = u[0][1],u[0][2],u[0][3]
	ux[0],uy[0] = {nx,ny,uxy}, {nx,ny,uxy}
	local sx,sy = {},{}
	if nx~=#x or ny~=#y then 
		print('dimension of x or y does not match u in grad'); return ux,uy
	end
	if uxy==2 then nx,ny,x,y = ny,nx,y,x end
	for i=1,nx-1 do sx[i] = x[i+1] - x[i] end
	for j=1,ny-1 do sy[j] = y[j+1] - y[j] end
	if ord==2 then -- 2nd derivative calculation
		for j=1,ny do
			uy[j],row,uj = {},{},u[j]
			alf = sx[2]/sx[1]; row[1] = 2*(alf*uj[1] - (1+alf)*uj[2] + uj[3])/(sx[2]*(sx[1]+sx[2]))
			alf = sx[nx-1]/sx[nx-2]; row[nx] = 2*(alf*uj[nx-2] - (1+alf)*uj[nx-1] + uj[nx])/(sx[nx-1]*(sx[nx-2]+sx[nx-1]))
			for i=2,nx-1 do alf = sx[i]/sx[i-1]; row[i] = 2*(alf*uj[i-1] - (1+alf)*uj[i] + uj[i+1])/(sx[i-1]*(sx[i-1]+sx[i])) end
			ux[j] = row
		end
		for i=1,nx do
			alf = sy[2]/sy[1]; uy[1][i] = (alf*u[1][i] - (1+alf)*u[2][i] + u[3][i])/(sy[2]*(sy[2]+sy[1]))
			alf=sy[ny-1]/sy[ny-2]; uy[ny][i]=(alf*u[ny-2][i]-(1+alf)*u[ny-1][i]+u[ny][i])/(sy[ny-1]*(sy[ny-1]+sy[ny-2]))
			for j=2,ny-1 do	
			alf = sy[j]/sy[j-1]; uy[j][i] = (alf*u[j+1][i] - (1+alf)*u[j][i] + u[j-1][i])/(sy[j]*(sy[j]+sy[j-1])) end
		end
		if uxy==1 then return ux,uy else return uy,ux end -- return 2nd derivative
	end
	for j=1,ny do -- 1st derivative calculation 
		uy[j],row,uj = {}, {}, u[j];
		alf=sx[2]/sx[1]; row[1] = (-uj[1]*(2+alf)*alf + uj[2]*(1+alf)^2 - uj[3])/(alf*(sx[1]+sx[2]))
		alf=sx[nx-1]/sx[nx-2]; row[nx] = (uj[nx-2]*alf^2-uj[nx-1]*(1+alf)^2+uj[nx]*(1+2*alf))/(alf*(sx[nx-1]+sx[nx-2]))
		for i=2,nx-1 do alf = sx[i]/sx[i-1]; row[i] = (uj[i+1] + uj[i]*(alf^2-1) - alf^2*uj[i-1])/(alf*(sx[i]+sx[i-1])) end
		ux[j] = row
	end
	for i=1,nx	do
		alf = sy[2]/sy[1]; uy[1][i] = (-u[1][i]*(2+alf)*alf + u[2][i]*(1+alf)^2 - u[3][i])/(alf*(sy[2]+sy[1]))
		alf=sy[ny-1]/sy[ny-2]; uy[ny][i]=(u[ny-2][i]*alf^2-u[ny-1][i]*(1+alf)^2+u[ny][i]*(1+2*alf))/(alf*(sy[ny-1]+sy[ny-2]))
		for j=2,ny-1 do	
		alf = sy[j]/sy[j-1]; uy[j][i] = (u[j+1][i] + u[j][i]*(alf^2-1) - alf^2*u[j-1][i])/(alf*(sy[j]+sy[j-1])) end
	end	
	if uxy==1 then return ux,uy else return uy,ux end -- return 1st derivative
end
setfenv(grad,{table=table})
 
setxy = function(xv,yv,Nx,Ny)
	if Nx==nil then Nx = 100 end
	Ny = Ny or Nx
	local x,y,xymn,xymx = {}, {}, xv[1], xv[2]
	for i=1,Nx+1 do x[i] = xymn + (i-1)*(xymx-xymn)/Nx end
	xymn,xymx = yv[1],yv[2]
	for j=1,Ny+1 do y[j] = xymn + (j-1)*(xymx-xymn)/Ny end
	return x,y
end
setfenv(setxy,{})

intp2d = function(xy,u,xyp)
	local x,y,xp,yp = xy[1], xy[2], xyp[1], xyp[2]
	local u1,u2,u3,u4,k,dx,dy
	nx,ny,nu = #x, #y, #u
	if nu==ny then uxy = true else uxy = false end
	if y[jp]<yp then 
		while jp<ny do if y[jp]>=yp then break else jp = jp+1 end end
		jp = jp-1
	else 
		while jp>1 do if y[jp]<=yp then break else jp = jp-1 end end
	end
	if x[ip]<xp then
		while ip<nx do if x[ip]>=xp then break else ip = ip+1 end end
		ip = ip-1
	else 
		while ip>1 do if x[ip]<=xp then break else ip = ip-1 end end
	end
	if uxy then u1,u2,u3,u4 = u[jp][ip], u[jp][ip+1], u[jp+1][ip], u[jp+1][ip+1]
	else k = ip+(jp-1)*nx; u1,u2,u3,u4 = u[k], u[k+1], u[k+nx], u[k+nx+1] end
	xp,yp = (xp-x[ip])/(x[ip+1]-x[ip]), (yp-y[jp])/(y[jp+1]-y[jp])
	return u1 + (u2-u1)*xp + (u3-u1)*yp + (u1+u4-u2-u3)*xp*yp
end
setfenv(intp2d,{table=table,jp=1,ip=1,print=print})

usave = function(u1,u2,nx,ny)
	local nk,neq,up = #u1,#u2, {}
	for j=1,ny do -- Save as 2D array for return table
		up[j] = {}; for i=1,nx do k = i+(j-1)*nx; up[j][i] = u2[k] end
	end
	up[0] = {nx,ny,1}
	u1[nk+1] = up
	return u1
end
setfenv(usave,{table=table})
	
utsave = function(u1,u2,x,y,t,px)
	u1 = u1 or {}
	local nk,nt = #u1, #px
	if nk~=nt+1 then for i=1,nt+1 do u1[i] = {} end end
	nk = #u1[1]+1
	u1[1][nk] = t
	for i=1,nt do u1[i+1][nk] = intp2d({x,y},u2,px[i]) end
	return u1
end
setfenv(utsave,{table=table,intp2d=intp2d})
		
mnorm = function(a,b,n) -- Normalize lines of matrix equations -- probably not needed
	n = n or #b
	local td
	if type(a[1][1])=='table' then td = true else td = false end
	for j=1,n do
		row = a[j]
		if td then 
			an = row[j][1]+row[j][2]
			row[j][1],row[j][2] = row[j][1]/an,row[j][2]/an
		else an = row[j]; row[j] = row[j]/an end
		for i,v in pairs(row) do
			if i~=j then row[i] = row[i]/an end
		end
		b[j] = b[j]/an
	end
end
Spmat = Spmat
spsolve = function(a,b,u,lan) -- Single step of simple over relaxation
	local n,epsm,an,bb,row,eps,jm,im = #b,0.0
	local um = 0.0
	local ua,un = {},{}--; for m=1,n do ua[m] = u[m] end 
	lan = lan or 1
	for j=1,n do 
		row,eps = a[j],b[j]
		for i,v in pairs(row) do
			if i==j then an = v 
			else eps = eps - v*u[i] end
		end
		eps = eps/an; un[j] = (1-lan)*u[j] + lan*eps
		if abs(un[j]-u[j])>abs(epsm) then epsm,jm,im = eps,j,i end
		um = max(um,abs(u[j]))
	end
	return un, epsm, um, jm, im
end
setfenv(spsolve,{table=table,type=type,abs=math.abs,max=math.max,pairs=pairs,print=print})

	
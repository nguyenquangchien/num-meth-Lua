-- File pde2fe.lua -- Functions for 2D finite element analysis

require"gauss"; require"spgauss"; require'pde2bv'
local NHF=13; local calcsh, deriv1tri, deriv2tri, ltoxy
local usave,utsave
setup2feeqs = function(eqs,nts,u) -- Set up 2D equations for finite elements
	-- u is of form u[k], nts = {nodes[], elements[], sides[]} nodes, elements and sides
	local u,bxext,byext = u or {}, 0.0, 0.0
	local uxx,ux,uy,uyy,uxa,uya,ut
	local nds,tri,sds = nts[1], nts[2], nts[3]
	local nnds,nel,nsds = #nds, #tri, #sds
	local fctuxx,fctux,fctu,fctuy,fctuyy = 1.0,FACT,FACT,FACT,1.0
	local a, b = Spmat.new(nnds,-nnds), Spmat.new(nnds,-nnds)
	local fv,fx,fu,fy,xi,yi,ui,xn,yn,mk = {},{},{},{},{},{},{},{},{},{}
	local fxx,fyy,mxx,vx,vy,he,trik,n1,n2,n3,m1
	local lipts = {{1/3,1/3,1/3},{.6,.2,.2},{.2,.6,.2},{.2,.2,.6}} -- Four point integration
	local lkpts = {{1/3,.6,.2,.2},{1/3,.2,.6,.2},{1/3,.2,.2,.6}} -- k=l1,l2,l3 at four points
	local fipts = {-27/48,25/48,25/48,25/48} -- Factors for integratopn points
	if #u~=nnds then -- Make sure u is a proper table
		local uxl = {}; for k=1,nnds do uxl[k] = 0.0 end; u = uxl 
	end
	mxx = calcsh(tri,nds); vx = sqrt(1/mxx) 
	mxx,vy = 0.0, vx -- Find maximum solution value
	for k=1,nnds do mxx = max(mxx, abs(u[k])) end
	if mxx~=0.0 then -- Scale probe factors to maximum solution value
		fctu = FACT*mxx
		fctux,fctuy = fctu*vx, fctu*vy
	end
	uxa,uya = deriv1tri(tri,nds,u) -- Get first derivative arrays
	eq,eqb = eqs[1], eqs[2] -- Differential, boundary equations
	for k=1,nnds do a[k],b[k] = {}, 0 end
	for k=1,nel do --Loop over triangles
		trik = tri[k] -- one by one
		n1,n2,n3,m1 = trik[1], trik[2], trik[3], trik[12]
		ui[1],ui[2],ui[3] = u[n1], u[n2], u[n3] -- Three node solutions
		xi[1],xi[2],xi[3] = nds[n1][1], nds[n2][1], nds[n3][1] -- arrays of x,y values
		yi[1],yi[2],yi[3] = nds[n1][2], nds[n2][2], nds[n3][2]
		ns = {n1,n2,n3}
		ux,uy,uxx,uyy = uxa[k],uya[k],0.0,0.0 -- Local derivatives
		he = trik[NHF] -- shape factors for triangle
		area,fxx,fyy = he[0], 0.0, 0.0 -- area factor
		bxext,byext = 0.0,0.0
		if truefd then -- extra first derivative terms from da/dx and da/dy -- d(a du/dx)/dx form 
			for j=1,3 do -- loop over three triangle nodes
				xt,yt,ut = xi[j],yi[j],ui[j]
				lpts = {0,0,0}; lpts[j] = 1
				fvv = eq(xt,yt,uxx,ux,ut,uy,uyy,k,m1,ns,lpts)
				bxext=bxext-he[j][2]*(eq(xt,yt,uxx+fctuxx,ux,ut,uy,uyy,k,m1,ns,lpts)-fvv)/fctuxx
				byext=byext-he[j][3]*(eq(xt,yt,uxx,ux,ut,uy,uyy+fctuyy,k,m1,ns,lpts)-fvv)/fctuyy
			end -- Now have da/dx and da/dy over triangle -- approximated as constant
		end 
		for j=1,4 do -- loop over 4 integration points
			lpts,fwt = lipts[j], fipts[j]; l1,l2,l3 = lpts[1], lpts[2], lpts[3]
			ut = ui[1]*l1 + ui[2]*l2 + ui[3]*l3 -- Values at integration points
			xt,yt = xi[1]*l1 + xi[2]*l2 + xi[3]*l3, yi[1]*l1 + yi[2]*l2 + yi[3]*l3
			fvv = eq(xt,yt,uxx,ux,ut,uy,uyy,k,m1,ns,lpts) -- Get partial derivatives
			fxx = fxx - fwt*(eq(xt,yt,uxx+fctuxx,ux,ut,uy,uyy,k,m1,ns,lpts)-fvv)/fctuxx -- Sum ax
			fx[j] = fwt*((eq(xt,yt,uxx,ux+fctux,ut,uy,uyy,k,m1,ns,lpts)-fvv)/fctux+bxext) -- bx terms
			fu[j] = fwt*(eq(xt,yt,uxx,ux,ut+fctu,uy,uyy,k,m1,ns,lpts)-fvv)/fctu -- c terms
			fy[j] = fwt*((eq(xt,yt,uxx,ux,ut,uy+fctuy,uyy,k,m1,ns,lpts)-fvv)/fctuy+byext) -- by terns
			fyy = fyy - fwt*(eq(xt,yt,uxx,ux,ut,uy,uyy+fctuyy,k,m1,ns,lpts)-fvv)/fctuyy -- Sum ay terms
			fv[j] = fwt*(fvv + bxext*ux + byext*uy) -- Fo terms
		end
		fxx,fyy = fxx*area, fyy*area -- common area weighting
		for i=1,3 do -- loop over triangle nodes
			nnd,lk = trik[i], lkpts[i] -- primary node number
			fb,fbx,fby,fc = 0.0, 0.0, 0.0, 0.0 -- b matrix factor, ux and uy factors
			for j=1,4 do -- Step over itegration points
				fb = fb + lk[j]*fv[j] -- b matrix weighting
				fbx,fby = fbx + lk[j]*fx[j], fby + lk[j]*fy[j] -- ux,uy weighting
			end
			fbx,fby = fbx*area, fby*area -- common area weithting
			arow = a[nnd] -- Row of a matrix for inserting
			hx,hy = he[i][2], he[i][3] -- h factor derivatives
			b[nnd] = b[nnd] - fb*area - ux*hx*fxx - uy*hy*fyy 
			for j=1,3 do -- step over 3 shape functions
				nc, hji,fc ,lj= trik[j], he[j], 0.0, lkpts[j]
				lx,ly = hji[2], hji[3] -- derivatives of shape functions
				fa = hx*lx*fxx + hy*ly*fyy -- uxx and uyy factors
				if fa~=0.0 then arow[nc] = (arow[nc] or 0) + fa end 
				if fbx~=0.0 then arow[nc] = (arow[nc] or 0) + fbx*lx end -- ux factors
				if fby~=0.0 then arow[nc] = (arow[nc] or 0) + fby*ly end -- uy factors
				for k=1,4 do fc = fc + fu[k]*lk[k]*lj[k] end -- sum u factors
				if fc~=0.0 then arow[nc] = (arow[nc] or 0) + fc*area end -- u factors
			end
		end
	end
	for k,sda in pairs(sds) do -- Loop over boundary nodes
		arow,ut,ux = a[k], u[k], 0.0
		n1,n2,n3 = k,sda[1][2],sda[2][2]
		m1,mk[1],mk[2] = nds[k][3],nds[n2][3],nds[n3][3] -- boundary condition markers
		for j=1,2 do -- loop over two sides
			s1 = sda[j]; n1 = s1[2]
			if sda[1][3]<0 then nsgn = -1 else nsgn = 1 end
			xt,yt = nds[k][1] - nds[n1][1], nds[k][2] - nds[n1][2]
			fu[j] = sqrt(xt^2 + yt^2)
			xn[j],yn[j] = -yt*nsgn, xt*nsgn -- normal vector for sides * length
		end
		lt = fu[1] + fu[2]
		for j=1,2 do 
			s1 = sda[j]; if s1[3]<0 then ntri = s1[4] else ntri = s1[3] end
			ux = ux + (uxa[ntri]*xn[j]+uya[ntri]*yn[j])/lt -- weighted normal field
		end
		fvv = eqb(k,ut,ux,m1,k)
		fb = (eqb(k,ut,ux+fctux,m1,k)-fvv)/fctux -- normal derivtaive factor
		if fb==0.0 then -- No derivative term present in boundary condition
			arow = {}; a[k] = arow
			arow[k] = (eqb(k,ut+fctu,ux,m1,k)-fvv)/fctu; b[k] = -fvv
		else -- Mixed or normal boundary condition
			if m1==mk[1] then mk[2] = m1 -- Test for equal boundary markers
			elseif m1==mk[2] then mk[1] = m1 end -- else isolated boundary node
			for j=1,2 do -- Now get mixed boundary conditions
				s1 = sda[j]; n1 = s1[2]
				if s1[3]<0 then ntri = s1[4] else ntri = s1[3] end
				ux = (uxa[ntri]*xn[j]+uya[ntri]*yn[j])/fu[j] -- mormal field
				fvv =  eqb(k,ut,ux,mk[j],k)
				fyy = (eqb(k,ut+fctu,ux,mk[j],k)-fvv)/fctu
				fxx = (eqb(k,ut,ux+fctux,mk[j],k)-fvv)/fctux
				if fyy~=0.0 then
					arow[k],arow[n1]=arow[k]-fyy*fu[j]/(fxx*3),arow[n1]-fyy*fu[j]/(fxx*6)
				end
				b[k] = b[k] + fvv*fu[j]/(fxx*2) - ux*fu[j]/2 + (fyy*u[k]*fu[j]/(fxx*3) + fyy*u[n1]*fu[j]/(fxx*6))
			end
		end
	end
	return a,b
end
setfenv(setup2feeqs,{FACT=1.e-4,NHF=NHF,table=table,Spmat=Spmat,sqrt=math.sqrt,
max=math.max,abs=math.abs,calcsh=calcsh,deriv1tri=deriv1tri,deriv2tri=deriv2tri,he=he,
pairs=pairs,truefd=false,print=print,table=table,whatis=whatis})

pde2fe = function(feqs,nts,u,tpsola) -- 2 spatial dimensions PDE Solver
	u = u or {}
	if type(tpsola)~='table' then tpsola = {tpsola} end
	local tpsol,rx,ry = tpsola[1] or SPM, tpsola[2], tpsola[3] -- Default to sparse matrix solution
	local umx,errm,krrm,a,b,uxy = 0.0
	local uold,ua,n,itmx = {}, {}, #nts[1], 0
	if #u~=n then for k=1,n do u[k] = u[k] or 0.0 end end
	if linear==1 then nnmx=1 else nnmx=NMX end -- One Newton cycle for linear eqns
	for k=1,n do uold[k],umx,ua[k] = u[k], max(umx,u[k]), 0.0 end
	for int=1,nnmx do
		a,b = setup2feeqs(feqs,nts,uold)
		if tpsol==SPM then -- Solve with spgauss, sparse matrix solver
			spgauss(a,b); ua = b -- b is new solution
		elseif tpsol==COE then -- Solve with Chebychev odd/even SOR
			if rx==nil then rx = cos(pi/sqrt(n)) end
			ua = pde2bvcoe(a,b,ua,rx,umx)
		elseif tpsol==SOR then -- Solve with SOR
			if rx==nil then rx = cos(pi/sqrt(n)) end
			ua = pde2bvsor(a,b,ua,rx,umx)
		else
			print('Unknown solution method specified in pde2fe.'); break
		end
		errm,krrm,umx,itmx = 0.0,0,0.0,itmx+1
		for k=1,n do errm,uold[k] = max(errm,abs(ua[k])), uold[k]+ua[k] end
		for k=1,n do umx,ua[k] = max(umx,abs(uold[k])), 0.0 end
		if nprint~=0 then print('Completed Newton iteration',int,'with correction',errm)
			io.flush() end
		if abs(errm)<ERR*umx then itmx = int; break end
	end
	if itmx==NMX then print('Maximum number of iterations exceeded in pde2fe!!')
		io.flush() end
	return uold,errm,itmx
end
setfenv(pde2fe,{spgauss=spgauss,max=math.max,cos=math.cos,pi=math.pi,sqrt=math.sqrt,
abs=math.abs,table=table,setup2feeqs=setup2feeqs,pde2bvcoe=pde2bvcoe,
pde2bvsor=pde2bvsor,linear=0,SPM=1,COE=2,SOR=3,NMX=20,ERR=1.e-4,print=print,
type=type,nprint=0,io=io,Spmat=Spmat,})

pde1stp2fe1t = function(feqs,tvals,nts,ua,tpsola,ut,utt) -- 2 spatial dimensions plus time PDE Solver
	local j, neq, t, h, h2,h2sq,hs,hx,hy,hz -- Local variables for function
	local un,unn,unew = {},{},{}
	local jfirst,nit,nitt,errm = 0,0,0 -- Number of iterations
	local neq = #nts[1]
	local tmin,tmax,ntval = tvals[1],tvals[2],tvals[3]
	local u = ua or {} -- Single arrays 	
	ut = ut or {}; utt = utt or {}
	-- Functions to add next time values and time derivatives 
	local fpde,fbv = feqs[1], feqs[2]
	feq = { -- Local functions to add time and time derivatives
		function(x,y,uxx,ux,up,uy,uyy,nt,m,nds,lds) -- General spatial point
		local ut,utt,n1,n2,n3,l1,l2,l3 -- Local coordinates and nodes
			n1,n2,n3,l1,l2,l3 = nds[1],nds[2],nds[3],lds[1],lds[2],lds[3]
			ut = (up - (un[n1]*l1+un[n2]*l2+un[n3]*l3))/h2 -- Interpolated un and unn
			utt = (up - (unn[n1]*l1+unn[n2]*l2+unn[n3]*l3))/h2sq
			return fpde(x,y,t,uxx,ux,up,uy,uyy,ut,utt,nt,m,nds,lds) -- With added time
		end,
		function(nd,u,und,nt,m) -- Boundary equations
			local ut,utt = (u - un[nd])/h2, (u - unn[nd])/h2sq
			return fbv(nd,t,u,und,ut,utt,nt,m) -- With added time
		end
	}
	t = tmin -- Initial t value
	hs = (tmax - t)/ntval -- Equal increments in t used, no adjusting step size
	-- If initial derivative not specified, use Backwards differencing for first 4 points
	if #ut~=neq then for m=1,neq do ut[m] = ut[m] or 0 end end -- Complete initial derivatives 
	if #utt~=neq then for m=1,neq do utt[m] = 0 end
		jfirst,h = 0,0.25*hs; h2,h2sq,hy,hx,hz = h,h*h,h,0,0 -- Set to BD parameters
	else -- Use TP or BD parameters
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
			u,errm,nitt = pde2fe(feq,nts,u,tpsola) -- Solve PDE at time t
			if nitt>nit then nit = nitt end -- Monitor maximun number of iterations
			-- New derivative values, same function as in feq() 
			for m=1,neq do ut[m],utt[m] = (u[m] - un[m])/h2,(u[m] - unn[m])/h2sq end
		until jfirst>=4 -- End of first interval repeat using Backwards difference
		if k==1 then 
		if bd~=false then jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = h,h*h,h,0,0 -- BD parameters
			else jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = hs/2,h*h/4,h,h/2,h*h/4 end -- TP parameters
		end
		if nprint~=0 then print('Completed time =',t,' with correction',errm); io.flush() end
	end -- End of main loop on t, now return solution array
	return u,errm,nit
end -- End of pde1stp2fe1t
setfenv(pde1stp2fe1t,{table=table,pde2fe=pde2fe,print=print,io=io,bd=false,nprint=0})

pde2fe1t = function(feqs,tvals,nts,uar,tpsola) -- Multi Time Step Solver
	local tpts,tsave = nil
	local ua,ub,uc,tl,ns = {},{},{}
	local ntp,ni,u,ut,utt
	local NMAX,ND = getfenv(pde2bv).NMX,10,0
	local neq,j = #nts[1]
	local nit,tmin,ntkeep,ntval,dtkeep = 0
	if tvals[5]==nil then tpts = false 
	else
		if type(tvals[5])~='table' then 
			tpts = false; print('Must specify array of x,y points for time data')
		else tpts = true end
	end
	tsave = tvals[5]
	uar = uar or {} -- No table may be passed
	if #uar<neq then u,ut,utt = uar[1], uar[2] or {}, uar[3] or {}
	else u,ut,utt = uar, {}, {} end
	u = u or {}
	if #u~=neq then for k=1,neq do u[k] = 0.0 end end
	for k=1,neq do ub[k] = u[k] end
	ua = usave(ua,ub) -- Save initial values
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
				k = k+1; t2[k] = tmin + j*(tmax - tmin)/(nsav)
				t3[k] = tvals[4][i] or NTVAL 
			end
			tmin = tmax
		end
		tvals[2],tvals[3] = t2, t3
	end
	ntkeep = max(#tvals[3]-1,1); tmin = tvals[2][1]
	if tpts then uta = utsave(uta,u,nts,tmin,tsave) end -- save initial values
	for i=1,ntkeep do -- Step over time increments
		ub,errm,ni = pde1stp2fe1t(feqs,{tmin,tvals[2][i+1],tvals[3][i]},nts,ub,tpsola,ut,utt) 
		tmin = tvals[2][i+1]
		if ni==NMAX then 
			print("Error: Maximum number of iterations exceeded in pde2bv")
			print("Results may not be accurate! Time parameters =",tl,tvals[2][i],ns);io.flush()
		end
		nit = max(nit,ni); ua = usave(ua,ub)
		if tpts then uta = utsave(uta,ub,nts,tmin,tsave) end -- save if requested
		if nprint~=0 then print('Completed time =',tmin,' with correction',errm); io.flush() end
		if seeplot~=nil then print('Plot at time =',tmin); io.flush()
		if type(seeplot)=='table' then -- Possible popup plot
			sxy = ltoxy(ub,seeplot[1],seeplot[2],nts[1],nts[2]);	splot(sxy) 
		end end
	end
	if tpts then return ua,uta,errm,nit else return ua,errm,nit end
end
setfenv(pde2fe1t,{table=table,type=type,print=print,io=io,splot=splot,seeplot=nil,max=math.max,
NKEEP=10,NTVAL=10,pde1stp2fe1t=pde1stp2fe1t,getfenv=getfenv,nprint=0,cplot=cplot,ltoxy=ltoxy})

pde2fe1tqs = function(feqs,tvals,nts,uar,tpsola) -- Quick scan log time solver
	local nps,npts,tpts -- save per decade, additional steps per save
	local ttvals,nl,nu,fact = {}
	local nt,i,j,neq = #tvals,2,0
	local ua,ub,uta,nit,ni = {},{},{},0 
	local u,ut,utt,ntv,tmp,sxy
	local neq = #nts[1]
	local nprint = nprint or getfenv(pde2bv).nprint
	if tvals[5]==nil then tpts = false 
	else
		if type(tvals[5])~='table' then 
			tpts = false; print('Must specify array of x,y points for time data')
		else tpts = true end
	end
	local tsave = tvals[5]
	if #uar<neq then u,ut,utt = uar[1], uar[2] or {}, uar[3] or {}
	else u,ut,utt = uar, {}, {} end
	if nt<2 then print('Error, must specify two times in pd2vfe1tqs')
		return end
	nps,npts = floor(tvals[3] or NPS),floor(tvals[4] or NPTS)
	if #u~=neq then for k=1,neq do u[k] = 0.0 end end
	for k=1,neq do ub[k] = u[k] end; ua = usave(ua,ub) -- Save initial values
	ntv,nl = nps*npts,tvals[1]; tmp = (tvals[2][1]-tvals[1])/ntv
	nu = nl + tmp
	if tpts then uta = utsave(uta,ub,nts,nl,tsave) end -- save initial values
	for i=1,ntv do -- loop over initial interval
		ub,errm,ni = pde1stp2fe1t(feqs,{nl,nu,1,1},nts,ub,tpsola,ut,utt)
		nit = max(nit,ni) 
		if tpts then uta = utsave(uta,ub,nts,nu,tsave) end -- save time values
		if nprint~=0 then print('In quick scan, initial time =',nu);io.flush() end
		nl,nu = nu,nu+tmp
	end
	ua = usave(ua,ub,nx,ny) -- Save values after initial interval
	if seeplot~=nil then if type(seeplot)=='table' then
		print('Plot at time =',nl); io.flush() -- Possible popup plot
		sxy = ltoxy(ub,seeplot[1],seeplot[2],nts[1],nts[2]); splot(sxy)  
	end end
	nps,npts = floor(tvals[3] or NPS),floor(tvals[4] or NPTS)
	fact = 10^(1/(nps*npts)) -- Factor between steps
	nl = 10^(floor(log10(tvals[2][1])))
	nu = 10^(ceil(log10(tvals[2][2])))*1.00000001/fact
	while nl<=nu do -- Step over log time spacings
		ub,errm,ni = pde1stp2fe1t(feqs,{nl,nl*fact,1},nts,ub,tpsola,ut,utt)
		nit,j,nl = max(nit,ni), j+1, nl*fact
		if tpts then uta = utsave(uta,ub,nts,nl,tsave) end -- save if requested
		if nprint~=0 then print('In quick scan, t =',nl); io.flush() end
		if j==npts then -- Time to save array
			i,j = i+1,0; ua = usave(ua,ub) -- Track saved number
			if seeplot~=nil then if type(seeplot)=='table' then
				print('Plot at time =',nl); io.flush() -- Possible popup olot
				sxy = ltoxy(ub,seeplot[1],seeplot[2],nts[1],nts[2]); splot(sxy)  
			end end
			if nprint~=0 then 
				print('Saved at Time = ',nl,' Number of iterations in quick scan = ',nit,'\n')
				io.flush(); nit = 0
			end
		end
	end
	if tpts then return ua,uta,errm,nit else return ua,errm,nit end
end
setfenv(pde2fe1tqs,{print=print,type=type,floor=math.floor,max=math.max,log10=math.log10,
ceil=math.ceil,nprint=0,io=io,NPS=1,NPTS=10,getfenv=getfenv,seeplot=nil,usave=usave,
pde1stp2fe1t=pde1stp2fe1t,ltoxy=ltoxy,splot=splot,whatis=whatis})

function areatr(t,n) -- Area of triangle t at nodes n = {...{n1,n2,n3}...}
	local nd1, nd2, nd3 = n[t[1]], n[t[2]], n[t[3]]
	local x1,x2,x3,y1,y2,y3 = nd1[1], nd2[1], nd3[1], nd1[2], nd2[2], nd3[2]
	return math.abs(0.5*((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)))
end

function hfunc(tri,nd,tn) -- Evaluation of shape functions
	local t1,h = tri[tn], {}
	local x1,x2,x3 = nd[t1[1]][1], nd[t1[2]][1], nd[t1[3]][1]
	local y1,y2,y3 = nd[t1[1]][2], nd[t1[2]][2], nd[t1[3]][2]
	local a2 = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) --2*area
	h[0] = a2/2 --math.abs(a2)
	h[1] = {(x2*y3-x3*y2)/a2,(y2-y3)/a2,(x3-x2)/a2}
	h[2] = {(x3*y1-x1*y3)/a2,(y3-y1)/a2,(x1-x3)/a2}
	h[3] = {(x1*y2-x2*y1)/a2,(y1-y2)/a2,(x2-x1)/a2}
	return h
end		
	
function hsfunc(tri,nd,tn) -- Evaluation of shape functions
	t1 = tri[tn]
	x,y = {},{}; ha = {}
	for i=1,3 do x[i],y[i] = nd[t1[i]][1], nd[t1[i]][2] end
	s2 = 'U(x,y) = '; s3 = 'U(x,y) = '
	for nd = 1,3 do
		a = {{1,x[1],y[1]},{1,x[2],y[2]},{1,x[3],y[3]}}
		if nd==1 then b = {1,0,0}
		elseif nd==2 then b = {0,1,0}
		else b = {0,0,1} end
		gauss(a,b)
		ha[nd] = {b[1],b[2],b[3]}
		if nprint~=0 then 
			s1 = tostring(b[1])..' + ('..tostring(b[2])..')*x + ('..tostring(b[3])..')*y'
			s4 = 'h['..tostring(tn)..']['..tostring(t1[nd])..']'
			s2 = s2..'U['..tostring(t1[nd])..']*{'..s1..'}'
			s3 = s3..'U['..tostring(t1[nd])..']*'..s4
			if nd<3 then s2,s3 = s2..' + ', s3..' + ' end
		end
	end
	return ha, s2,s3
end
setfenv(hsfunc,{tostring=tostring,gauss=gauss,nprint=0})
	
function tript(nds,tri,pt) -- Find triangle for pt = {x,y}
	local ti,n1,n2,n3,x1,x2,x3,xp,y1,y2,y3,yp,den,l1,l2
	local nnds = #tri
	nlst = min(nnds-2,max(2,nlst)) -- start at last found
	local nl1,nl2 = nlst,nlst+1
	while nl1>0 or nl2<=nnds do -- step over triangles
		if nl1>0 then -- Down from last found
			ti = tri[nl1]; n1,n2,n3 = nds[ti[1]], nds[ti[2]], nds[ti[3]] 
			x1,x2,x3,xp = n1[1], n2[1], n3[1], pt[1]
			y1,y2,y3,yp = n1[2], n2[2], n3[2], pt[2]
			den = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)
			l1 = ((xp-x3)*(y2-y3)- (yp-y3)*(x2-x3))/den
			if l1>=LMIN and l1<=LMAX then
				l2 = ((yp-y3)*(x1-x3) - (xp-x3)*(y1-y3))/den
				if l2>=LMIN and l2<=LMAX then 
					if l1+l2<=LMAX then nlst = nl1; return nlst,l1,l2 end -- Found it!
				end
			end
			nl1 = nl1-1
		end
		if nl2<=nnds then -- up from last found
			ti = tri[nl2]; n1,n2,n3 = nds[ti[1]], nds[ti[2]], nds[ti[3]] 
			x1,x2,x3,xp = n1[1], n2[1], n3[1], pt[1]
			y1,y2,y3,yp = n1[2], n2[2], n3[2], pt[2]
			den = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)
			l1 = ((xp-x3)*(y2-y3)- (yp-y3)*(x2-x3))/den
			if l1>=LMIN and l1<=LMAX then
				l2 = ((yp-y3)*(x1-x3) - (xp-x3)*(y1-y3))/den
				if l2>=LMIN and l2<=LMAX then 
					if l1+l2<=LMAX then nlst = nl2; return nlst,l1,l2 end -- Found it!
				end
			end
			nl2 = nl2+1
		end
	end
	return nil, l1, l2 -- Not found!!
end	
setfenv(tript,{table=table,min=math.min,max=math.max,nlst=2,LMIN=-1.e-7,LMAX=1+1.e-7})

intptri = function(u,nds,tri,pt) -- Interpolate solution over triangular grid	
	local trpt,l1,l2,tr = tript(nds,tri,pt) -- Find triangle for point
	if trpt~=nil then
		tr = tri[trpt]
		return u[tr[1]]*l1 + u[tr[2]]*l2 + u[tr[3]]*(1-l1-l2), trpt
	end
	return trpt, trpt -- Point not in structure
end 

function normvec(nd,sds,nds) -- Find normal surface vectors around node nd
	local found,nv,n1,n2,ne = 0, {}
	for i=1,#sds do
		sd = sds[i] -- test sides in turn 
		if sd[1]==nd or sd[2]==nd then
			if sd[3]==-1 or sd[4]==-1 then
				found = found+1
				if sd[1]== nd then n1,n2 = nd,sd[2] 
				else n1,n2 = sd[1],nd end
				if sd[3]==-1 then dir,ne = -1,sd[4] else dir,ne = 1,sd[3] end
				x1,y1,x2,y2 = nds[n1][1],nds[n1][2],nds[n2][1],nds[n2][2]
				dy,dx = y2-y1, x2-x1; nv[found] = {dir*dy, -dir*dx, math.sqrt(dx^2 + dy^2),ne,nd}
			end
		end
		if found==2 then break end
	end
	if found~=2 then print('Error, two boundary surfaces not found for node:',nd) end
	return nv
end	
function normvecs(sds,nds) -- Find all normal surface vectors 
	local nsnd,nvs = 0,{}
	for i=1,#nds do -- Test all nodes to see if on surface
		if nds[i][3]~=0 then -- On surface
			nsnd = nsnd+1; nvs[nsnd] = normvec(i,sds,nds)
		end
	end
	return nvs
end	

intgtri = function(fxy,tri,nds) -- Quadratic integration over triangle
	local nd1, nd2, nd3 = nds[tri[1]], nds[tri[2]], nds[tri[3]]
	local x1,x2,x3,y1,y2,y3 = nd1[1], nd2[1], nd3[1], nd1[2], nd2[2], nd3[2]
	local ara = areatr(tri,nds)
	return ara*(fxy(0.5*(x1+x2),0.5*(y1+y2)) + fxy(0.5*(x1+x3),0.5*(y1+y3)) +
		fxy(0.5*(x2+x3),0.5*(y2+y3)))/3, ara
end
intgtric = function(fxy,tri,nds) -- Cubic integration over triangle
	local nd1, nd2, nd3 = nds[tri[1]], nds[tri[2]], nds[tri[3]]
	local x1,x2,x3,y1,y2,y3 = nd1[1], nd2[1], nd3[1], nd1[2], nd2[2], nd3[2]
	local ara = areatr(tri,nds)
	return ara*((fxy(.6*x1+.2*(x2+x3),.6*y1+.2*(y2+y3)) +
		fxy(.6*x2+.2*(x1+x3),.6*y2+.2*(y1+y3)) +
		fxy(.6*x3+.2*(x1+x2),.6*y3+.2*(y1+y2)))*25/48 -
		fxy((x1+x2+x3)/3,(y1+y2+y3)/3)*27/48), ara
end
deriv1tri = function(tri,nds,u) -- First derivatives ux,uy over triangles 
	local ux,uy = {},{}
	local nt = #tri -- Number of triangles
	local trii,n1,n2,n3,hfx
	for i=1,nt do -- Step over triangles
		trii = tri[i]
		n1,n2,n3 = trii[1], trii[2], trii[3] -- Nodes of triangle
		hfx = trii[NHF] or hfunc(tri,nds,i)
		ux[i] = u[n1]*hfx[1][2] + u[n2]*hfx[2][2] + u[n3]*hfx[3][2]
		uy[i] = u[n1]*hfx[1][3] + u[n2]*hfx[2][3] + u[n3]*hfx[3][3]
	end
	return ux,uy
end; derivtri = deriv1tri
deriv2tri = function(tri,nds,u) -- Second derivative uxx, uyy over triangles
	local a,b,uxx,uyy = {},{},{},{}; b[5] = 0.0
	local nt1,nt2
	local nbnds,bnds,tadj = 0, {},{}
	nt = #tri -- Number of triangles
	for i=1,nt do -- Step over triangles
		trii = tri[i]
		n1,n2,n3 = trii[1], trii[2], trii[3] -- Nodes of triangle
		ndk,nex = 3,{n1,n2,n3} -- Number of adjacent nodes
		for j=4,6 do -- Find nodes of adjacent triangles
			trix = trii[j] -- Adjacent triangle
			if trix~=-1 and trii[12]==tri[trix][12] then 
				for k=1,3 do -- Nodes of adjacent triangle
					nxj = tri[trix][k]
					if nxj~=n1 and nxj~=n2 and nxj~=n3 then
						ndk = ndk+1; nex[ndk],tadj[ndk-3] = nxj, trix
					end
				end
			end
		end
		if ndk>5 then -- have 6 nodes		
			for j=1,ndk do -- Now have 4 to 6 nodes
				nd = nds[nex[j]]; xn,yn = nd[1],nd[2] -- node number
				a[j] = {1, xn, yn, xn^2, yn^2, xn*yn}
				b[j] = u[nex[j]]
			end
			gauss(a,b,ndk)
			uxx[i],uyy[i] = 2*b[4], 2*b[5]
		else -- only one adjacent node
			nbnds = nbnds+1; bnds[nbnds] = {i, ndk, tadj[1], tadj[2]}
		end
	end
	for i,ntp in pairs(bnds) do -- triangles with one or two adjacent triangle
		trii,ndk,nt1,nt2 = ntp[1], ntp[2], ntp[3], ntp[4]
		if ndk==5 then
			if uxx[nt1]==nil then
				if uxx[nt2]~=nil then 
				uxx[trii],uyy[trii],bnds[i] = uxx[nt2],uyy[nt2],nil 
				else print('both are nil') end
			elseif uxx[nt2]==nil then 
			uxx[trii],uyy[trii],bnds[i] = uxx[nt1],uyy[nt1],nil 
			else 
			uxx[trii],uyy[trii],bnds[i] = 0.5*(uxx[nt1]+uxx[nt2]), 0.5*(uyy[nt1]+uyy[nt2]),nil end
		else
			if uxx[nt1]~=nil then uxx[trii],uyy[trii],bnds[i] = uxx[nt1],uyy[nt1],nil end
		end
	end 
	for i,ntp in pairs(bnds) do print('still have nils at triangle',bnds[1]) end
	for i,ntp in pairs(bnds) do -- do again if needed
		trii,ndk,nt1,nt2 = ntp[1], ntp[2], ntp[3], ntp[4]
		if ndk==2 then
			if uxx[nt1]==nil then
				if uxx[nt2]~=nil then uxx[trii],uyy[trii],bnds[i] = uxx[nt2],uyy[nt2],nil end
			elseif uxx[nt2]==nil then uxx[trii],uyy[trii],bnds[i] = uxx[nt2],uyy[nt2],nil 
			else uxx[trii],uyy[trii],bnds[i] = 0.5*(uxx[nt1]+uxx[nt2]), 0.5*(uyy[nt1]+uyy[nt2]),nil end
		else
			if uxx[nt1]~=nil then uxx[trii],uyy[trii],bnds[i] = uxx[nt1],uyy[nt1],nil end
		end
	end
	for i,ntp in pairs(bnds) do print('still have nils at triangle',bnds[1]) end
	return uxx,uyy
end	

function calcsh(tri,nds) -- Calculate shape functions for all triangles
	local ara,hi = 0.0
	for i=1,#tri do 
		hi = hfunc(tri,nds,i); tri[i][NHF] = hi; ara = ara+hi[0]*.5 end
	return ara -- return total area
end

function toxysol(u,xa,ya,nodes,triangles)
	local k,sol,sxy= 1, {{},{},{}},{}
	local nx,ny = #xa, #ya
	local x,y,val,row
	if triangles== nil then nodes,triangles = nodes[1],nodes[2] end
	for j=1,ny do y = ya[j]
		for i=1,nx do
			x = xa[i]; val = intptri(u,nodes,triangles,{x,y})
			if val~=nil then sol[1][k],sol[2][k],sol[3][k] = x,y,val; k = k+1 end 
		end
	end
	for i=1,nx do x = xa[i]
		row = {}
		for j=1,ny do
			y = ya[j]; val = intptri(u,nodes,triangles,{x,y})
			if val~=nil then sol[1][k],sol[2][k],sol[3][k] = x,y,val; k = k+1 end
			if val~=nil then row[j] = val end
		end
		sxy[i] = row
	end
	return sol,sxy
end
ltoxy = toxysol

function read_nts(file)
	local nts = {}
	local fx = io.open(file,'r')
	local nc,i,jj,k,kk = 25,-1,0
	for lines in fx:lines() do -- Read line by line
		k,i,jj = 0,i+1,0
		for j=1,nc do -- Data separators may be single tab, commas or space
			kk = string.find(lines,":",k+1) -- possible ':' separator
			if kk==nil then -- ':' is not a seperator?
				kk = string.find(lines," ",k+1) -- possible space separator
				if kk==nil then break end -- Valid separator not found
			end
			nmb = tonumber(string.sub(lines,k,kk-1)); k = kk+1
			if nmb~=nil and j>1 then nts[i] = nts[i] or {}; jj = jj+1 nts[i][jj] = nmb end 
		end 	
		if jj<nc then nmb = tonumber(string.sub(lines,k))
			if nmb~=nil then jj = jj+1; nts[i] = nts[i] or {}; nts[i][jj] = nmb end 
		end
	end
	fx:close(); return nts -- return data array
end
function readnts(file,scale)
	local row
	local nds = read_nts(file..'.n')
	local tri = read_nts(file..'.e')
	local sd = read_nts(file..'.s')
	local jsd,sds = 0,{}
	for i=1,#nds do
		if nds[i][3]>100 then nds[i][3] = 0 end
	end
	for i=1,#tri do
		row = tri[i]
		for j=1,9 do if row[j]~=-1 then row[j] = row[j]+1 end end
	end
	if scale~=nil then 
		for i=1,#nds do
			row = nds[i]; row[1],row[2] = scale*row[1],scale*row[2]
		end
		for i=1,#tri do
			row = tri[i]; row[10],row[11] = scale*row[10],scale*row[11]
		end
	end
	for i=1,#sd do -- Save only boundary sides
		row = sd[i]
		if row[3]==-1 or row[4]==-1 then
			for j=1,4 do if row[j]~=-1 then row[j] = row[j]+1 end end
			n1,n2 = row[1],row[2]
			if sds[n1]==nil then sds[n1] = {} end
			if sds[n1][1]==nil then sds[n1][1] = {row[1],row[2],row[3],row[4]}
			else sds[n1][2] = {row[1],row[2],row[3],row[4]} end
			if sds[n2]==nil then sds[n2] = {} end
			if sds[n2][1]==nil then sds[n2][1] = {row[1],row[2],row[3],row[4]}
			else sds[n2][2] = {row[1],row[2],row[3],row[4]} end
		end
	end
	for i,_ in pairs(sds) do -- Reverse directions if needed
		s1 = sds[i][1]; s2 = sds[i][2]
		if s1[1]~=i then s1[1],s1[2],s1[3],s1[4] = s1[2],s1[1],s1[4],s1[3] end
		if s2[1]~=i then s2[1],s2[2],s2[3],s2[4] = s2[2],s2[1],s2[4],s2[3] end
	end
	return nds,tri,sds
end	

usave = function(u1,u2)
	local nk,neq,up = #u1,#u2, {}
	for k=1,neq do up[k] = u2[k] end
	u1[nk+1] = up
	return u1
end
setfenv(usave,{table=table,print=print,whatis=whatis})
	
utsave = function(u1,u2,nts,t,px)
	u1 = u1 or {}
	local nk,nt = #u1,#px
	if nk~=nt+1 then for i=1,nt+1 do u1[i] = {} end end
	nk = #u1[1]+1
	u1[1][nk] = t
	for i=1,nt do u1[i+1][nk] = intptri(u2,nts[1],nts[2],px[i]) end
	return u1
end
setfenv(utsave,{table=table,intptri=intptri,whatis=whatis})

local glpt = {-math.sqrt(3/5), 0, math.sqrt(3/5)}
local wt = {5/9, 8/9, 5/9} -- parameters for Gauss integ, limits -1 to 1, both dimensions
intg2d = function(fxy,lx,ly) -- two dimensional Gaussian integration
	local lxl,lxu = lx[1], lx[2]
	local lyl,lyu = ly[1], ly[2]
	local a0,a1 = 0.5*(lxu + lxl), 0.5*(lxu - lxl)
	local b0,b1 = 0.5*(lyu + lyl), 0.5*(lyu - lyl)
	local ans = 0.0
	for i=1,3 do
		for j=1,3 do
			ans = ans + wt[i]*wt[j]*fxy(a0+a1*glpt[i], b0+b1*glpt[j])
		end
	end
	return ans*(lxu-lxl)*(lyu-lyl)/4
end


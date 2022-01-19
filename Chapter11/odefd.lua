-- /* File odefd.lua */

require"Matrix" -- Needed for matrix algebra
ode2bvfd = function(eqs,x,u)
	local neq,s,upp,up,uppi,upi,ui = #u, {}, {}, {}, {}, {}, {} -- neq = #equations
	a,b,c,d = Matrix.new(neq,neq),Matrix.new(neq,neq),Matrix.new(neq,neq),Matrix.new(neq,1)
	local alfi,e,f = 0, {}, {}
	local feq ,nend = eqs[1], NMAX
	local fval,fu,fctupp,fctup,fctu = {}, {}, {}, {}, {}
	local uppi,upi,ui,cmax,imax,umax = {}, {}, {}, {}, {}, {}
	for m=1,neq do fctupp[m],fctup[m],fctu[m] = FACT,FACT,FACT end
	
	nx = #x; nxm1 = nx-1 -- number of x values -- associated with index 'i'
	for i=1,nxm1 do s[i] = x[i+1] - x[i] end
		
	bound = function(nb,nxb) -- Function to evaluate specified boundary conditions
		xi = x[nxb]
		for m=1,neq do upi[m],ui[m] = up[m][nxb],u[m][nxb] end
		eqs[nb](fval,ui,upi)
		if nb==2 then c1,c2 = -1/s[1],1/s[1] else c1,c2 = 1/s[nx-1],-1/s[nx-1] end
		for m=1,neq do 
			d[m][1] = fval[m]
			upi[m] = upi[m]+fctup[m]; eqs[nb](fu,ui,upi); upi[m] = upi[m]-fctup[m]
			for n=1,neq do -- Probe up factor
				fjk = (fu[n] - fval[n])/fctup[m]
				b[n][m] = c1*fjk
				if nb==2 then c[n][m] = c2*fjk else a[n][m] = c2*fjk end
			end -- Then probe u factor below
			ui[m] = ui[m]+fctu[m]; eqs[nb](fu,ui,upi,ui); ui[m] = ui[m]-fctu[m] 
			for n=1,neq do b[n][m] = b[n][m] + (fu[n] - fval[n])/fctu[m] end
		end
	end -- Return boundary equations involving a,b,c,d
			
	for m=1,neq do upp[m],up[m] = {},{} end -- Derivative arrays
	for ni=1,NMAX do -- ni = iteration number count
		for m=1,neq do fctupp[m],fctup[m],fctu[m] = 0,0,0 end 
		for i=2,nx-1 do -- Calculate second derivative and first derivative for present itteration
			si,sisi,alfi = s[i],(s[i]+s[i-1]),s[i]/s[i-1]
			c1,c2,c3 = 2/(si*sisi),-2*(alfi+1)/(si*sisi),2*alfi/(si*sisi)
			for m=1,neq do
				fctu[m] = fctu[m] + abs(u[m][i]*si) 
				duu = c1*u[m][i+1] + c2*u[m][i] + c3*u[m][i-1]
				upp[m][i],fctupp[m] = duu,fctupp[m]+abs(duu*si) -- Second derivative array
			end
			c1,c2,c3 = 1/(alfi*sisi),(alfi-1)/si,-alfi/sisi
			for m=1,neq do
				duu = c1*u[m][i+1] + c2*u[m][i] + c3*u[m][i-1]
				up[m][i],fctup[m] = duu,fctup[m]+abs(duu*si) -- First derivative array
			end
		end
		alfi = s[2]/s[1]
		for m=1,neq do -- Special treatment for lower end point
			upp[m][1] = upp[m][2]
			up[m][1] = (-u[m][3] + u[m][2]*(1+alfi)^2 - u[m][1]*alfi*(2+alfi))/(alfi*(s[2]+s[1]))
		end
		alfi = s[nxm1]/s[nx-2]
		for m=1,neq do -- Special treatment for upper end point
			upp[m][nx] = upp[m][nxm1]
			up[m][nx] = (u[m][nx]*(1+2*alfi) - u[m][nxm1]*(1+alfi)^2 + 
				u[m][nx-2]*alfi^2)/(alfi*(s[nxm1]+s[nx-2]))
		end
		for m=1,neq do -- Add in end points to protect against large values with small changes
			fctupp[m] = fctupp[m] + abs(upp[m][1]) + abs(upp[m][nx])
			fctup[m] = fctup[m] + abs(up[m][1]) + abs(up[m][nx])
			fctu[m] = fctu[m] + abs(u[m][1]) + abs(u[m][nx])
		end
		for m=1,neq do -- Average values of variables and derivatives
			fctupp[m],fctup[m],fctu[m] = FACT*fctupp[m], FACT*fctup[m],FACT*fctu[m]
			if fctupp[m]==0 then fctupp[m] = FACT end -- Try to protect against zero values
			if fctup[m]==0 then fctup[m] = FACT end
			if fctu[m]==0 then fctu[m] = FACT end
		end
		if umin[1]~=nil then -- limit fctu values based upon umin values
			for m=1,neq do fctu[m] = max(fctu[m],FACT*umin[m]) end
		end
		
		bound(2,1) -- Evaluate lower boundary conditions -- Returns b,c and d coefficients as matrices
		gi = b^-1 -- Matrix algebra for e[] and f[]
		if type(gi)=='number' then 
			printf('Error in left boundary values\nCheck boundary equations\n') 
			return end
		e[2],f[2] = -gi*c,-gi*d
							
		for i=2,nx-1 do -- Set up a,b,c,d arrays and save e and h as arrays -- do not need to save a,b,c,d
			xi,si,sisi,alfi = x[i],s[i],(s[i]+s[i-1]),s[i]/s[i-1]
			for m=1,neq do -- Set up arrays for local values of second and first derivatives
				uppi[m],upi[m],ui[m] = upp[m][i],up[m][i],u[m][i]
				for n=1,neq do a[n][m],b[n][m],c[n][m] = 0,0,0 end -- Zero a,b,c arrays from previous iteration
			end
			feq(fval,xi,ui,upi,uppi,i) -- Evaluate equations -- returns difference from desired zero value
			for m=1,neq do -- increment each variable in order - - k = variable number
				d[m][1] = fval[m] -- Set d[] array value
				c1,c2,c3 = 2*alfi/(si*sisi), -2*(alfi+1)/(si*sisi), 2/(si*sisi)
				uppi[m] = uppi[m] + fctupp[m]
				feq(fu,xi,ui,upi,uppi,i) -- Probe upp factor for a,b,c components
				for n=1,neq do -- Now collect changes over equations -- j = equation number
					fjk = (fu[n] - fval[n])/fctupp[m] -- Update a,b,c values due to upp
					a[n][m],b[n][m],c[n][m] = a[n][m]+c1*fjk,b[n][m]+c2*fjk,c[n][m]+c3*fjk 
				end
				c1,c2,c3 = -alfi/sisi, (alfi-1)/si, 1/(alfi*sisi)
				uppi[m],upi[m] = uppi[m] - fctupp[m], upi[m] + fctup[m]
				feq(fu,xi,ui,upi,uppi,i) -- probe up factor for a,b,c components
				for n=1,neq do
					fjk = (fu[n] - fval[n])/fctup[m] -- Update a,b,c values due to up
					a[n][m],b[n][m],c[n][m] = a[n][m]+c1*fjk,b[n][m]+c2*fjk,c[n][m]+c3*fjk 
				end
				upi[m],ui[m] = upi[m] - fctup[m], ui[m] + fctu[m]
				feq(fu,xi,ui,upi,uppi,i) -- Probe u factor and update b values due to u
				for n=1,neq do b[n][m] = b[n][m] + (fu[n]-fval[n])/fctu[m] end  
				ui[m] = ui[m] - fctu[m]
			end
			-- Solve tridagonal matrix equations using matrix algebra
			gi = (a*e[i] + b)^-1-- Matrix algebra used here to calculate e[] and f[]
			e[i+1] = -gi*c; f[i+1] = -gi*(d + a*f[i])
		end -- Now have [i] array of e[] and f[] Matrix factors
		
		bound(3,nx) -- Evaluate upper boundary condition -- Calculates a,b and d coefficients
		gi = (a*e[nx] + b)^-1
		if type(gi)=='number' then 
			printf('Error in right boundary values\nCheck boundary equations\n')
			return end
		d = -gi*(a*f[nx] + d)
		
		for m=1,neq do -- Zero error factors
			cmax[m],imax[m],umax[m]  = 0,0,0 
		end
		for i=nx,1,-1 do -- Now go from point nx to point 1 calculating corrections to equations -- in d array
			for m=1,neq do 
				du = d[m][1] -- Each correction is taken in turn
				c1 = abs(du); if c1>umax[m] then umax[m] = c1 end
				c2 = abs(u[m][i]) + c1; if c2~=0.0 then c1 = c1/c2 end
				if c1>cmax[m] then cmax[m]=c1; imax[m]=i end
				u[m][i] = u[m][i] + du -- Update solutions -- This is what we work for 
			end
			if i==1 then break end
			d = e[i]*d + f[i] -- Matrix algebra for next correction at point i-1
		end -- Back for next point, i-1
		
		if nprint~=0 then -- Print iteration data
			printf("-- %i-- Iteration number, Maxumum relative, absolute corrections are: \n",ni)
			for m=1,neq do	printf("(%i) %e, %e, at %i ; ", m, cmax[m], umax[m],imax[m]) end
			printf("\n"); io.flush()
		end
		
		c1 = 1-- Now see if solution meets accuracy criteria
		for m=1,neq do	if cmax[m]>ERROR then c1 = 0 end end -- Relative accuracy criteria
		if c1==1 then nend=ni; break end -- Relative error condition met
		if umin[1]~=nil then -- Absolute accuracy limits specified
			c1 = 1
			for m=1,neq do 
				if umin[m]~=nil then -- Limit for variable m specified
					if umin[n]~=0 then if umax[m]>umin[m] then c1 = 0 end -- Test absolute accuracy
					else if umax[m]>fctu[m]*ERROR/FACT then c1 = 0 end -- Default absolute accuracy
				end end
			end
			if c1==1 then nend=ni; break end -- Absolute error condition met
		end
	end
	return nend, cmax, umax, upp, up -- Solution values returned in u calling array
	-- Derivatives returned in case user needs values, they don't have to be recalculated
end
setfenv(ode2bvfd,{type=type,abs=math.abs,max=math.max,Matrix=Matrix,
	ERROR=1.e-5,umin={},FACT=1.e-6,NMAX=50,nprint=0,printf=printf,io=io,
	math=math,table=table,unpack=unpack,ode2bvfd=ode2bvfd})

odefd = function(eqs,x,u) -- Simple FD interface
	local nx,neq,x1,x2 = #x, #u
	local u1,u2 = {},{}
	if type(u[1])~='table' then neq, nuv, u = 1, neq, {u} end
	local nuv = #u[1]
	for i=1,neq do u1[i],u2[i] = u[i][1],u[i][nuv] end 
	if nx<=3 then -- Generate linear step distribution
		x1,x2,nx = x[1],x[2],x[3] or 2001; nx = 2*math.floor(nx/2)+1
		local dx = (x2-x1)/(nx-1)
		for i=1,nx do x[i] = x1+(i-1)*dx end
	end -- If u array not specified, generate below
	if nuv~=nx then 
		for j=1,neq do for i=1,nx do u[j][i] = u1[j]+(i-1)/(nx-1)*(u2[j]-u1[j]) end end
	end	
	if type(eqs)=='function' then eqs = {eqs} end -- Boundary functions not given
	eqs[2]  = eqs[2] or function(fbv,u,up) 
		for i=1,neq do fbv[i] = u[i] - u1[i] end
	end
	eqs[3] = eqs[3] or function(fbv,u,up) 
		for i=1,neq do fbv[i] = u[i] - u2[i] end
	end
	return {x,unpack(u)},ode2bvfd(eqs,x,u)
end	
setfenv(odefd,getfenv(ode2bvfd))
odeerror = function(si1,si2) -- ODE error estimator for two solutions 
	require"intp"
	local err,ny,s1,s2 = {}, #si1
	local nd1,nd2 = #si1[1], #si2[1]
	local fac = 3
	if odebiv==odebrk then fac = 15 end
	for i=1,ny do err[i] = {} end
	if nd1>nd2 then s1,s2 = si1,si2 
	else nd1,nd2,s1,s2 = nd2,nd1,si2,si1 end
	if nd1~=(nd2-1)*2+1 then print('Error: Incorrect files input to odeerror, nd1,nd2 = ',
		nd1,nd2) end
	err[1] = s1[1] -- Set x values for arrays
	for k=2,ny do for i=1,nd1 do
			err[k][i] = (intp(s2[1],s2[k],s1[1][i]) - s1[k][i])/fac
	end end
	return err
end
odefde = function(eqs,x,u) -- Simple FD interface with error estimate return
	local s1,n1,jj = odefd(eqs,x,u)
	local xx,uu,sp,up = {},{},s1[1]
	local nx,nu = #sp, #s1
	for i=1,nu-1 do uu[i] = {} end
	jj = 1; for i=1,nx,2 do xx[jj] = sp[i]; jj = jj+1 end
	for i=2,nu do
		sp,up = s1[i],uu[i-1]
		jj = 1; for j=1,nx,2 do up[jj] = sp[j]; jj = jj+1 end
	end
	local s2,n2 = odefd(eqs,xx,uu); uu = {s1[1]}
	for i=2,nu do
		sp,up = s1[i],s2[i]
		uerr = odeerror({s1[1],sp},{s2[1],up}); uu[i] = uerr[2]
	end
	return s1,uu,math.max(n1,n2)
end
setfenv(odefde,getfenv(ode2bvfd));getfenv(odefde).odefd = odefd
getfenv(odefde).odeerror=odeerror

require"intp" -- needed for lintp()
xad = function(x,u)
	local ts,nx,nu,tss = {0}, #x, #u
	local um,umm,xn = {},{0},{}
	local xm,fac = x[nx]-x[1],.2
	for i=1,nu do um[i] = 0 end
	for j=2,nx do
		for i=1,nu do um[i] = um[i] + math.abs(u[i][j]-u[i][j-1]) end
	end -- um[] now contains maximum change across structure
	for j=2,nx do
		tss = 0 -- maximum percent change
		for i=1,nu do tss = math.max(tss, math.abs(u[i][j]-u[i][j-1])/um[i]) end
		tss = math.max(tss, fac*(x[j]-x[j-1])/xm) -- For flat regions of solution
		ts[j] = ts[j-1] + tss
	end
	tss = ts[nx] -- Maximum percentage change
	for j=2,nx do ts[j] = ts[j]*(nx-1)/tss end -- Normalize to 0 to nx-1
	for j=1,nx do xn[j] = lintp(ts,x,j-1) end
	for i=1,nu do 
		um[i] = {} -- Arrays for new solutions
		for j=1,nx do um[i][j] = lintp(x,u[i],xn[j]) end
	end
	return xn,um
end
		

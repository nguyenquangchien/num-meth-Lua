-- /* File odebvfd.lua */

odebv1fd = function(eqs,x,u)
	local s,up,upp,e,f,du = {},{},{},{},{},{}
	local fval,fu,fup,fupp,uppi,upi,ui,xi,duu
	local fctupp,fctup,fctu = FACT,FACT,FACT
	local feq,nx = eqs[1], #x
	local nxm1,nend = nx-1, NMAX
	for i=1,nx-1 do s[i] = x[i+1] - x[i] end
	
	bound = function(nb,nxb) -- Function to evaluate specified boundary conditions
		upi,ui,xi = up[nxb],u[nxb],x[nxb] -- boundary values
		fval = eqs[nb](ui,upi)
		fup = (eqs[nb](ui,upi+fctup) - fval)/fctup
		fu = (eqs[nb](ui+fctu,upi) - fval)/fctu
	end
	
	for j=1,NMAX do -- Major loop for iterative solution
		fctupp,fctup = 0,0
		for i=2,nx-1 do -- Calculate second derivative and first derivative for present itteration
			si,sisi,alfi = s[i],(s[i]+s[i-1]),s[i]/s[i-1]
			duu = 2*(u[i+1] - (alfi+1)*u[i] + alfi*u[i-1])/(si*sisi)
			upp[i],fctupp = duu, fctupp+abs(duu)
			duu = (u[i+1] + (alfi^2-1)*u[i] -alfi^2*u[i-1])/(alfi*sisi)
			up[i],fctup = duu, fctup+abs(duu)
		end
		alfi = s[2]/s[1] -- Handle end points as special cases, lower boundary
		upp[1] = upp[2] - (upp[3]-upp[2])/alfi 
		up[1] = (-u[3] + u[2]*(1+alfi)^2 - u[1]*(2/alfi+1))/(alfi*(s[2]+s[1]))
		alfi = s[nxm1]/s[nx-2] -- Upper boundary
		upp[nx] = upp[nxm1] + (upp[nxm1]-upp[nxm1-1])*alfi
		up[nx] = (u[nx]*(1+2*alfi) - u[nxm1]*(1+alfi)^2 + u[nx-2]*alfi^2)/(alfi*(s[nxm1]+s[nx-2]))
		fctupp,fctup = FACT*fctupp/nx, FACT*fctup/nx
		if fctupp==0 then fctupp = FACT end -- Try to protect against 0 values
		if fctup==0 then fctup = FACT end
		bound(2,1) -- Evaluate lower boundary conditions -- fval,fu,fup calculated for boundary
		duu = fup - fu*s[1] -- Determines first values of e and f
		e[2],f[2] = fup/duu, fval*s[1]/duu
		for i=2,nx-1 do -- Set up a,b,c,d arrays and save e and f as arrays -- do not need to save a,b,c,d
			uppi,upi,ui,xi = upp[i],up[i],u[i],x[i]
			fval = feq(xi,ui,upi,uppi)
			fupp = (feq(xi,ui,upi,uppi+fctupp) - fval)/fctupp -- Probe upp factor
			fup = (feq(xi,ui,upi+fctup,uppi) - fval)/fctup -- Probe up factor
			fu = (feq(xi,ui+fctu,upi,uppi) - fval)/fctu -- Probe u factor
			si,sisi,alfi = s[i],(s[i]+s[i-1]),s[i]/s[i-1]
			ai = (2*fupp/si - fup)*alfi/sisi
			bi = fu - (2*(alfi+1)*fupp/si - (alfi^2-1)*fup/alfi)/sisi
			ci,di = (2*fupp/s[i] + fup/alfi)/sisi, fval
			-- Forward reduction of tridiagonal system
			gi = 1/(ai*e[i] + bi)
			e[i+1],f[i+1] = -gi*ci, -gi*(di + ai*f[i])
		end
		
		bound(3,nx) -- Evaluate upper boundary conditions -- fval,fu,fup calculated for boundary
		-- Now back substitute for correction values
		du[nx] = -(fval*s[nxm1] - f[nx]*fup)/(fu*s[nxm1] + fup*(1-e[nx])) -- End point correction from boundary conditions
		for i=nx,2,-1 do  -- Calculate correction values at all points using e and h arrays
			du[i-1] = e[i]*du[i] + f[i]
		end
		-- Now update solution and check for desired accuracy
		cmax,umax,imax,fctu = 0,0,0,0
		for i=1,nx do
			c1 = abs(du[i]); if c1>umax then umax = c1 end
			c2 = abs(u[i]) + c1; if c2~=0.0 then c1=c1/c2 end
			if c1>cmax then cmax=c1; imax=i end
			u[i] = u[i] + du[i] -- Add corection to previous solution
			fctu = fctu+abs(u[i])
		end
		fctu = fctu/nx; if fctu==0 then fctu = FACT end
		if nprint~=0 then 
			printf("Iteration number %i, Maximum  relative, absolute correction = %e, %e at %i\n",j,cmax,umax,imax) 
			io.flush()
		end
		if cmax<ERROR then nend=j; break end
		if FABS~=0 then if umax<fctu*ERROR then nend=j; break end end
		fctu = FACT*fctu
	end
	return nend, cmax, umax -- Solution values returned in u calling array
end
setfenv(odebv1fd,{abs=math.abs,ERROR=1.e-6,FABS=1,FACT=1.e-3,EPXX=1.e-14,NMAX=50,nprint=0,io=io,printf=printf,
table=table})

odeerror = function(si1,si2) -- ODE error estimator for two solutions 
	require"intp"
	local err,ny,s1,s2 = {},#si1
	local nd1,nd2 = #si1[1], #si2[1]
	local fac = 3
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

ode1fd = function(eqs,x,u) -- Simple FD interface
	local nx,nu,x1,x2 = #x,#u
	local u1,u2 = u[1],u[nu]
	if nx<=3 then -- Generate linear step distribution
		x1,x2,nx = x[1],x[2],x[3] or 2001 
		nx = 2*math.floor(nx/2)+1 -- Make odd integer
		local dx = (x2-x1)/(nx-1) -- Even number of intervals
		for i=1,nx do x[i] = x1+(i-1)*dx end
	end -- If u array not specified, generate below
	if nu~=nx then for i=1,nx do u[i] = u1+(i-1)/(nx-1)*(u2-u1) end end
	if type(eqs)=='function' then eqs = {eqs} end -- Boundary functions not given
	eqs[2]  = eqs[2] or function(u,up) return u - u1 end
	eqs[3] = eqs[3] or function(u,up) return u - u2 end
	return {x,u},odebv1fd(eqs,x,u)
end	
ode1fde = function(eqs,x,u) -- Simple FD interface with error estimate return
	local s1,n1 = ode1fd(eqs,x,u)
	local xx,uu,j = {},{},1
	for i=1,#s1[1],2 do 
		xx[j],uu[j] = s1[1][i],s1[2][i]
		j = j+1
	end
	local s2,n2 = ode1fd(eqs,xx,uu)
	return s1,odeerror(s1,s2),math.max(n1,n2)
end

xgp = function(xmin,xmax,a,n)
	n = n or 2001; n = 2*math.floor(n/2)+1 -- Make odd integer
	local s = xmax - xmin
	require'newton'
	local fgn = function(r) return r-(s*(r-1)/a+1)^(1/n) end
	local r = newton(fgn,2)
	local step,x = a,{}; x[1] = xmin
	for i=2,n-1 do
		step = step*r
		x[i] = x[i-1] + step
	end
	x[n] = xmax; return x
end
xlg = function(xmin,xmax,a,n)
	n = n or 2001; n = 2*math.floor(n/2)+1 -- Make odd integer
	local fac = ((xmax-xmin)/a)^(1/(n-2))
	local step,x = a,{}; x[1],x[2] = xmin, step
	local i=2
	while 1 do
		step = step*fac
		i = i+1; x[i] = step
		if 2*x[i]-x[i-1]>xmax then break end
	end
	x[n] = xmax; return x
end
xll = function(xmin,xmax,a,n)
	n = n or 2001; n = 2*math.floor(n/2)+1 -- Make odd integer
	local dx = xmax-xmin
	local fac = math.log(dx/a)/math.log(n-1)
	local x = {}; x[1] = xmin
	for i=2,n-1 do x[i] = xmin + a*(i-1)^fac end
	x[n] = xmax; return x
end

	

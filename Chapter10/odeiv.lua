-- /* File odeiv.lua */
-- Programs to integrate first order differential equations
require"nsolv"; local atend

odebiv = function(feqs,tvals,u,up,upp)
	local j, neq, t, h, h2,hs,hx -- Local variables for function
	local sol, yn, jfirst = {}, {}, 0
	local nit,nitt = 0,0 -- Number of iterations
	local neq = #u
	local tmin,tmax,ntval = tvals[1],tvals[2],tvals[3]
	-- Function to calculate next time values using nsolv()
	local fnext = function(eqs,u) 
		for m=1,neq do
			up[m] = (u[m] - yn[m])/h2 -- h2 = h/2 for TP rule, h for BD
		end
		feqs(eqs,t,u,up) -- Now call user defined differential equations
	end
	-- Use input value for number of intervals or set at default value
	up = up or {}
	for m=1,neq+1 do -- Array to return solution values with t value first
		sol[m] = {}
		if m==1 then sol[m][1] = tmin else sol[m][1] = u[m-1]  end
	end
	t = tmin -- Initial t value
	hs = (tmax - t)/ntval -- Equal increments in t used
	-- If initial derivative not specified, use BD for first 4 points
	if #up~=neq then -- up array not input
		for m=1,neq do up[m] = 0 end
		jfirst,h,h2,hx = 0,hs/4,hs/4,0 -- Set to BD parameters
	else jfirst,h,h2,hx = 4,hs,hs/2,hs/2 end -- Set to TP parameters
	for k=1,ntval do -- Main loop for incrementing independent variable (t)
		repeat -- Use backwards differencing with 4 sub intervals of size h/4
			jfirst = jfirst+1
			-- Set up yn array and get ready for next solution value
			for m=1,neq do
				yn[m] = u[m] + hx*up[m] -- hx = 0 or h/2
				u[m] = u[m] + h*up[m] -- Predicted value of u array
			end
			t = t + h -- Now increment t to next t value
			-- Calculate new u values at next time step, u returns new values
			nitt = nsolv(fnext,u,step) 
			if nitt>nit then nit = nitt end -- Monitor maximun # of iterations
			-- New Derivative values, same function as in fnext
			for m=1,neq do up[m] = (u[m] - yn[m])/h2 end
		until jfirst>=4 -- End of first interval repeating using BD
		if k==1 then h,h2,hx = hs,hs/2,hs/2 end -- Set to TP parameters
		sol[1][k+1] = t -- Save calculated values
		for m=1,neq do sol[m+1][k+1] = u[m]	end
	end -- End of main loop on t, now return solution array
	sol[1][ntval+1] = tmax; return sol,nit -- Solution and maximim number of iterations
end -- End of odebiv
setfenv(odebiv,{nsolv=nsolv,step=nil})
--[[ /****************************************************************************
	odebiv solves a system of first order differential equations of the form:
Fm(up,u,t) = 0. where m is the number of functions and equations to be solved.  The
defining equations may consist of any nonlinear or linear combinations of the variables
and derivatives and the functions may have mixed variables and derivatives.  On input, 
feqs is the name of a function returning the zero valued functions, u is the initial values
of the variables, up is the initial derivatives of the variables and tvals is an array of
{tmin,tmax,ntval} specifying the range of time for the solution. On output, u and up
the last values of variables and derivatives, so that the function can be called again and
the solution will pick up where it left off if an increased time is specified. 
*******************************************************************************/]]

odeiv = function(feqs,tvals,u,up,upp) -- Multi Time Step Solver
	local sa,sb,tl,ns = {}
	local ntp,ni,nit 
	local NMAX,ND = getfenv(nsolv).NMAX,1000 -- 1000 default intervals
	up,upp = up or {},upp or {}
	if type(tvals)=='number' then tvals = {0,tvals,0} end
	local j = #tvals
	if j==1 then tvals = {0,{tvals[1]},2*ND} end
	if j==2 then tvals[3] = 2*ND end
	if type(tvals[2])=='number' then tvals[2] = {tvals[2]} end
	ntp = #tvals[2]
	if type(tvals[3])=='table' then ns = tvals[3][1] else ns = tvals[3] end
	nit = 0
	for i=1,ntp do
		if i>1 then tl = tvals[2][i-1] else tl = tvals[1] end
		if type(tvals[3])=='table' then ns = tvals[3][i] or ns end
		sb,ni = odebiv(feqs,{tl,tvals[2][i],ns},u,up,upp)
		if ni==NMAX then 
			print("Error: Maximum number of iterations exceeded in nsolv")
			print("Results may not be accurate!"); print(tl,tvals[2][i],ns)
		end
		if ni>nit then nit = ni end
		if i>1 then sa = atend(sa,sb) else sa = sb end
	end
	return sa,nit
end

odeivqs = function(feqs,tvals,u,up) -- ODE Quick Scan function
	local NPTS,NPS = 20,2
	local ttvals,nl,nu,fact = {},10
	local nt,j = #tvals,1
	if nt<2 then print('Error, must specify two times in obeivqs'); return end
	NPTS = floor(((tvals[3] or NPTS)+1)/2)
	nl = 10^(floor(log10(tvals[2][1])))
	nu = 10^(ceil(log10(tvals[2][2])))*1.000001
	fact = 10^(1/NPTS)
	ttvals[1],ttvals[2],ttvals[3],nl = nl,{},NPS,nl*fact
	while nl<= nu do ttvals[2][j],nl,j = nl,nl*fact,j+1 end
	ttvals[2][#ttvals[2]] = 10^(ceil(log10(tvals[2][2]))) -- Exact value at end
	u,up = u or {},up or {} -- Linear steqs for first interval
	odeiv(feqs,{tvals[1],ttvals[1],NPTS},u,up) -- NPTS points
	return odeiv(feqs,ttvals,u,up)
end
setfenv(odeivqs,{floor=math.floor,ceil=math.ceil,log10=math.log10,odeiv=odeiv})

odebrk = function(feqs,tvals,u,up,upp) -- Basic Runge-Kutta ODE integration code
	local j, neq, t, h, h2 -- Local variables for function
	local sol, m1,m2,m3,m4,u1 = {}, {}, {}, {}, {}, {}
	local nit,nitt = 0,0 -- Number of iterations
	local neq = #u -- Number of equations 
	local tmin,tmax,ntval = tvals[1],tvals[2],tvals[3]
	
	fderiv = function(eqs,up) -- Function for calculating derivatives
		feqs(eqs,t,u1,up)
	end
	-- Use input value for number of intervals or set at default value
	if type(up)~='table' then up = {} end
	for m=1,neq+1 do -- Array to return solution values with t value first
		sol[m] = {}
		if m==1 then sol[m][1] = tmin else	sol[m][1] = u[m-1]  end
	end
	t = tmin; hs = (tmax - t)/ntval -- Equal increments in t used
	for m=1,neq do up[m] = 0 end
	h,h2 = hs,hs/2 -- Set to RK parameters
	for k=1,ntval do -- Main loop for incrementing independent variable (t)
		for m=1,neq do u1[m] = u[m] end
		nitt = nsolv(fderiv,up,step) -- Update up values
		for m=1,neq do m1[m] = h*up[m]; u1[m] = u[m] + m1[m]/2 end
		t = t+h2; nitt = nsolv(fderiv,up,step)+nitt -- next up value
		for m=1,neq do m2[m] = h*up[m]; u1[m] = u[m] + m2[m]/2 end
		nitt = nsolv(fderiv,up,step)+nitt -- next up value
		for m=1,neq do m3[m] = h*up[m]; u1[m] = u[m] +m3[m] end
		t = t+h2; nitt = nsolv(fderiv,up,step)+nitt -- next up value
		for m=1,neq do m4[m] = h*up[m] end
		for m=1,neq do u[m] = u[m] + (m1[m]+2*m2[m]+2*m3[m]+m4[m])/6 end
		sol[1][k+1] = t -- Save calculated values
		for m=1,neq do sol[m+1][k+1] = u[m] end
		if nitt>nit then nit = nitt end -- Monitor maximun # of iterations
	end -- End of main loop on t, now return solution array
	return sol,nit/4 -- Solution and maximim number of iterations
end -- End of odebrk
setfenv(odebrk,{type=type,nsolv=nsolv,step=nil})	

odeivs = function(feqs,tvals,u,up) -- Adaptive Step Size Solver
	local ttvals,sa,upp,t,s1,s2 = {},{},{}
	local k,ni,nit,err,fac,relerr,relmin,t1,t2,t3,h,h2,hmax,hmin,tmax
	local neq,NMAX,abs = #u, getfenv(nsolv).NMAX, math.abs
	local u1,u2,up1,up2,upp1,upp2 ={},{},{},{},{},{} -- Saved arrays
	local NTMIN,NTMAX,TFAC,HFAC,FAC,RELERR,fe,fd,tt,nhmax = 1000,25000,1.e-6,1.e-12,.8,1.e-5,.5,1,1+1.e-12,0
	if odebiv==odebrk then fe,fd = 0.25,0.2 end -- Set parameters for RK algorithm
	up = up or {} -- Below is for different tvals formats
	if type(tvals)=='number' then tvals = {0,tvals} end
	if #tvals==1 then tvals = {0,tvals[1]} end
	if type(tvals[2])=='number' then tvals[2] = {tvals[2]} end
	t,tmax = tvals[1],tvals[2][1]
	hmin,hmax = tvals[2][2] or (tmax-t)*HFAC, tvals[2][3] or (tmax-t)/NTMIN 
	relerr = tvals[3] or RELERR*neq; relmin = relerr/5
	nit,k,h = 0,1,(tmax-t)*TFAC; h2 = h/2 -- Use TFAC for initial step size
	for i=1,neq+1 do sa[i] = {} end; sa[1][1] = t -- Set initial solution values
	for i=1,neq do t1,t2 = u[i],up[i]; sa[i+1][1],u1[i],u2[i],up1[i],up2[i] = t1,t1,t1,t2,t2 end
	while 1 do -- Major time step loop
		while 1 do -- Adjust step size downward until local relative error condition met
			ttvals[1],ttvals[2],ttvals[3] = t,t+h,1
			s1,nx = odebiv(feqs,ttvals,u1,up1,upp1) -- One step solution
			ttvals[3] = 2; s2,ni = odebiv(feqs,ttvals,u2,up2,upp2) -- Repeat using two steps
			err = 0 -- Evaluate maximum relative error
			for i=1,neq do
				fac =  fd*abs(u2[i]-u1[i])/(abs(u1[i]) + abs(u2[i]))
				if fac>err then err = fac end
			end
			if h==hmin then break end -- Just accept, nothing else useful to do 
			if err<relerr then break end -- Accept error met
			if nx==NMAX then -- Didn't converge try again at half step size
				if h==hmax then hmax = hmax/2 end
				h,h2 = h/2,h2/2 
			elseif err==1 then h,h2 = h/2,h2/2 -- Try again at half step size
			else h = (relerr/err)^fe*FAC*h; h2 = h/2 end -- Adjust step size downward and try again 
			if abs(h)<abs(hmin) then h,h2 = hmin,hmin/2 end -- Use mimimum step size
			for i=1,neq do t1,t2,t3 = u[i],up[i],upp[i]
			u1[i],u2[i],up1[i],up2[i],upp1[i],upp2[i] = t1,t1,t2,t2,t3,t3 end -- Set to old values 
		end -- loop back if relerr criteria not met
		if ni==NMAX and err>relerr then -- Print warning message
			print("Error at t =" ,t," : Maximum number of iterations exceeded in nsolv")
			print("     Results are probably not accurate!")
		end
		if ni>nit then nit = ni end
		for i=1,2 do -- Save valid solution values contained in s2 array -- 2 time points
			k,t = k+1,t+h2; sa[1][k] = t; for j=2,neq+1 do sa[j][k] = s2[j][i+1] end 
		end
		if k>NTMAX then -- Limit solution to NTMAX data points (100,000 default)
			print("Number of adaptive data points exceeds ",NTMAX)
			print("     Best effort at solution returned!"); break
		end
		for i=1,neq do t1,t2,t3 = u2[i],up2[i],upp2[i]; u[i],u1[i],up[i],up1[i],upp[i],upp1[i] = t1,t1,t2,t2,t3,t3 end
		if h>0 then if t>=tmax then break end -- Exit if finished, else reset parameters
		elseif t<=tmax then break end -- Negative time step
		if err<relmin then h = h*1.4 end -- Adjust step size upward by factor of 1.4
		if abs(h)>abs(hmax) then 
			h = hmax; nhmax = nhmax+1
			if nhmax>10 then nhmax,hmax = 0,hmax*1.4 end 
		else nhmax = 0 end 
		if h>0 then if t+h+h>tmax then h = tmax*tt - t  end 
		elseif t+h+h<tmax then h = tmax*tt - t end 
		h2 = h/2 -- Finalize h and h/2
	end -- loop back for next time step
	sa[1][#sa[1]] = tmax; return sa,nit -- Set limit to exactly tmax
end	

odeb12 = function(feqs,tvals,u,up,upp)
	local j, neq, t, h, h2,h2sq,hs,hx,hy,hz -- Local variables for function
	local sol,ynn,yn,jfirst = {},{},{},0
	local nit,nitt = 0,0 -- Number of iterations
	local neq = #u
	local tmin,tmax,ntval = tvals[1],tvals[2],tvals[3]
	-- Function to calculate next time values using nsolv()
	local fnext = function(eqs,u) 
		for m=1,neq do
			up[m] = (u[m] - yn[m])/h2 -- h2 = h/2 for TP rule, h for BD
			upp[m] = (u[m] - ynn[m])/h2sq -- h2sq = (h/2)^2 for TP, (h)^2 for BD
		end
		feqs(eqs,t,u,up,upp) -- Now call user defined differential equations
	end
	-- Use input value for number of intervals or set at default value
	up,upp = up or {}, upp or {}
	for m=1,neq+1 do -- Array to return solution values with t value first in array
		sol[m] = {}
		if m==1 then sol[m][1] = tmin else sol[m][1] = u[m-1]  end
	end
	t = tmin -- Initial t value
	hs = (tmax - t)/ntval -- Equal increments in t used, no adjusting step size
	-- If initial derivative not specified, use Backwards differencing for first 4 points
	if #up~=neq then for m=1,neq do up[m] = up[m] or 0 end end -- Complete initial derivatives 
	if #upp~=neq then for m=1,neq do upp[m] = 0 end
		jfirst,h = 0,0.25*hs; h2,h2sq,hy,hx,hz = h,h*h,h,0,0 -- Set to BD parameters
	else jfirst,h = 4,hs; h2,h2sq,hy,hx,hz = hs/2,h*h/4,h,h/2,h*h/4 end -- TP parameters
	for k=1,ntval do -- Main loop for incrementing independent variable (t)
		repeat -- Use backwards differencing for first interval with 4 sub intervals of size h/4
			jfirst = jfirst+1
			-- Set up yn, and ynn arrays and get ready to solve equations
			for m=1,neq do
				yn[m] = u[m] + hx*up[m] -- hx = 0 or h/2
				ynn[m] = u[m] + hy*up[m] + hz*upp[m] -- hy=h, hz=0 or (h/2)^2
				u[m] = u[m] + h*(up[m] + 0.5*h*upp[m]) -- Predicted value of u array
			end
			t = t + h -- Now increment t to next t value
			-- Calculate new u values at next time step, new values are returned in u array
			nitt = nsolv(fnext,u,step) 
			if nitt>nit then nit = nitt end -- Monitor maximun number of iterations
			-- New derivative values, same function as in fnext 
			for m=1,neq do up[m],upp[m] = (u[m] - yn[m])/h2,(u[m] - ynn[m])/h2sq end
		until jfirst>=4 -- End of first interval repeat using Backwards difference
		if k==1 then h = hs; h2,h2sq,hy,hx,hz = h/2,h*h/4,h,h/2,h*h/4 end -- TP parms.
		sol[1][k+1] = t; for m=1,neq do sol[m+1][k+1] = u[m]	end -- Save values
	end -- End of main loop on t, now return solution array
	sol[1][ntval+1] = tmax; return sol,nit -- Solution and maximim number of iterations
end -- End of odeb12
setfenv(odeb12,{nsolv=nsolv,step=nil})
--[[ /****************************************************************************
	odeb12() solves a system of second order differential equations of the form:
Fm(t,u,up,upp) = 0. where m is the number of functions and equations to be solved.  The
defining equations may consist of any nonlinear or linear combinations of the variables
and derivatives and the functions may have mixed variables and derivatives.  On input, 
feqs is the name of a function returning the zero valued functions, u is the initial values
of the variables, upp, and up are the initial derivatives of the variables and tvals is an array of
{tmin,tmax,ntval} specifying the range of time for the solution. On output, u, upp and up
are the last values of variables and derivatives, so that the function can be called again 
and the solution will pick up where it left off if an increased time is specified. 
*******************************************************************************/]]

odeerror = function(si1,si2) -- ODE error estimator for two solutions 
	require"intp"
	local err,ny,s1,s2,sx,sy = {},#si1
	local nd1,nd2 = #si1[1], #si2[1]
	local fac,j,rev = 3,1,false
	if getfenv(odeiv).odebiv==odebrk or odebiv==odebrk then fac = 15 end
	for i=1,ny do err[i] = {} end
	if nd1>nd2 then s1,s2 = si1,si2 
	else nd1,nd2,s1,s2 = nd2,nd1,si2,si1 end
	if nd1~=(nd2-1)*2+1 then print('Error: Incorrect files input to odeerror, nd1,nd2 = ',
		nd1,nd2) end
	err[1] = s1[1] -- set x values for errors
	if s1[1][nd1]<s1[1][1] then rev = true end -- reversed direction
	if rev then sx = {}; for i=1,nd2 do sx[i] = s2[1][nd2+1-i] end 
	else sx = s2[1] end
	for k=2,ny do  -- interpolate for error values
		if rev then sy = {}; for i=1,nd2 do sy[i] = s2[k][nd2+1-i] end
		else sy = s2[k] end
		for i=1,nd1 do err[k][i],j = (intp(sx,sy,s1[1][i]) - s1[k][i])/fac, j+1 end  
		j = 1 
	end
	for k=2,ny do for i=2,nd1-1,2 do -- Interpolate for even data points
		err[k][i] = err[k][i-1]+(err[k][i+1]-err[k][i-1])*(err[1][i+1]-err[1][i])/(err[1][i+1]-err[1][i-1])
	end end
	return err
end
errstat = function(yt) -- ODE error statistics
	require"prob"
	local nd = #yt
	local ym,pm = 0,0
	local mean,std = stats(yt)
	for i=1,nd do
		if math.abs(yt[i]) > ym then ym,pm = math.abs(yt[i]),i end
	end
	return std,ym,pm,mean
end
maxvalue = function(yt)
	local nv,nd = #yt,#yt[1]
	local mx,tmp = {}
	for i=1,nv do
		tmp = 0
		for j=1,nd do
			if math.abs(yt[i][j])>tmp then tmp=math.abs(yt[i][j]) end
		end
		mx[i] = tmp
	end
	return mx
end
odeivqse = function(feqs,tvals,u,up) -- ODE Quick Scan with error function
	up = up or {}
	local ul,ulp = {},{} 
	for i=1,#u do ul[i],ulp[i] = u[i],up[i] end
	local s1,n1 = odeivqs(feqs,tvals,u,up)
	tvals[3] = tvals[3]/2
	local s2,n2 = odeivqs(feqs,tvals,ul,ulp) 
	return s1,odeerror(s1,s2),max(n1,n2),s2
end
setfenv(odeivqse,{max=math.max,odeivqs=odeivqs,odeerror=odeerror})
odeive = function(feqs,tvals,u,up) -- ODE Multi-Step TP solver with error estimation
	up = up or {}
	local ul,ulp = {},{}
	local ttvals = {tvals[1],tvals[2],{}}
	for i=1,#u do ul[i],ulp[i] = u[i],up[i] end
	local s1,n1 = odeiv(feqs,tvals,u,up)
	if type(tvals[3])=='number' then ttvals[3][1] = tvals[3]/2
	else for i=1,#tvals[3] do ttvals[3][i] = tvals[3][i]/2 end end
	local s2,n2 = odeiv(feqs,ttvals,ul,ulp)
	return s1,odeerror(s1,s2),max(n1,n2),s2
end	
setfenv(odeive,{max=math.max,type=type,odeiv=odeiv,odeerror=odeerror})
odeivse =  function(feqs,tvals,u,up) -- Adaptive Step Size Solver with error estimation
	up = up or {}
	local u1,up1,k = {},{},1
	local ttvals = {tvals[1],{},NS}
	for i=1,#u do u1[i],up1[i] = u[i],up[i] end -- Save initial values
	local s1,n1 = odeivs(feqs,tvals,u,up) -- Real solution at 2h step sizes
	for i=3,#s1[1],2 do ttvals[2][k],k = s1[1][i],k+1 end
	local s2,n2 = odeiv(feqs,ttvals,u1,up1) -- Comparison at h step sizes
	local err,j,nk,nd1,nd2 = odeerror(s1,s2),1,#s1,#s1[1],#s2[1]
	local er2 = {}
	for k=1,#s1 do 
		er2[k],s1[k],j = {},{},1 
		for i=1,nd2,2 do s1[k][j],er2[k][j],j = s2[k][i],err[k][i], j+1 end
	end 
	return s1,er2,min(n1,n2),s2 -- Return solution and errors
end
setfenv(odeivse,{odeivs=odeivs,odeiv=odeiv,odeerror=odeerror,min=math.min,NS=4})
	
atend = function(a,b)
	local nc,m1,m2 = #a,#a[1],#b[1]
	for i=1,nc do
		for j=2,m2 do a[i][j-1+m1] = b[i][j] end
	end
	return a
end
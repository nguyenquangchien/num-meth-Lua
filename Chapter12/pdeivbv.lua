-- /* File pdeivbv.lua */
-- Programs to solve partial differential equations of the form initial value, boundary value in two variables
require"odefd"; require"intp"
pdebivbv = function(eqsub,tvals, x,u,ut,utt) -- Basic time step solver for partial differential equation
	local neq,nx,nut,nutt,t,h,h2,h2sq,hs,att,btt,at = #u, #x, #ut, #utt
	local nit,nitt,eqsu = 0,0, eqsub[1] -- Number of iterations
	local fu,xt,xtt,uti,utti = {},{},{},{},{}
	eqtr = function(fu,x,u,ux,uxx,i) -- Map of time differentials into spatial values
		for m=1,neq do 
			uti[m] = (u[m] - xt[m][i])/h2 -- h2 = h/2 for trapezoidal rule, h for Backwards differencing
			utti[m] = (u[m] - xtt[m][i])/h2sq -- Second derivative
		end; eqsu(fu,x,t,u,ux,uxx,uti,utti) -- Call defining equation set, returning fu values
	end
	eqbl = function(fu,u,ux) -- map left boundary condition
		for m=1,neq do
			uti[m] = (u[m] - xt[m][1])/h2 -- Left time derivative
		end; eqsub[2](fu,u,ux,t,uti)
	end
	eqbr = function(fu,u,ux) -- map right boundary condition
		for m=1,neq do
			uti[m] = (u[m] - xt[m][nx])/h2 -- Left time derivative
		end; eqsub[3](fu,u,ux,t,uti)
	end
	local eqsx = {eqtr,eqbl,eqbr}
	local tmin,tmax,ntval = tvals[1],tvals[2],tvals[3]
	t,hs = tmin, (tmax-tmin)/ntval  -- Initial t value
	if nutt~=neq then for m=1,neq do utt[m] = {} end end
	if nut~=neq then for m=1,neq do ut[m] = {} end end
	for m=1,neq do
		if #utt[m]~=nx then nutt = 0; for k=1,nx do utt[m][k] = 0 end end
		if #ut[m]~=nx then nut = 0; for k=1,nx do ut[m][k] = 0 end end
	end
	if nutt~=neq then -- utt aray not input, use backwards difference initially
		jfirst,h,h2,h2sq,att,at,btt=0,hs/4,hs/4,(hs/4)^2,hs/4,0,0
	else jfirst,h,h2,h2sq,att,at,btt=4,hs,hs/2,(hs/2)^2,hs,hs/2,(hs/2)^2 end
	for m=1,neq do xt[m],xtt[m],fu[m] = {},{},0 end 
	for k=1,ntval do -- Beginning of major time loop -- Heart of solution
		repeat -- Use Backwards differencing for first interval with 4 sub intervals
			jfirst = jfirst+1
			for i=1,nx do -- Set up xx arrays for use in replacing time derivatives
				for m=1,neq do
					xt[m][i] = u[m][i] + at*ut[m][i]
					xtt[m][i] = u[m][i] + att*ut[m][i] + btt*utt[m][i]
				end
			end -- These are used by eqtr() function when called
			t = t + h -- Now increment t to next t value
			nitt = ode2bvfd(eqsx,x,u) -- Calculate new u values, returned in u array
			if nitt>nit then nit = nitt end -- Monitor maximum number of iterations
			if nprint~=0 then 
				printf('Time = %e Number of iterations in pdebivbv = %d \n\n' ,t,nitt)
				io.flush()
			end
			for i=1,nx do -- Now calculate derivative values at new time point
				for m=1,neq do
					ut[m][i] = (u[m][i] - xt[m][i])/h2 
					utt[m][i] = (u[m][i] - xtt[m][i])/h2sq
				end
			end
		until jfirst>=4 -- End of first interval with Backwards differencing
		if k==1 then h,h2,h2sq,att,at,btt=hs,hs/2,(hs/2)^2,hs,hs/2,(hs/2)^2 end
	end -- End of major loop, Go back and do for another time
	return nit, tmax, tmin, ntval
end	-- End of function	
setfenv(pdebivbv,{ode2bvfd=ode2bvfd,printf=printf,type=type,nprint=0,io=io,plot=plot})

local usave,utsave
pdeivbv = function(eqsub,tvals,x,u,ut,utt) -- Simple initial-value, boundary-value interface  
	local tmin,ntkeep,ntval,dtkeep,nit
	local usol = {x} -- Array to hold solutions
	if type(u[1])~='table' then u = {u} end -- For one dimensional case
	local neq,nx = #u, #x -- Number of solution variables
	ntkeep,ntval = tvals[3] or NKEEP, tvals[4] or NTVAL -- Set number of t intervals
	if type(tvals[2])=='table' then ntkeep,ntval = #tvals[2], tvals[3] or NTVAL end
	if type(tvals[2])=='number' then -- Simple input format
		if type(tvals[3])~='table' then tvals[3] = {} end --Ensure table 
		dtkeep = (tvals[2] - tvals[1])/ntkeep; tvals[2] = {}
		for j=1,ntkeep do -- Set up time arrays for simple input
			tvals[2][j] = tvals[1] + j*dtkeep
			tvals[3][j] = tvals[3][j] or ntval
		end
	elseif type(tvals[3])=='number' then 
		tvals[3] = {};
		for j=1,ntkeep do tvals[3][j] = ntval end 
	end
	if type(utt)~='table' then utt = {} end -- Ensure tables for derivitives 
	if type(ut)~='table' then ut = {} end
	local nutt,nut = #utt, #ut
	if nut~=neq then for m=1,neq do ut[m] = {} end end
	for m=1,neq do
		if #ut[m]~=nx then for k=1,nx do ut[m][k] = 0 end end
	end
	tmin = tvals[1]
	for i=1,ntkeep do -- Step over time increments
		usol = usave(usol,u) -- Save solution values
		nit = pdebivbv(eqsub,{tmin,tvals[2][i],tvals[3][i]},x,u,ut,utt) -- Now solve it !
		tmin = tvals[2][i] -- Increment initial time value
		if nprint~=0 then 
			printf('Time = %e Number of iterations in pdeivbv = %d \n\n' ,tmin,nit)
			io.flush()
		end
		if nprint==2 then for j=1,neq do plot(x,u[j]) end end
	end
	usol = usave(usol,u) -- Save final values
	return usol
end
setfenv(pdeivbv,{pdebivbv=pdebivbv,printf=printf,plot=plot,type=type,
nprint=0,io=io,NKEEP=10,NTVAL=10,usave=usave})

pdeivbvt = function(eqsub,tvals,x,u,ut,utt) -- Simple initial-value, boundary-value interface 
	if tvals[5]==nil then return pdeivbv(eqsub,tvals,x,u,ut,utt) end
	local tmin,ntkeep,ntval,dtkeep,nit,nitt
	local tmax,tloc,dt
	local nprint = nprint or getfenv(pdeivbv).nprint
	local usol, sta = {x}, {} -- Array to hold solutions
	if type(u[1])~='table' then u = {u} end -- For one dimensional case
	local neq,nx = #u, #x -- Number of solution variables
	ntkeep,ntval = tvals[3] or NKEEP, tvals[4] or NTVAL -- Set number of t intervals
	if type(tvals[2])=='table' then ntkeep,ntval = #tvals[2], tvals[3] or NTVAL end
	if type(tvals[2])=='number' then -- Simple input format
		if type(tvals[3])~='table' then tvals[3] = {} end --Ensure table 
		dtkeep = (tvals[2] - tvals[1])/ntkeep; tvals[2] = {}
		for j=1,ntkeep do -- Set up time arrays for simple input
			tvals[2][j] = tvals[1] + j*dtkeep
			tvals[3][j] = tvals[3][j] or ntval
		end
	elseif type(tvals[3])=='number' then 
		tvals[3] = {};
		for j=1,ntkeep do tvals[3][j] = ntval end 
	end
	if type(utt)~='table' then utt = {} end -- Ensure tables for derivitives 
	if type(ut)~='table' then ut = {} end
	local nutt,nut = #utt, #ut
	if nut~=neq then for m=1,neq do ut[m] = {} end end
	for m=1,neq do
		if #ut[m]~=nx then for k=1,nx do ut[m][k] = 0 end end
	end
	if type(tvals[5])=='number' then
		nts = tvals[5]; tvals[5] = {}; tmp = x[nx] - x[1]
		for i=1,nts+1 do tvals[5][i] = (i-1)*tmp/nts + x[1] end
	end
	tmin = tvals[1]
	sta = utsave(sta,u,x,tmin,tvals[5]) -- Save initial time values
	for i=1,ntkeep do -- Step over time increments
		tmax = tvals[2][i]; ntval = tvals[3][i]
		dt = (tmax - tmin)/ntval; tloc = tmin + dt
		usol = usave(usol,u) -- Save solution values
		nitt = 0
		for j=1,ntval do -- step over each spatial save interval
			nit = pdebivbv(eqsub,{tmin,tloc,1},x,u,ut,utt) -- one time interval
			if nit>nitt then nitt = nit end
			sta = utsave(sta,u,x,tloc,tvals[5]) -- Save time data
			tmin = tloc; tloc = tloc + dt
		end -- Now have one spatial save intveral
		tmin = tvals[2][i] -- Increment initial time value
		if nprint~=0 then 
			printf('Time = %e Number of iterations in pdeivbvt = %d \n\n' ,tmin,nitt)
			io.flush()
		end
		if nprint==2 then for j=1,neq do plot(x,u[j]) end end
	end
	usol = usave(usol,u) -- Save final values
	return usol, sta
end
setfenv(pdeivbvt,{pdebivbv=pdebivbv,printf=printf,plot=plot,type=type,getfenv=getfenv,
nprint,io=io,NKEEP=10,NTVAL=10,usave=usave,pdeivbv=pdeivbv,utsave=utsave})

pdeivbvqs = function(eqsub,tvals,x,u,ut,utt) -- Quick scan log time solver
	local nps,npts -- save per decade, additional steps per save
	local ttvals,nl,nu,fact = {}
	local nt,neq,j = #tvals, #u, 0
	local nitt,nit,sol = 0 
	local nprint = nprint or getfenv(pdeivbv).nprint
	ut = ut or {}; utt = utt or {}
	if type(u[1])~='table' then u = {u} end -- For single equation
	if nt<2 then printf('Error, must specify two times in pdeivbvqs')
		return end
	nps,npts = math.floor(tvals[3] or NPS),math.floor(tvals[4] or NPTS)
	fact = 10^(1/(nps*npts)) -- Factor between steps
	nl = 10^(math.floor(math.log10(tvals[2][1])))
	nu = 10^(math.ceil(math.log10(tvals[2][2])))*1.00000000001/fact
	sol = pdeivbv(eqsub,{tvals[1],nl,1,nps*npts},x,u,ut,utt) -- Initial interval
	while nl<=nu do -- Step over log time spacings
		nit = pdebivbv(eqsub,{nl,nl*fact,1},x,u,ut,utt)
		if nit>nitt then nitt = nit end
		nl = nl*fact; j = j+1
		if j==npts then 
			sol = usave(sol,u); j = 0 -- Save solutions
			if nprint~=0 then 
				printf('Time = %e Number of iterations in pdeivbvqs = %d \n\n' ,nl,nitt)
				io.flush(); nitt = 0
			end
		end
	end
	return sol
end
setfenv(pdeivbvqs,{pdebivbv=pdebivbv,pdeivbv=pdeivbv,printf=printf,type=type,
math=math,nprint,io=io,NPS=1,NPTS=10,usave=usave,getfenv=getfenv})

pdeivbvqst = function(eqsub,tvals,x,u,ut,utt) -- Quick scan log time solver
	if tvals[5]==nil then return pdeivbvqs(eqsub,tvals,x,u,ut,utt) end
	local nps,npts -- save per decade, additional steps per save
	local ttvals,nl,nu,fact = {}
	local nt,neq,j = #tvals, #u, 0
	local nitt,nit,sol = 0,0,{x}
	local nprint = nprint or getfenv(pdeivbv).nprint
	local nx,sta,k,nts,tmp = #x, {}, 1
	ut = ut or {}; utt = utt or {}
	if type(u[1])~='table' then u = {u} end -- For single equation
	if nt<2 then printf('Error, must specify two times in pdeivbvqs')
		return end
	nps,npts = math.floor(tvals[3] or NPS),math.floor(tvals[4] or NPTS)
	if type(tvals[5])=='number' then
		nts = tvals[5]; tvals[5] = {}; tmp = x[nx] - x[1]
		for i=1,nts+1 do tvals[5][i] = (i-1)*tmp/nts + x[1] end
	end
	nts,nl = nps*npts,tvals[1]; tmp = (tvals[2][1]-tvals[1])/nts
	nu = nl + tmp
	sol,sta = usave(sol,u),utsave(sta,u,x,nl,tvals[5]) -- save initial values
	for i=1,nts do -- loop over initial interval
		nit = pdebivbv(eqsub,{nl,nu,1,1},x,u,ut,utt)
		if nit>nitt then nitt = nit end
		sta = utsave(sta,u,x,nu,tvals[5]) -- save time values
		nl,nu = nu,nu+tmp
	end
	sol = usave(sol,u) -- Save at end of first interval
	if nprint~=0 then
		printf('Time = %e Number of iterations in pdeivbvqs = %d \n\n' ,nl,nitt)
		io.flush(); nitt = 0
	end
	fact = 10^(1/nts) -- Factor between steps
	nl = 10^(math.floor(math.log10(tvals[2][1])))
	nu = 10^(math.ceil(math.log10(tvals[2][2])))*1.00000000001/fact
	while nl<=nu do -- Step over log time spacings
		nit = pdebivbv(eqsub,{nl,nl*fact,1},x,u,ut,utt)
		if nit>nitt then nitt = nit end
		nl = nl*fact; j = j+1
		sta = utsave(sta,u,x,nl,tvals[5]) -- save time dependent solutions
		if j==npts then 
			sol = usave(sol,u); j = 0 -- Save solutions
			if nprint~=0 then 
				printf('Time = %e Number of iterations in pdeivbvqs = %d \n\n' ,nl,nitt)
				io.flush(); nitt = 0
			end
		end
	end
	return sol, sta
end
setfenv(pdeivbvqst,{getn=table.getn,pdebivbv=pdebivbv,table=table,
printf=printf,type=type,math=math,nprint,io=io,NPS=1,NPTS=10,
usave=usave,utsave=utsave,getfenv=getfenv,pdeivbvqs=pdeivbvqs})

pdeivbvqs1 = function(eqsub,tvals,x,u,ut,utt) -- Quick scan log time solver
	local NPS,NPTS,nps,npts = 1,10
	local ttvals,nl,nu,fact = {}
	local nt,j = #tvals, 1
	local sol
	ut = ut or {}; utt = utt or {}
	if nt<2 then print('Error, must specify two times in pdeivbvqs')
		return end
	nps = math.floor(tvals[3] or NPS)
	npts = math.floor(tvals[4] or NPTS)
	nl = 10^(math.floor(math.log10(tvals[2][1])))
	nu = 10^(math.ceil(math.log10(tvals[2][2])))
	fact = 10^(1/nps) -- Factor between steps
	ttvals[1],ttvals[2],ttvals[3] = tvals[1],{},{}
	while nl<=nu do ttvals[2][j],ttvals[3][j],nl,j = nl,npts,nl*fact,j+1 end
	sol = pdeivbv(eqsub,ttvals,x,u,ut,utt)
	return sol
end

usave = function(u1,u2)
	local nk,neq,nx = #u1+1, #u2, #u2[1]
	for j=1,neq do u1[nk]={}; for k=1,nx do u1[nk][k] = u2[j][k] end; nk = nk+1 end
	return u1
end
utsave = function(u1,u2,x,t,px)
	u1 = u1 or {}
	local nk,neq,nx,nt = #u1, #u2, #x, #px
	if nk~=neq*nt+1 then for i=1,neq*nt+1 do u1[i] = {} end end
	nk = #u1[1]+1
	u1[1][nk] = t
	for i=1,nt do for k=1,neq do 
	u1[neq*(i-1)+1+k][nk] = intp(x,u2[k],px[i]) end end
	return u1
end


	

	

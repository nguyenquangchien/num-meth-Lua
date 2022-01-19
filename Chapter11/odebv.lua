-- /* File odebv.lua */
-- Shooting(ST) method for boundary value problems

require"odeiv"; require"nsolv"
local loc_odebiv = odebiv

odebvst = function(feqs,tvals,u,up) -- Basic shooting boundary value solver
	local s,ns,neq,nu,ni,fbound,err -- If no boundary function, return initial vlaue solution below
	if type(feqs)=='table' then feqs,fbound = feqs[1],feqs[2] 
	else return odeiv(feqs,tvals,u,up) end
	local inbvs,uL,upL = {},{},{}
	nu = #u
	-- Define local function for Newton's method
	local fbeval = function(beqs,bvs) -- bvs[] is array of boundary values
		for i=1,nu do u[i] = bvs[i]; uL[i] = bvs[i] end -- Set initial values
		for i=nu+1,neq do up[i-nu] = bvs[i]; upL[i-nu] = bvs[i] end -- Set initial derivatives
		odeiv(feqs,tvals,u,up) -- Solve initial value problem, only returned u and up are important
		fbound(beqs,uL,u,upL,up) -- Evaluate errors in boundary equations, return in beqs
	end
	-- End of function for Newton's method
	fbound(inbvs,u,u,up,up) -- probe fbound to get number of boundary equations
	neq = #inbvs -- Number of boundary value equations
	if neq~=nu then odebiv = odeb12 end -- Second degree equations present 
	for i=1,nu do inbvs[i] = u[i] end -- Initial values
	for i=nu+1,neq do inbvs[i] = up[i-nu] end -- Initial derivative values
	ni,err = nsolv(fbeval,inbvs) -- Solve boundary value problem subject to fbeval() equations
	for i=1,nu do u[i] = inbvs[i] end -- Final starting values
	for i=nu+1,neq do up[i-nu] = inbvs[i] end -- Final starting derivatives
	s,ns = odeiv(feqs,tvals,u,up) -- Now evaluate final solution for shooting method
	odebiv = loc_odebiv -- Always leave with odebiv() reset
	return s,ni,ns,err
end	

odebvste = function(feqs,tvals,u,up) -- ST Multi-Step BV solver with error estimation
	local nu,ul,ulp = #u,{}
	local ttvals = {tvals[1],tvals[2],{}}
	for i=1,nu do ul[i] = u[i] end
	if up~=nil then ulp = {}; for i=1,nu do ulp[i] = up[i] end end
	local s1,n1 = odebvst(feqs,tvals,u,up)
	if type(tvals[3])=='number' then ttvals[3][1] = tvals[3]/2
	else for i=1,#tvals[3] do ttvals[3][i] = tvals[3][i]/2 end end
	local s2,n2 = odebvst(feqs,ttvals,ul,ulp)
	return s1,odeerror(s1,s2),math.max(n1,n2)
end	

-- Shooting method for boundary value plus eigenvalue problems
odebvev = function(feqs,tvals,ui,upi) -- Basic shooting eigen value solver
	odebiv = odeb12 -- Second order equations
	local nr,s,ns,ni,fbound,E,feq,feqse,err = 'false'
	feq,E,fbound = feqs[1],feqs[2],feqs[3]
	if type(E)=='function' then nr,E,fbound = 'true',fbound,E end
	local u,up,nu = {},{},#ui
	-- Define local default function for boundary values
	local defb = function(beqs,uF) -- Default function for fbound if not input
		for i=1,ns do beqs[i] = uF[i] end -- Zero boundary conditions
	end
	fbound = fbound or defb; ns = #E
	-- Define local function for Newton's method
	local fbeval = function(beqs,E) -- E[] is array of eigenvalues
		feqse = function(eqs,x,du,dup,dupp) -- Local calling function
			feq(eqs,E,x,du,dup,dupp) -- Adds E to calling argument
		end -- Now use in calling equation solver
		for i=1,nu do u[i],up[i] = ui[i],upi[i] end -- Set initial values
		odeiv(feqse,tvals,u,up) -- Solve initial value problem
		fbound(beqs,u,up) -- Update boundary errors
	end -- End of function for Newton's method
	ni,err = nsolv(fbeval,E) -- Solve boundary value problem subject to fbeval() equations
	for i=1,nu do u[i],up[i] = ui[i],upi[i] end
	s,ns = odeiv(feqse,tvals,u,up) -- Update with final solution values
	odebiv = loc_odebiv -- Always leave with odebiv() reset
	if nr then nr,feqs[2],feqs[3] = 'false',fbound,E end -- Reverse
	return s,ns,ni,err,E
end	

odebveve = function(feqs,tvals,u,up) -- ST Multi-Step BV solver with error estimation
	local E1,nu,ul,ulp = {},#u,{}
	local ttvals = {tvals[1],tvals[2],{}}
	for i=1,nu do ul[i] = u[i] end
	if up~=nil then ulp = {}; for i=1,nu do ulp[i] = up[i] end end
	local s1,n1,n2,err,E = odebvev(feqs,tvals,u,up)
	for i=1,#E do E1[i] = E[i] end -- Save eigenvalues
	if type(tvals[3])=='number' then ttvals[3] = tvals[3]/2
	else for i=1,#tvals[3] do ttvals[3][i] = tvals[3][i]/2 end end
	local s2,_,_,_,E = odebvev(feqs,ttvals,ul,ulp)
	for i=1,#E do E[i],E1[i] = E1[i],(4*E1[i]-E[i])/3 end
	return s1,odeerror(s1,s2),E1,math.max(n1,n2)
end	

sqnorm = function(s,area) -- Normalize square of function to given area
	local sn,nx,atot = {}, #s[1], 0
	local x,y = s[1],s[2]
	for i=1,nx do sn[i] = y[i]^2 end -- Square wave function
	for i=1,nx-1 do atot = atot + (x[i+1]-x[i])*(sn[i+1]+sn[i])/2 end
	atot = math.sqrt(area/math.abs(atot))
	for i=1,nx do sn[i] = y[i]*atot end
	return {x,sn},sprob({x, sn}),atot
end

sprob = function(s) -- 
	local sq,n,y = {},#s[1],s[2]
	for i=1,n do sq[i] = y[i]^2 end
	return {s[1],sq}
end

nzeros = function(s,nex)
	local sv,nz,nx,nxx = s[2],0,#s[1]
	nex = nex or math.max(nx/10,10)
	nxx = math.floor(nex)
	for i=nxx,nx-nxx,2 do
		if sv[i]*sv[i+2] < 0 then nz = nz+1 end
	end
	return nz
end
		
	
	
			




-- /* File spline.lua */  Spline function for interpolation
spline = function(xp,yp,x,deriv,d2x)
	if type(deriv)=='table' then deriv,d2x = d2x,deriv end
	deriv = deriv or 0
	local n,nd2x = #xp, true
	local a,b,c,d = {},{},{},{}
	if d2x==nil then d2x = {} -- Table of coefficients input?
	else if #d2x==n then nd2x = false end end
	if n<3 then
		print('Insufficient number of data points in spline function')
		return nil
	end
	-- Set up equation coefficients, a[i]dx2[i-1]+b[i]dx2[i]+c[i]dx2[i+1] = d[i]
	if nd2x then -- Need coeficients?
		a[1],b[1],d[1] = 0,1,0
		b[n],c[n],d[n] = 1,0,0
		if stype==2 then c[1],a[n] = -1, -1 
		else c[1],a[n] = 0, 0 end
		for i=2,n-1 do
			a[i],b[i],c[i] = xp[i]-xp[i-1], 2*(xp[i+1]-xp[i-1]), xp[i+1]-xp[i]
			d[i] = 6/(xp[i+1]-xp[i])*(yp[i+1]-yp[i]) + 6/(xp[i]-xp[i-1])*(yp[i-1]-yp[i])
		end
		-- Solve tridiagonal system of equations
		for i=2,n do -- Forward reduction of off diagonal elements
			ei = a[i]/b[i-1]
			b[i],d[i] = b[i] - ei*c[i-1], d[i] - ei*d[i-1]
		end
		d2x[n] = d[n]/b[n] -- Last coefficient now calculated
		for i=n-1,1,-1 do -- Back substitution for final solution
			d2x[i] = (d[i]-c[i]*d2x[i+1])/b[i] 
		end
	end -- End of coefficients, if needed
	while x>xp[il] and il<n do il = il+1 end
	while x<xp[il-1] and il>2 do il = il-1 end
	if il~=ilold then
		c0 = xp[il]-xp[il-1]
		c1,c2 = d2x[il-1]/(6*c0), d2x[il]/(6*c0)
		c3,c4 = yp[il-1]/c0-d2x[il-1]*c0/6, yp[il]/c0-d2x[il]*c0/6
		ilold = il
	end
	if deriv==0 then 
		return c1*(xp[il]-x)^3+c2*(x-xp[il-1])^3+c3*(xp[il]-x)+c4*(x-xp[il-1]) end
	if deriv==1 then return 3*(-c1*(xp[il]-x)^2+c2*(x-xp[il-1])^2)-c3+c4 end
	if deriv==2 then return 6*(c1*(xp[il]-x)+c2*(x-xp[il-1])) end
end
setfenv(spline, {type=type,print=print,il=2,ilold=0,stype=1})

splinef = function(xd,yd,nd) -- Embed data and convert to f(x)
	local xdd,ydd,d2x = {},{},{} -- Store data locally
	nd = nd or #xd
	if xd[nd]<xd[1] then for i=1,nd do
		xdd[i],ydd[i] = xd[nd+1-i],yd[nd+1-i] end
	else for i=1,nd do 
		xdd[i],ydd[i] = xd[i],yd[i] end end
	return (function(x,drv) return spline(xdd,ydd,x,drv,d2x) end)
end
		

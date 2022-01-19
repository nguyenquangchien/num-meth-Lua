--/*  File intp.lua */

intp = function(xp,yp,x,deriv,ypd) -- LCB -- Local cubic interpolation
	if type(deriv)=='table' then deriv,ypd = ypd, deriv end
	deriv = deriv or 0
	local n = #xp
	local x0,x1,x2,x3,y0,y1,y2,y3,c1,c2,c3,yp1,yp2
	if il>n-2 then il=n-2 end
	if n<3 then fprintf(stderr,
		'Number of points must be > 2 in intp.\n'); return nil
	end
	if xp[n]<xp[1] then 
		fprintf(stderr,"Array of x values in intp must be monotonic increasing\n")
		fprintf(stderr,"Change x values to negative values and try again\n")
	end
	-- Find appropriate interval for requested interpolation
	while x>xp[il] and il<n-1 do il = il+1 end
	while x<xp[il] and il>1 do il = il-1 end
	if ypd and ypd[il] then yp1 = ypd[il] 
	else -- must evaluate derivative at lower point
		if il==1 then 
			x0,x1,x2 = xp[1],xp[2],xp[3]
			y0,y1,y2 = yp[1],yp[2],yp[3]
			c1,c2,c3 = x1-x0,x2-x1,x2-x0
			yp1 = -y0*(c1+c3)/(c3*c1)+y1*c3/(c1*c2)-y2*c1/(c2*c3)
			x1,x2,y1,y2 = x0,x1,y0,y1
		else -- general point
			x0,x1,x2 = xp[il-1],xp[il],xp[il+1]
			y0,y1,y2 = yp[il-1],yp[il],yp[il+1]
			c1,c2,c3 = x1-x0,x2-x1,x2-x0
			yp1 = -y0*c2/(c1*c3)+y1*(c2-c1)/(c1*c2)+y2*c1/(c2*c3)
		end
		if ypd then ypd[il] = yp1 end -- Save if table available
	end
	if ydp and ypd[il+1] then yp2 = ypd[il+1]
	else -- must evaluate derivative at upper point
		if il==n-1 then 
			x1,x2,x3 = xp[n-2],xp[n-1],xp[n]
			y1,y2,y3 = yp[n-2],yp[n-1],yp[n]
			c1,c2,c3 = x3-x2,x2-x1,x3-x1
			yp2 = y1*c1/(c2*c3)-y2*c3/(c1*c2)+y3*(c1+c3)/(c1*c3)
			x1,x2,y1,y2 = x2,x3,y2,y3
		else -- general point
			x1,x2,x3 = xp[il],xp[il+1],xp[il+2]
			y1,y2,y3 = yp[il],yp[il+1],yp[il+2]
			c1,c2,c3 = x3-x2,x2-x1,x3-x1
			yp2 = -y1*c1/(c2*c3)+y2*(c1-c2)/(c1*c2)+y3*c2/(c1*c3)
		end
		if ypd then ypd[il+1] = yp2 end
	end
	c1,c2 = x2-x1,x-x1
	b = (3*(y2-y1)/c1 - 2*yp1 - yp2)/c1
	c = (yp2 + yp1 - 2*(y2-y1)/c1)/(c1*c1)
	if deriv==0 then return y1+((c*c2+b)*c2+yp1)*c2 end
	if deriv==1 then return yp1+(3*c*c2+2*b)*c2 end
	if deriv==2 then return 2*b+6*c*c2 end
end
setfenv(intp,{il=1,fprintf=fprintf,type=type})

intpf = function(xd,yd,nd) -- Embed data and convert to f(x)
	local xdd,ydd,ydrv = {},{},{} -- Store data locally
	nd = nd or #xd
	if xd[nd]<xd[1] then for i=1,nd do -- reverse data
		xdd[i],ydd[i] = xd[nd+1-i],yd[nd+1-i] end
	else for i=1,nd do 
		xdd[i],ydd[i] = xd[i],yd[i] end end
	return (function(x,drv) return intp(xdd,ydd,x,drv,ydrv) end)
end

intpa = function(xyd,np) -- Interpolate on arrays of data
	local na,ndd = #xyd, #xyd[1]
	local xmin,xmax = xyd[1][1],xyd[1][ndd]
	if type(np)=='table' then xp,np = {np},#np
	else
		xp,np = {{}},np+1
		for i=1,np do xp[1][i] = xmin + (i-1)*(xmax-xmin)/(np-1) end
	end
	for j=2,na do
		xp[j] = {}
		for i=1,np do xp[j][i] = intp(xyd[1],xyd[j],ndd,xp[1][i]) end
	end
	return xp
end

intp4 = function(xp,yp,x,deriv) -- Four point cubic interpolation
	deriv = deriv or 0
	local n = #xp
	local x0,x1,x2,x3,y0,y1,y2,y3,c10,c20,c30,c31,c32,c21,d0,d1,d2,d3
	if il>n-2 then il=n-2 end
	if n<3 then fprintf(stderr,
		'Number of points must be > 2 in intp.\n'); return nil
	end
	if xp[n]<xp[1] then 
		fprintf(stderr,"Array of x values in intp must be monotonic increasing\n")
		fprintf(stderr,"Change x values to negative values and try again\n")
	end
	-- Find appropriate interval for requested interpolation
	while x>xp[il] and il<n-1 do il = il+1 end
	while x<xp[il] and il>1 do il = il-1 end
	if il==1 then 
		x0,x1,x2,x3 = xp[1],xp[2],xp[3],xp[4]
		y0,y1,y2,y3 = yp[1],yp[2],yp[3],yp[4]
	elseif il==n-1 then 
		x0,x1,x2,x3 = xp[n-3],xp[n-2],xp[n-1],xp[n]
		y0,y1,y2,y3 = yp[n-3],yp[n-2],yp[n-1],yp[n]
	else
		x0,x1,x2,x3 = xp[il-1],xp[il],xp[il+1],xp[il+2]
		y0,y1,y2,y3 = yp[il-1],yp[il],yp[il+1],yp[il+2]
	end
--print('x,il =',x,il)	
	c10,c20,c30,c31,c32,c21 = x1-x0,x2-x0,x3-x0,x3-x1,x3-x2,x2-x1
	d0,d1,d2,d3 = c10*c20*c30,c10*c21*c31,c20*c21*c32,c30*c31*c32
	x0,x1,x2,x3 = x-x0,x-x1,x-x2,x-x3 -- Now shift to displacements
	if deriv==0 then return -y0*x1*x2*x3/d0+y1*x0*x2*x3/d1-y2*x0*x1*x3/d2+y3*x0*x1*x2/d3 end
	if deriv==1 then return -y0*(x2*x3+x1*x3+x1*x2)/d0+y1*(x2*x3+x0*x3+x0*x2)/d1-
		y2*(x0*x1+x0*x3+x1*x3)/d2+y3*(x0*x2+x0*x1+x1*x2)/d3 end
	if deriv==2 then return 2*(-y0*(x1+x2+x3)/d0+y1*(x0+x2+x3)/d1-y2*(x0+x1+x3)/d2+
		y3*(x0+x1+x2)/d3) end
end
setfenv(intp4,{il=1,fprintf=fprintf,print=print})

intp4f = function(xd,yd) -- Embed data and convert to f(x)
	local xdd,ydd = {},{} -- Store data locally
	nd = #xd
	if xd[nd]<xd[1] then for i=1,nd do -- reverse data
		xdd[i],ydd[i] = xd[nd+1-i],yd[nd+1-i] end
	else for i=1,nd do 
		xdd[i],ydd[i] = xd[i],yd[i] end end
	return (function(x,drv) return intp4(xdd,ydd,x,drv) end)
end

lintp = function(x,y,xu,deriv) -- Linear interpolation with first derivative
	n = #x
	while xu>x[il+1] and il<n-1 do il = il+1 end
	while xu<x[il] and il>1 do il = il-1 end
	y1,y2,x1,x2 = y[il],y[il+1],x[il],x[il+1]
	if deriv==1 then return (y2-y1)/(x2-x1) end
	return y1 + (y2-y1)*(xu-x1)/(x2-x1)
end		
setfenv(lintp,{il=1})	

lintpf = function(xd,yd)
	local xdd,ydd = {},{} -- So user can't change data
	local nd = #xd
	if xd[nd]<xd[1] then for i=1,nd do
		xdd[i],ydd[i] = xd[nd+1-i],yd[nd+1-i] end
	else for i=1,nd do 
		xdd[i],ydd[i] = xd[i],yd[i] end end
	return function(x,drv) return lintp(xdd,ydd,x,drv) end
end


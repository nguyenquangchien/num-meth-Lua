-- /* File Data_Fit.lua */
-- Program to fit data with noise to a smooth curve for plotting
-- Number of knots specified by user or determined by program <= 10 floating knots
require"nlstsq"; require"intp"
local checkorder

datafit = function(yx,nc)
	local nd,c,xc = #yx[1] -- yx array contains {yd, xd}
	yx = checkorder(yx,nd) -- increasing x array values?
	local xmin,xmax = yx[2][1], yx[2][nd] -- x min and max values 
	local actv,dx,del = {}
	local fdf = function(yx,c)  -- Interpolation function for datafit()
		return intp(xc,c,yx[2])-yx[1] -- Interpolation function with set of knots
	end
	if type(nc)=='table' then -- arrays for xc and fixed knots?
		xc,c = nc[1] or {}, nc[2] or {} -- floating knots and c?
		nc = #xc -- All xc values must be passed
	else 
		if (nc==nil) or (nc==0) then -- Use default number of knots 
			nc = math.min(nd/10,15); nc = math.max(nc,3)
		end -- End result should be between 3 and 15
		nc,xc,c = math.ceil(nc), {}, {}
	end
	dx = (xmax-xmin)/(nc-1) -- Use equal spacing on x-axis for knots
	for i=1,nc do -- Set initial x and y values of floating knots
		 xc[i] = xc[i] or (xmin + (i-1)*dx) -- Take passed values or calculate 
		 if c[i] ~=nil then actv[i] = 0 else actv[i] = 1 end -- value or float?
		c[i] = c[i] or 0 -- linear problem, initial guesses can be zero
	end
	local nlenv = getfenv(nlstsq); local lcp,lyl = nlenv.nprint, nlenv.ylinear
	getfenv(nlstsq).nprint = 0; getfenv(nlstsq).ylinear = 1 -- Linear in y , no print
	del,err = nlstsq(yx, fw, fdf, c, actv) -- call data fit
	nlenv.nprint, nlenv.ylinear = lcp, lyl -- restore values in nlstsq
	return function(x,df) return intp(xc,c,x,df) end, {xc, c, del}, err -- f(x)
end

refitdata = function(yx,intpf,xcc) -- Add two additional knots for each call
	local j,nd = 0,#yx[1]
	local xmin,xmax = 0,0
	local c,nc = xcc[2],#xcc[2]
	local xc,del -- Needed for interpolation function
	local fdf = function(yx,c)  -- Interpolation function for datafit()
		return intp(xc,c,yx[2])-yx[1] -- Interpolation function with set of knots
	end
	xc = xcc[1] -- Global value for fdf
	for i=2,nc-1 do -- Find point of maximum second derivative change
		xmin = math.abs(intpf(xc[i]*(1.001),2) - intpf(xc[i]*(.9999),2))
		if xmin>xmax then j = i; xmax = xmin end
	end
	table.insert(xc,j+1,(xc[j+1] + xc[j])*.5) -- Add first point
	table.insert(xc,j,(xc[j] + xc[j-1])*.5)  -- Two additional points now in xc 
	table.insert(c,j+1,0); table.insert(c,j,0)
	nc = nc + 2
	del,err = nlstsq(yx,fw,fdf,c) -- call data fit with 2 extra knots
	return function(x,df) return intp(xc,c,x,df) end, {xc, c, del}, err -- f(x)
end
datafitn = function(yx,na,nc)
	local intpf,xcc,err = datafit(yx,nc)
	for i=1,na do
		intpf,xcc,err = refitdata(yx,intpf,xcc)
	end
	return intpf,xcc,err
end
checkorder = function(yx,nd) -- Check order of x array
	local ord,yx2 = 1, yx[2]
	for i=1,nd-1 do if yx2[i+1]<yx2[i] then ord = 0; break end end
	if ord==1 then return yx end -- Already ordered array
	local xx,yy = {},{} -- Set up ordered x array with corresponding y array
	for i=1,nd do yy[i],xx[i] = yx[1][i],yx2[i] end
	local xmin,imin
	for i=1,nd do -- Swap values to order x with corresponding y values
		xmin,imin = xx[i],i
		for j=i,nd do
			if xx[j]<xmin then xmin,imin = xx[j],j end
		end
		yy[i],yy[imin],xx[i],xx[imin] = yy[imin],yy[i],xx[imin],xx[i]
	end
	return {yy,xx}
end


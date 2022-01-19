--/* File intg.lua */ -- Integration with Romberg enhancement
intg = function(xmin,xmax,f,...)
	local nend,n = RMAX, NINIT -- Initially set to 100 slices or intervals
	local dx,iend,p2,nend = (xmax - xmin)/n, 0.0 -- Special end treatment if infinite
	fxx = tostring(f(xmin,...)) -- if statements handle infinite function values at boundaries.
	if fxx==inf or fxx==und then -- Infinite lower limit value of function, inf = 1/0
		dx = dx/NDX
		iend = fend(f(xmin+dx,...),f(xmin+dx/2,...),dx)
		xmin = xmin+dx; dx = (xmax-xmin)/n
	end
	fxx = tostring(f(xmax,...))
	if fxx==inf or fxx==und then -- Infinite upper limit value of function
		dx = dx/NDX
		iend = iend + fend(f(xmax-dx,...),f(xmax-dx/2,...),dx)
		xmax = xmax-dx; dx = (xmax-xmin)/n
	end
	local a1 = 0.5*(f(xmin,...)+f(xmax,...)) -- End point values for trapezoidal rule
	for i=1,n-1 do a1 = a1 + f(xmin+i*dx,...) end -- Sum internal points
	local a = {dx*a1} -- Initial area estimate from trapezoidal rule
	for j=1,RMAX do -- Romberg cycles -- 12 maximum
		n,dx = 2*n, dx/2 -- Double number of slices and half spacing of slices
		for i=1,n-1,2 do a1 = a1 + f(xmin+i*dx,...) end -- Sum over additional slices
		a[j+1] = dx*a1 -- Updated area with additional slices
		p2 = 2; alast=a[1]
		for k=j,1,-1 do -- Romberg convegence accelerations
			a[k] = (p2*a[k+1] - a[k])/(p2-1)
			p2 = 2*p2
		end
		p2 = abs(a[1]-alast)/abs(a[1]+alast) -- Error estimate for this cycle
		if p2<ERROR then nend=j; break end
	end
	return a[1]+iend,nend,abs(a[1]-alast) -- Add area of ends and return with number of cycles
end

local fend = function(f1,f2,dx) -- Approximate values of ends for infinite f(x)
	local n = log(f2/f1)/log(2) -- Approximates with power law on x near ends, f1*x^(-n)
	if n>=1 then print("Value of integral appears to be infinite") end
	return f1*dx/(1-n)
end	
setfenv(intg,{abs=math.abs,RMAX=15,NINIT=100,NDX=1000,ERROR=1.e-11,tostring=tostring,inf=tostring(math.huge),und=tostring(0/0),fend=fend})
setfenv(fend,{log=math.log,print=print})
	
intg_inf = function(f,...) -- evaluate integral from 0 to infinity, using folded interval
	return intg(0,1,(function(x,...) return f(x,...)+f(1/x,...)/x^2 end),...)
end


--/* File deriv.lua */ -- Derivative with Richardson's extrapolation
deriv = function(f, x,...)  -- Central derivative of f(x)
	local dx = DXX*x -- DXX = .01 to begin
	if dx==0.0 then dx = DXX end
	local a = {(f(x+dx,...) - f(x-dx,...))/(2*dx)}
	local p2, nj
	for j=1,NMAX do -- Richardson's improvement
		dx = dx/2; p2 = 4 -- Second order error
		a[j+1] = (f(x+dx,...) - f(x-dx,...))/(2*dx)
		for k=j,1,-1 do
			a[k],p2 = (p2*a[k+1] - a[k])/(p2-1), 4*p2
		end
		nj = j 
		p2 = abs(a[1]-a[2])/abs(a[1]+a[2])
		if p2<err then break end
	end
	return a[1],nj,p2
end
setfenv(deriv,{abs=(Complex or math).abs,NMAX=10,DXX=1.e-2,err=1.e-14})
deriv2 = function(f, x,...) -- Function for second derivative of f(x)
	return deriv((function(x,...) return deriv(f, x,...) end), x,...)
end
rderiv = function(f, x,...) -- Right derivative of f(x)
	local dx = abs(DXX*x) -- DXX = .01 to begin
	if dx==0.0 then dx = DXX end
	local x0 = x + dx*DXD -- A little right of point
	local a = {(f(x0+dx,...)-f(x0,...))/dx} -- One sided derivative
	local p2, nj
	for j=1,NMAX do -- Richardson's improvement
		dx = dx/2; p2 = 2 -- First order error
		x0 = x + dx*DXD
		a[j+1] = (f(x0+dx,...)-f(x0,...))/dx
		for k=j,1,-1 do
			a[k],p2 = (p2*a[k+1] - a[k])/(p2-1), 2*p2
		end
		nj = j 
		p2 = abs(a[1]-a[2])/abs(a[1]+a[2])
		if p2<err then break end
	end
	return a[1], nj, p2
end
setfenv(rderiv,{abs=(Complex or math).abs,NMAX=10,DXX=1.e-2,DXD=1.e-6,err=1.e-14})
rderiv2 = function(f, x,...) -- Right second derivative of f(x)
	return rderiv((function(x,...) return rderiv(f, x,...) end), x,...)
end
lderiv = function(f, x,...) -- Left derivative of f(x)
	local dx = abs(DXX*x) -- DXX = .01 to begin
	if dx==0.0 then dx = DXX end
	local x0 = x - DXD*dx -- A little left of point
	local a = {(f(x0,...) - f(x0-dx,...))/dx}
	local p2, nj
	for j=1,NMAX do -- Richardson's improvement
		dx = dx/2; p2 = 2 -- First order error
		x0 = x - dx*DXD
		a[j+1] = (f(x0,...) - f(x0-dx,...))/dx
		for k=j,1,-1 do
			a[k],p2 = (p2*a[k+1] - a[k])/(p2-1), 2*p2
		end
		nj = j 
		p2 = abs(a[1]-a[2])/abs(a[1]+a[2])
		if p2<err then break end
	end
	return a[1], nj, p2
end
setfenv(lderiv,getfenv(rderiv))
lderiv2 = function(f, x,...) -- Left second derivative of f(x)
	return lderiv((function(x,...) return lderiv(f, x,...) end), x,...)
end
sderiv = function(f,x,...)
	local dx = abs(DXX*x)
	if dx==0.0 then dx = DXX end
	return (f(x+dx,...)-f(x-dx,...))/(2*dx)
end
setfenv(sderiv,{DXX=1.e-6,abs=(Complex or math).abs})
function pderiv(f,x,i)
	local xi,fv = x[i]
	local dx = xi*FACT
	if dx==0 then dx = FACT end
	x[i] = xi + dx; fv = f(x)
	x[i] = xi - dx; fv = (fv - f(x))/(2*dx)
	x[i] = xi; return fv
end
setfenv(pderiv,{FACT=1.e-6})
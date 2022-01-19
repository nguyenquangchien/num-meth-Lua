-- /* File pderiv.lua */ -- Partial derivative 
function pderiv(f,x,i)
	local xi,fv = x[i]
	local dx = xi*FACT
	if dx==0 then dx = FACT end
	x[i] = xi + dx; fv = f(x)
	x[i] = xi - dx; fv = (fv - f(x))/(2*dx)
	x[i] = xi; return fv
end
setfenv(pderiv,{FACT=1.e-6})


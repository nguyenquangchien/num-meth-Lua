-- /* File list5_7.lua */ -- Constrained optimization

require"deriv"; require"nsolv"

xio = {1,2,-3,4,-5} -- minimum points
f = function(x) -- Multi-variable function
	local fv = 0
	for i=1,nx do fv = fv + (x[i] - xio[i])^2 end
	return math.exp(fv)
end
gs = { -- Table of constraint functions
	function(x) return x[1]-x[2]-1 end, -- #1
	function(x) return x[3]+x[4] end, -- #2
	function(x) return x[4]+x[5] end -- #3
}

nx,nl = #xio, #gs -- number variables and constraints
x = {0,0,0,0,0,0,0,0} -- Initial guess at minimum
dx = -.2; step = {dx,dx,dx,dx,dx,0,0,0} -- Limit steps

eqs = function(y,x) -- Equations to be solved
	local yv
	for i=1,nx do 
		yv = pderiv(f,x,i) 
		for j=1,nl do yv = yv + x[nx+j]*pderiv(gs[j],x,i) end
		y[i] = yv
	end
	for i=1,nl do y[nx+i] = gs[i](x) end
end

print(nsolv(eqs,x,step)) -- Now solve them
print('Solution values are:')
table.foreach(x,print) -- Print minimum point
print('Value at solution =',f(x))
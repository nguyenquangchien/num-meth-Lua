-- File list5_6.lua -- Minimum problem using pderiv()

require"deriv"; require"nsolv"

xio = {1,2,-3,4,-5} -- minimum points
fuv = function(x) -- Multi-variable function
	local fv = 0
	for i=1,#x do fv = fv + (x[i] - xio[i])^2 end
	return math.exp(fv)
end

x = {0,0,0,0,0} -- Initial guess at minimum
print('At initial guess')
print('Partial 1 derivative =',pderiv(fuv,x,1)) -- Just check on partial derivatives
print('Partial 3 derivative =',pderiv(fuv,x,3))

eqs = function(y,x) -- Equatuions to be solved
	for i=1,#x do 
		y[i] = pderiv(fuv,x,i) -- Force partial derivatives to zero
	end
end

print(nsolv(eqs,x)) -- Now solve them
print('Minimum occurs at')
table.foreach(x,print) -- Print minimum point
print('Minimum value =',fuv(x))
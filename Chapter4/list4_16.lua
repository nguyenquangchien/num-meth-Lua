-- /* File list4_16.lua */ -- Simple examples of polynomial operations

Polynomial = {}; Polynomial_mt = {} -- Polynomial table and metatable

Polynomial.new = function(...) -- Creates a new Polynomial type with associated metatable
	return setmetatable({...},Polynomial_mt) 
end
Polynomial_mt.__add = function(p1,p2)
	local n,sum = math.max(#p1, #p2), {}
	for i=1,n do
		sum[i] = (p1[i] or 0) + (p2[i] or 0)
	end
	return setmetatable(sum,Polynomial_mt)
end
Polynomial_mt.__unm = function(p)
	local pp = {}
	for i=1,#p do pp[i] = -p[i] end
	return setmetatable(pp,Polynomial_mt)
end
Polynomial_mt.__sub = function(p1,p2)
	p2 = -p2; return p1 + p2
end
Polynomial_mt.__tostring = function(p) -- Returns a string in standard form a + bx +cx^2 + --
	s = tostring(p[1])
	for i=2,#p do
		s = s.." + ("..tostring(p[i])..")*x^"..tostring(i-1)
	end
	return s
end
Polynomial_mt.__mul = function(p1,p2) -- Multiplies two polynomials p1 and p2
	local n1,n2,pp,k,fact = #p1, #p2, {}
	for i=1,n1 do
		k,fact = i, p1[i]
		for j=1,n2 do
			pp[k],k = (pp[k] or 0) + fact*p2[j], k+1
		end
	end
	return setmetatable(pp,Polynomial_mt)
end	

-- Now some examples of Polynomial operations
p1 = Polynomial.new(1,2,3,4)
p2 = Polynomial.new(4,3,2,1)
print(p1); print(p2)
p3 = p1*p2; print(p3)
print(p1-p3)

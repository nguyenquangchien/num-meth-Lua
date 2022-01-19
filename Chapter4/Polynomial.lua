-- /* File Polynomial.lua */
-- Polynomial Class and associated functions

require"Complex"; require"nsolv"; require"newton"

Polynomial = {} -- Define class known as Polynomial, with null metatable
setmetatable(Polynomial,{})
-- Define metatable operations 
local Polynomial_mt = {__index = Polynomial, __type = 'Polynomial'} -- Polynomial
local PolynomialR_mt = {__index = Polynomial, __type = 'PolyR'} -- Ratio of polynomials

Polynomial.new = function(p,...) -- Creates a new Polynomial type with associated metatable
	local p2,t1,t2 = select(1,...),type(p)
	if p2==nil then p2,t2 = 1,nil else t2 = type(p2) end
	if t1=='number' then return setmetatable({p,...},Polynomial_mt) end
	if p2==1 then 
		return setmetatable({unpack(p)},getmetatable(p) or Polynomial_mt) 
	end
	if t2=='number' then p2 = setmetatable({...},Polynomial_mt) end
	if t2=='table' then p2 = setmetatable(p2,Polynomial_mt) end
	t2 = type(p2); 
	if t1=='table' then 
		t1 = 'Polynomial'; setmetatable(p,Polynomial_mt) 
	end
	if t1=='Polynomial' and t2=='Polynomial' then -- PolyR type
		return setmetatable({p,p2},PolynomialR_mt)
	end
	return p/p2 -- Perform polynomial math and return new polynomial
end

Polynomial_mt.__add = function(p1,p2)
	local t1,t2 = type(p1), type(p2)
	if t1=='PolyR' or t2=='PolyR' then 
		if t1~='PolyR' then p1,p2 = p2,p1 end
		if type(p2)~='PolyR' then return Polynomial.new(p1[1]+p2*p1[2],p1[2])
		else return Polynomial.new(p1[1]*p2[2]+p2[1]*p1[2],p1[2]*p2[2]) end
	end
	if t1=='number' then p1,p2 = p2,p1 end
	if t2=='number' then p2 = {p2} end
	local n,sum = math.max(#p1, #p2), {}
	for i=1,n do
		sum[i] = (p1[i] or 0) + (p2[i] or 0)
	end
	return setmetatable(sum,Polynomial_mt)
end
PolynomialR_mt.__add=Polynomial_mt.__add

Polynomial_mt.__unm = function(p)
	local pp = {}
	for i=1,#p do
		pp[i] = -p[i]
	end
	return setmetatable(pp,Polynomial_mt)
end
PolynomialR_mt.__unm = function(p)
	return setmetatable({-p[1],p[2]},PolynomialR_mt) 
end

Polynomial_mt.__sub = function(p1,p2)
	p2 = -p2
	return p1+ p2
end
PolynomialR_mt.__sub = Polynomial_mt.__sub

Polynomial_mt.__tostring = function(p) -- Returns a string in standard form a + bx +cx^2 + --
	s = tostring(p[1])
	for i=2,#p do
		if p[i]~=0.0 then s = s.." + ("..tostring(p[i])..")*x^"..tostring(i-1) end
	end
	return s
end
PolynomialR_mt.__tostring = function(p)
	return '{'..tostring(p[1])..'}/{'..tostring(p[2])..'}'
end
Polynomial_mt.__mul = function(p1,p2) -- Multiplies two polynomials p1 and p2
	local t1,t2 = type(p1), type(p2)
	--if t2=='PolyR' then return p2*p1 end
	if t1=='PolyR' or t2=='PolyR' then 
		if t1~='PolyR' then p1,p2 = p2,p1 end
		if type(p2)~='PolyR' then return Polynomial.new(p2*p1[1],p1[2])
		else return Polynomial.new(p1[1]*p2[1],p1[2]*p2[2]) end
	end
	if t1=='number' then p1,p2 = p2,p1 end
	if type(p2)=='number' then p2 = {p2} end
	local n1,n2,pp,k,fact = #p1, #p2, {}
	for i=1,n1 do
		k,fact = i, p1[i]
		for j=1,n2 do
			pp[k],k = (pp[k] or 0) + fact*p2[j], k+1
		end
	end
	return setmetatable(pp,Polynomial_mt)
end	
PolynomialR_mt.__mul = Polynomial_mt.__mul

Polynomial_mt.__div = function(p1,p2) -- Division -- Polynomial ratios
	t1,t2 = type(p1), type(p2)
	if t1=='PolyR' or t2=='PolyR' then 
		if t2~='PolyR' then return Polynomial.new(p1[1],p1[2]*p2)
		elseif type(p1)~='PolyR' then return Polynomial.new(p1*p2[2],p2[1])
		else return Polynomial.new(p1[1]*p2[2],p1[2]*p2[1]) end
	end
	local pp = {}
	if type(p1)=='Polynomial' then
		if type(p2)=='number' then
			for i=1,#p1 do pp[i] = p1[i]/p2 end
			return setmetatable(pp,Polynomial_mt)
		end
	end
	return Polynomial.new(p1,p2) -- Make PolyR type
end
PolynomialR_mt.__div = Polynomial_mt.__div

Polynomial_mt.__pow = function(p,n) -- Power of a polynomial
	local fac = {1}
	for i=1,n do
		fac = fac*p
	end
	return Polynomial.new(fac)
end
PolynomialR_mt.__pow = function(p,n)
	return Polynomial.new(p[1]^n,p[2]^n)
end

Polynomial_mt.__call = function(p,x)
	return Polynomial.value(p,x)
end
PolynomialR_mt.__call = function(p,x)
	return Polynomial.value(p[1])/Polynomial.value(p[2])
end
function Polynomial.scale(p,f) -- Scales a polynomial by a factor f -- x -> x*f ; x^2 -> x^2*f^2
	local fac,pp = 1, {}
	for i=1,#p do
		pp[i],fac = p[i]/fac, fac*f
	end
	return setmetatable(pp,Polynomial_mt)
end
function Polynomial.div(p1,p2) -- Long division 
	n,m = #p1,#p2
	while p1[n]==0.0 do p1[n],n = nil,n-1 end -- Zeros for highest powers?
	while p2[m]==0.0 do p2[m],m = nil,m-1 end
	if m>n then return Polynomial.new{0}, p1 end 
	local a,b = {}, {}; fac =  p2[m]
	for i=1,n do a[i] = p1[i]/fac end; for i=1,m do b[i] = p2[i]/fac end
	for i=n,m,-1 do 
		for j=1,m-1 do a[i-j] = a[i-j] - a[i]*b[m-j] end
	end
	for i=1,m-1 do b[i] = a[i]*fac end; b[m] = nil 
	for j=1,n do if j>n-m+1 then a[j] = nil else a[j] = a[j+m-1] end end
	while b[m-1]==0.0 do b[m-1],m = nil,m-1 end
	return setmetatable(a,Polynomial_mt),setmetatable(b,Polynomial_mt)
end
function Polynomial.roots(p) -- Solves for roots of polynomial p
	if type(p)=='PolyR' then -- Ratios of polynomials
		return {Polynomial.roots(p[1]),Polynomial.roots(p[2])} end
	local NMAX = getfenv(nsolv).NMAX -- Maximum number of iterations in nsolv
	local roots,pr = {}, {} -- Collect roots, local copy of p
	local n, m, f1, pol = #p, 1
	local CMP,REL,typ,nmx,iroot = -1,1,1,1,0
	while p[n]==0 and n>2 do n=n-1 end -- Check for zero's at end
	while p[nmx]==0 and nmx<n-2 do -- Check for zeros at beginning
		iroot = iroot+1
		roots[iroot], nmx = 0, nmx+1 -- Zero roots
	end
	for i=nmx,n do pr[m],m = p[n-m+1],m+1 end -- Revese polynomial
	n = n-nmx+1
	
	local function Pfactor(p1,c1,c2)
		local pr = {} -- leave p1 unchanged
		pr[1],pr[2] = p1[1], p1[2]-c1*p1[1]
		for i=3,n-1 do pr[i]=p1[i]-c1*pr[i-1] - c2*pr[i-2] end
		pr[n] = (p1[n] - c2*pr[n-2])
		return pr -- residuals in n-1 and n terms
	end
	local function eqs(y,c) -- Define function to force to zero
		if typ==REL then pol = Pfactor(pr, c[1]+c[2],c[1]*c[2]) -- Real pair
		else pol = Pfactor(pr, 2*c[1],c[1]^2+c[2]^2) end  -- Complex pair
		y[1],y[2] = pol[n-1], pol[n] -- Force remainder to zero
	end
	
	while n>3 do -- solve for quadtatic factors
		c = {c0[1],c0[2]}; nmx = nsolv(eqs,c)
		if nmx==NMAX then 
			typ,c = -typ, {c0[1],c0[2]}
			nmx = nsolv(eqs,c) -- Try again, Reverse real & complex
			if nmx==NMAX then -- Didn't succeed at finding root pair
				print('Iterations exceed NMAX of ',NMAX, ' in poly_root')
				print('At root number', iroot+1)
				print('Try using "getfenv(nsolv).NMAX = ',2*NMAX,'"')
			end
		end
		if typ==REL then pr = Pfactor(pr, c[1]+c[2],c[1]*c[2]) -- Real pair
		else pr = Pfactor(pr,2*c[1],c[1]^2+c[2]^2) end  -- Complex pair
		n,iroot = n-2, iroot+1
		if typ==CMP then roots[iroot],roots[iroot+1] = -c[1]+j*c[2], -c[1]-j*c[2]
		else roots[iroot],roots[iroot+1] = -c[1], -c[2] end
		iroot = iroot+1; pr[n+2], pr[n+1] = nil, nil -- Reduce order of polynomial by 2
	end
	iroot = iroot + 1 -- Final roots from either quadratic factor or linear factor
	if n>2 then -- Quadratic factor left
		f1 = pr[2]^2 - 4*pr[1]*pr[3]
		if f1>=0 then f1=sqrt(f1) else f1 = j*sqrt(-f1) end
		if pr[2]>0 then f1 = -0.5*(pr[2] + f1) else f1 = -0.5*(pr[2] - f1) end
		roots[iroot] = f1/pr[1]
		iroot = iroot + 1; roots[iroot] = pr[3]/f1
	else roots[iroot] = -pr[2]/pr[1]  end -- Linear factor left
	return roots
end	
setfenv(Polynomial.roots,{sqrt=math.sqrt,nsolv=nsolv,j=Complex.j,getfenv=getfenv,table=table,
	Polynomial=Polynomial,print=print,type=type,c0={1,-1}})
	
function Polynomial.newton(p,x) -- Newton method to find single root of polynomial equation
	local function fx(x) -- function to be called by newton() method
		local sum,fac = 0, 1
		for i=1,#p do -- Calculate value of polynomial
			sum,fac = sum + p[i]*fac, fac*x
		end
		return sum -- Return value of polynomial
	end -- end of fx function
	if type(x)=='table' then -- Called with table of roots
		local rts = {}
		for i=1,#x do rts[i] = newton(fx,x[i]) end
		return rts
	else -- Called with one root value
		return newton(fx,x) -- newton returns the value of x that makes fx=0
	end
end
function Polynomial.deriv(p) -- Returns derivative of polynomial
	local dp = {}
	for i=1,#p-1 do
		dp[i] = p[i+1]*i
	end
	return setmetatable(dp, Polynomial_mt) 
end
function Polynomial.intg(p) -- Returns integral of polynomial, c[1]=0
	local ip = {}
	ip[1] = 0
	for i=1,#p do
		ip[i+1] = p[i]/i
	end
	return setmetatable(ip, Polynomial_mt) 
end

function Polynomial.value(p,x) -- Evaluate the value of a polynomial
	if type(p)=='PolyR' then return p[1]:value(x), p[2]:value(x) end
	if (type(p)=='Polynomial') or (type(p)=='table') then
		local sum, val = 0, 1
		for i=1,#p do sum, val = sum + p[i]*val, val*x end
		return sum
	end
end

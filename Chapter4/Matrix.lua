-- /* File Matrix.lua */
-- Package for Matrix algebra and functions -- Restricted to 2 dimensional matrices

require'spgauss' -- Gauss solution code for system of linear equations
local abs = (Complex or math).abs -- Possible complex numbers
Matrix = {}
setmetatable(Matrix,{})
local Matrix_mt = {__index = Matrix, __type = 'Matrix'}

-- Define metatable operations of +,-,*,/,^ and unary -
Matrix_mt.__add = function(mx,my) -- Add function
	if type(mx)~='Matrix' or type(my)~='Matrix' then
		print('Can only add two matrices of similar size'); return
	end
	local m,n = Matrix.size(mx)
	if m~=#my then 
		print('Unequal row numbers in Matrix.add'); return
	end
	if n~=#my[1] then
		print('Unequal column numbers in Matrix.add'); return
	end
	local ml = Matrix.new(m)
	for i=1,m do for j=1,n do ml[i][j] = mx[i][j] + my[i][j] end end
	return setmetatable(ml,Matrix_mt)
end
			
Matrix_mt.__sub = function(mx,my) -- Subtract function
	if type(mx)~='Matrix' or type(my)~='Matrix' then
		print('Can only subtract two matrices of similar size'); return
	end
	local m,n = Matrix.size(mx)
	if m~=#my then 
		print('Unequal row numbers in Matrix.sub'); return
	end
	if n~=#my[1] then
		print('Unequal column numbers in Matrix.sub'); return
	end
	local ml = Matrix.new(m,n)
	for i=1,m do for j=1,n do ml[i][j] = mx[i][j] - my[i][j] end end
	return setmetatable(ml,Matrix_mt)
end
			
Matrix_mt.__mul = function(mx,my) -- Multiply function
	if type(mx)~='Matrix' then mx,my = my,mx end
	local m1,n1 = Matrix.size(mx)
	local ml = Matrix.new(m1,n1) 
	local sum
	if type(my)~='Matrix' then
		if type(my)=='number' then 
			for i=1,m1 do for j=1,n1 do ml[i][j] = my*mx[i][j] end end
		end
	else
		local m2,n2 = Matrix.size(my)
		if n1~=m2 then 
			print('Inconsistent Matrix sizes in Matrix.mul'); return
		end
		for i=1,m1 do
			for j=1,n2 do
				sum = 0
				for k=1,n1 do sum = sum + mx[i][k]*my[k][j] end
				ml[i][j] = sum
			end
		end
	end
	return setmetatable(ml,Matrix_mt)
end

Matrix_mt.__pow = function(mx,pow) -- Power function 
	local m,n = Matrix.size(mx)
	if m~=n then 
		print('Can take power of only square matrix'); return
	end
	local ml = Matrix.new(m)
	if pow==0 then
		for i=1,m do
			for j=1,n do if i==j then ml[i][j] = 1 else ml[i][j] = 0 end end
		end
	elseif pow<=-1 then 
		ml = Matrix.inv(mx) 
		if type(ml)=='number' then
			print('Singular matrix found when taking inverse',ml)
			ml = Matrix.new(m) -- Try to return something useful - Unity matrix
			for i=1,m do
				for j=1,m do if i==j then ml[i][j] = 1 else ml[i][j] = 0 end end
			end
		end
		if pow<-1 then -- Handle large negative powers
			setmetatable(ml,Matrix_mt)
			return ml^(-pow) -- Let __pow() handle with positive power
		end
	else
		ml = Matrix.new(mx)
		while pow>1 do ml,pow = ml*mx,pow-1 end
	end
	return setmetatable(ml,Matrix_mt)
end

Matrix_mt.__unm = function (mx) -- Negative of a matrix
	local m,n = Matrix.size(mx)
	local ml = Matrix.new(m,n)
	local row,rowj
	for i=1,m do 
		row,rowj = ml[i],mx[i]
		for j=1,n do row[j] = -rowj[j] end
	end
	return setmetatable(ml,Matrix_mt)
end

function Matrix_mt.__tostring(mx) -- Convert to string for printing
	local m,n = Matrix.size(mx)
	local s = '{{'
	for i=1,m do
		for j=1,n-1 do 
		s = s..(tostring(mx[i][j]) or '..')..', ' end
		s = s..(tostring(mx[i][n]) or '..')..'}'
		if i==m then s = s..'}' else s = s..',\n{' end
	end
	return s -- Printable string
end
	
-- Matrix definition
function Matrix.new(a,m,n) -- Define new Matrix
	if type(a)~='Matrix' then m,n = a, m or a end
	local mt = {} -- Make new array for Matrix
	if type(a)=='Matrix' or type(a)=='table' then
		m,n = Matrix.size(a) -- number of rows, columns
		for i=1,m do
			mt[i] = {}
			for j=1,n do mt[i][j] = a[i][j] end
		end
	else
		for i=1,m do mt[i] = {} end -- Null matrix elements
	end
	return setmetatable(mt, Matrix_mt) -- Set Matrix metatable
end

function Matrix.size(m) -- Return size of i,j matrix
	return #m, #m[1]
end

function Matrix.LUdecompose(a,n) -- LU decomposition of a matrix
	local d,TINY = 1,1.e-100
	local imax,big,dum,sum,temp
	local vv,indx = {},{}
	n = n or #a
	for i=1,n do -- Loop over rows for scaling informatio
		big = 0
		for j=1,n do -- Find largest element in row j
			temp = abs(a[i][j])
			if temp > big then big = temp end
		end
		if big==0 then print("Singular martix in LUdecompose") end
		vv[i] = 1/big -- Scale factor
	end
	for j=1,n do -- This is the main loop for Crout's method
		for i=1,j-1 do
			sum = a[i][j]
			for k=1,i-1 do sum = sum - a[i][k]*a[k][j] end
			a[i][j] = sum
		end
		big = 0
		for i=j,n do
			sum = a[i][j]
			for k=1,j-1 do sum = sum - a[i][k]*a[k][j] end
			a[i][j],dum = sum, vv[i]*abs(sum)
			if dum>= big then big,imax = dum,i end
		end
		if j~=imax then -- Interchange rows
			a[imax],a[j] = a[j],a[imax] 
			vv[imax],d = vv[j],-d
		end
		indx[j] = imax
		if a[j][j]==0 then a[j][j] = TINY end -- Singular matrix, but continue
		if j~=n then -- Divide by pivot element
			dum = 1/a[j][j]
			for i=j+1,n do a[i][j] = a[i][j]*dum end
		end
	end
	return a,indx,d -- LU matrix, interchange table and sign
end

function Matrix.LUbsolve(a,indx,bb) -- Solve system of equations by LU
	local n,ii,ip,sum = #bb,0
	local b = {} -- bb can be a column matrix or a row vector
	if type(bb[1])=='table' then for i=1,n do b[i]  = bb[i][1] end
	else b = bb end
	for i=1,n do -- Main loop, Forward substitution
		ip = indx[i];	sum,b[ip] = b[ip],b[i]
		if ii~=0 then for j=ii,i-1 do sum = sum - a[i][j]*b[j] end
		else if sum~=0 then ii = i end end
		b[i] = sum
	end
	for i=n,1,-1 do -- Main loop, Backward substitution
		sum = b[i]
		for j=i+1,n do sum = sum - a[i][j]*b[j] end
		b[i] = sum/a[i][i]
	end
	return b -- Solution set
end

function Matrix.LUfunction(a,n) -- Return function to solve equations
	-- by LU decomposition -- stores LU martix locally
	local a,indx = Matrix.LUdecompose(a,n)
	return function(b) return Matrix.LUbsolve(a,indx,b) end
end

function Matrix.LUsolve(a,b) -- Solve by LU decomposition
	local a,indx = Matrix.LUdecompose(a)
	return Matrix.LUbsolve(a,indx,b)
end

function Matrix.inv(m)
	local n,b = #m, {}
	local mc,ai = Matrix.new(m),Matrix.new(n,n)
	local fi = Matrix.LUfunction(mc,n)
	for i=1,n do
		for j=1,n do b[j] = 0 end
		b[i] = 1
		b = fi(b)
		for j=1,n do ai[j][i] = b[j] end
	end
	return ai
end

function Matrix.det(a,n) -- Determinant of matrix.
	n = n or #a
	local a,indx,d = Matrix.LUdecompose(a,n)
	for i=1,n do d = d*a[i][i] end
	return d
end

function Matrix.gauss(a,bb,n) -- a and bb changed by function
	n = n or #bb
	local b = {} -- bb can be a column matrix or a row vector
	if type(bb[1])=='table' then for i=1,n do b[i]  = bb[i][1] end
	else b = bb end
	return spgauss(a,b,n) -- spgauss() changes a and b
end

function Matrix.trace(a,n) -- Trace of a matrix
	n = n or #a
	local tx = a[1][1]
	for i=2,n do tx = tx + a[i][i] end
	return tx
end

function Matrix.eigeneq(a,n) -- Set up eigenvalue equation
	n = n or #a
	local p,pm,px = {}, Matrix.new(a)
	p[n+1] = -1
	for i=1,n do
		px = Matrix.trace(pm)/i
		for j=1,n do pm[j][j] = pm[j][j] - px end
		p[n-i+1] = px
		pm = a*pm -- Matrix multiplication, a and pm both matrices
	end
	return p -- form a0 + a1*x + a2*x^2 ...
end

function Matrix.eigenvectors(ai,eval) -- Calculate Matrix Eigenvectors
	local n = #ai
	if type(ai)~='Matrix' then ai = Matrix.new(ai) end
	eval = eval or Matrix.eigenvalues(ai,n) -- Get eigenvalues if not input
	local evect,b,big,tmp = Matrix.new(n),{}
	b[n] = 1 -- Assume one value for solution
	for j=1,n do -- Loop over eigenvalues
		a = Matrix.new(ai) -- Make copy for use in gauss()
		for i=1,n-1 do
			b[i] = -a[i][n]
			a[i][i] = a[i][i] - eval[j] -- Subtract eigenvalue
		end
		gauss(a,b,n-1) -- Solve for one eigenvector, a and b changed in gauss()
		big = 0
		for i=1,n do -- Normalize so largest element is +/- 1.000
			tmp = abs(b[i])
			if tmp> big then big = tmp end
		end
		for i=1,n do evect[i][j] = b[i]/big end
	end
	return evect,eval -- Also return eigenvalues for possible use
end
	
function Matrix.eigenvalues(a,b,n) -- Set up equation and solve using roots()
	require"Polynomial" -- Only load if needed for roots()
	local d
	if type(a)=='table' then a = Matrix.new(a) end
	if type(b)=='table' then b = Matrix.new(b) end
	if b==nil then d = a 
	elseif type(b)=='number' then d,n = a,b 
	elseif type(b)=='Matrix' then d = b^-1*a 
	else print('b must be a matrix in eigenvalues') end
	local rts = Polynomial.roots(Matrix.eigeneq(d,n))
	table.sort(rts,function(a,b) 
		if type(a)=='Complex' or type(b)=='Complex' then return Complex.abs(a)>Complex.abs(b)
		else return a>b end end) -- sort eigenvalues
	return rts,d -- Sorted most positive to most negative
end
			

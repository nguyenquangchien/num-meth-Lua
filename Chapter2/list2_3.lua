-- File list2.3 -- Simple code for implementing complex number operations

Cmp = {} -- Complex table
mtcmp = { -- metatable for Complex numbers
	__add = function(c1,c2) -- Add two Complex numbers
		return Cmp.new(c1.r+c2.r, c1.i+c2.i)
	end,
	__sub = function(c1,c2) -- Subtract two Complex numbers
		return Cmp.new(c1.r-c2.r, c1.i-c2.i)
	end,
	__mul = function(c1,c2) -- Multiple two Complex numbers
		return Cmp.new(c1.r*c2.r-c1.i*c2.i, c1.i*c2.r+c2.i*c1.r)
	end,
	__div = function(c1,c2) -- Divide two Complex numbers
		local d = c2.r*c2.r+c2.i*c2.i
		return Cmp.new((c1.r*c2.r+c1.i*c2.i)/d, (c1.i*c2.r-c2.i*c1.r)/d)
	end,
	__unm = function(c1) -- Negative of complex number
		return Cmp.new(-c1.r, -c1.i)
	end,
	__tostring = function(c1)
		return '('..c1.r..') + j('..c1.i..')'
	end
} -- End metatable functions
Cmp.new = function(r,i) -- Define new complex number
	local c1 = {r = r, i = i}
	return setmetatable(c1, mtcmp)
end



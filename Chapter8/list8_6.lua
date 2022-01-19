-- /* File list8_6.lua */ -- Example of obtaining Student's t values

require"prob"

pbv = {.5,.9,.95,.99,.999} -- Define probabilities
n = {1,2,4,6,10,20,40,120,10000} -- Define DOF's
ts = {}
for i=1,#n do
	for j=1,#pbv do
		ts[j] = ttable(n[i],pbv[j]) -- Get t values for table
	end
	print(ts[1],ts[2],ts[3],ts[4])
end
 
 -- Following code is actually in prob.lua
ttable = function(n,pbv)
	local Atf = function(x)
		return Atstud(x,n) - pbv
	end
	return newton(Atf,0)
end

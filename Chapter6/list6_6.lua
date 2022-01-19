-- /* File list6_6.lua */ -- Normal curve of error function definition
require"intp"; require"intg"

ncerr = function() -- Error function using 22 data points
	local x1,y1,x2,y2,x2t,y2t = {},{},{},{},{},{}
	local f1 = function(x) return math.exp(-x^2/2) end -- Normal curve
	local f2 = function(x) return f1(1/x)/x^2 end -- Folded normal curve
	local i = 2
	x1[1],y1[1] = 0,0
	for xi=.1,1.01,.1 do -- Integral from 0 to 1
		x1[i],y1[i] = xi, intg(xi-.1,xi,f1)/math.sqrt(2*math.pi)+y1[i-1]
		i = i+1
	end
	i = 2
	x2[1],y2[1] = 1,y1[#x1]
	for xi=.9,-.01,-.1 do -- Integral from 1 to infinity
		x2[i],y2[i] = xi, intg(xi,xi+.1,f2)/math.sqrt(2*math.pi)+y2[i-1]
		i = i+1
	end
	local n = #x2 -- Need to reverse table entries
	for i=1,n do y1[i],y2[i] = .5*y1[i]/y2[n],.5*y2[i]/y2[n] end -- Exactly .5 max
	for i=1,n do x2t[i],y2t[i] = x2[n+1-i],y2[n+1-i] end
	x2,y2 = x2t,y2t -- Now have increasing x2 values
	return function(x) -- Return Normal cumulative error function
		if x>=0 then
			if x<=1 then return intp(x1,y1,x)+0.5
			else return intp(x2,y2,1/x)+0.5 end
		else 
			if x<-1 then return 0.5-intp(x2,y2,-1/x)
			else return 0.5-intp(x1,y1,-x) end
		end
	end
end
-- Now test function
i,x2t,y2t = 1,{},{} -- Now test function
errfcn = ncerr()
for x=-3,3,.5 do
	x2t[i],y2t[i],i = x, errfcn(x), i+1
	print(x,y2t[i-1])
end
write_data("list6_6.dat",x2t,y2t)

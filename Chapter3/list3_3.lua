--/* File 3_3.lua Successive substitution example */

require"Complex"; require"ssroot"

function g2(x) -- Negative root function
	return -((3-2*x+4*x^2)/5)^.25
end
function g3(x) -- Complex root
	return j*((3-2*x+4*x^2)/5)^.25
end
function g4(x)
	return math.atan(1/x)
end

print(ssroot(g2,0)) -- Real initial guess
print(ssroot(g3,1)) -- Real initial guess
print(ssroot(g4,1)) -- Real initial guess

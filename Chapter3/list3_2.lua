--/* File 3_2.lua Successive substitution code */

function g1(x) -- Positive root function
	return ((3-2*x+4*x^2)/5)^.25
end
function g2(x) -- Negative root function
	return -((3-2*x+4*x^2)/5)^.25
end

function ssroot(g,x) -- Successive Substitution function
	local xold = x
	local itt
	for i=1,NT do
		x = g(x)
		if abs((xold-x)/x)<err then itt=i; break end
		xold = x
	end
	return x,itt
end
setfenv(ssroot,{abs=(Complex or math).abs,NT=200,err=1.e-10})

root1,n1 = ssroot(g1,0)
print(root1,n1,(1-root1)/root1)
root2,n2 = ssroot(g2,0)
print(root2,n2,(-1.2326393307127-root2)/root2)


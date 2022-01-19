--/* File ssroot.lua Successive substitution code */

function ssroot(g,x) -- Successive Substitution function
	local xold = x
	local itt=NT
	for i=1,NT do
		x = g(x)
		if abs((xold-x)/x)<err then itt=i; break end
		xold = x
	end
	return x,itt
end
setfenv(ssroot,{abs=(Complex or math).abs,NT=200,err=1.e-10})


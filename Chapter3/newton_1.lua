--/* File newton_1.lua */
-- Solves a nonlinear equation f(x) = 0 for the solution value x

function newton(f,x)
	local nend,cx,fx = NTT
	for i=1,NTT do
		local dx = x*FACT
		if dx==0.0 then dx = FACT end -- Protect against zero
		fx = f(x)
		cx = fx - f(x+dx)
		cx = fx*dx/cx -- Correction to x value
		x = x + cx
		dx = abs(cx) + abs(x)
		if dx==0.0 or abs(cx/dx)<ERR then nend = i; break end
	end
	return x,nend,cx,fx
end
setfenv(newton, {abs=(Complex or math).abs,FACT=1.e-6,ERR=1.e-8,
	NTT=200}) -- Make all variables and constants local to newton

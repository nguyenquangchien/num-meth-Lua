--/* File newton.lua */
-- Solves a nonlinear equation f(x)=0 for the solution value x, with step limit per iteration

function newton(f,x,step) -- step may be omited in calling
	if step then if step>0.0 and step<1.0 then step = 1./step end end
	nend = NTT
	for i=1,NTT do
		dx = x*FACT
		if dx==0.0 then dx = FACT end -- Protect against zero
		fx = f(x)
		cx = fx - f(x+dx)
		cx = fx*dx/cx -- Correction to x value
		if step and step~=0.0 then -- Skip if no step size specified
			if step<0.0 then -- Limit absolute step size
				if abs(cx)>-step then cx = -step*(cx/abs(cx)) end
			else -- Limit percentage increase in solution value
				dx = 1. + cx/x
				if dx>step then 
					cx = (step-1)*x 
				elseif dx<1/step then 
					cx = (1/step-1)*x 
				end
			end
		end -- End of step adjustment
		x = x + cx
		if nprint then print("Iteration ",i,x,dx) end -- Print values at each iteration
		dx = abs(cx) + abs(x)
		if dx==0.0 or abs(cx/dx)<ERR then nend = i; break end
	end
	return x,nend,cx/dx
end
-- Make all variables and constants local to newtonl
setfenv(newton, {abs=(Complex or math).abs,FACT=1.e-6,ERR=1.e-9,NTT=200,nprint=nil,print=print})

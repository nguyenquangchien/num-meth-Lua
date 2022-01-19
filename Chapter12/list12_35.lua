-- File list12_35.lua --
-- Example of h-2h algorithm for PDEs

require"odeiv"

u80, u160 = {},{} -- Define arrays for solutions
x80, x160 = {}, {} -- Linear spatial arrays
for i=1,81 do 
	u80[i] = {}; x80[i] = (i-1)/81 
end
for i=1,161 do 
	u160[i] = {}; x160[i] = (i-1)/160
end
read_data('list12_34u_80.dat',u80)
read_data('list12_34u_160.dat',u160)
err = {}; maxerr = 0.0
for i=1,81 do 
	s80 = {x80, u80[i]}
	s160 = {x160, u160[1+(i-1)*2]}
	erri = odeerror(s80,s160)[2] -- [2] for error values
	for j=1,#erri do -- Maximum error
		if math.abs(erri[j])>math.abs(maxerr) then
			maxerr, im, jm  = erri[j], i, j
		end
	end
	err[i] = erri
end
splot(err); print('Maximum error =', maxerr, 'at i,j =', im,jm)
write_data('list12_35.dat',err) -- Save error array



-- /* list11.32a.lua */
-- Solution of nonlinear BV problem with non-uniform spatial grid
require"odeiv"; require"nlstsqi"; require"intp"
getfenv(nlstsqi).nprint=1 -- print iteration values

x,ct,temp = {},{},{}; read_data('list11_24.dat',x,ct,temp) -- Data points
plot(x,ct);plot(x,temp)
xx,cy,ty = {},{},{}; j,nx = 1,#x
for i=0,48,4 do 
	xx[j],cy[j],ty[j] = i,intp(x,ct,i),intp(x,temp,i)
	j = j+1
end
--x[j-1],cy[j-1],ty[j-1] = x[nx], ct[nx],temp[nx]
write_data('list11_32c_in.dat',xx,cy,ty)
nd = #xx
--math.random()
fc,ft = 0.04,0.002
fc,ft = fc*(cy[1]+cy[nd]),ft*(ty[1]+ty[nd])
for i=1,nd do 
	cy[i] = cy[i] + (math.random()-0.5)*fc
	ty[i] = ty[i] + (math.random()-0.5)*ft
	--ty[i] = ty[i] + (math.random()-0.5)*ft-100
	--cy[i] = cy[i] + (math.random()-0.5)*fc+.02

end
plot({xx,cy},{x,ct});plot({xx,ty},{x,temp})
write_data('list11_32b_in.dat',xx,cy,ty)


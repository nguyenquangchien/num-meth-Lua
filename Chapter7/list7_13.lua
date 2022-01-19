-- /* File list7_13.lua */ -- Fit to two gaussian peaks with background

require"nlstsq"

g1,x1 = {},{}
read_data('gauss3.dat',g1,x1) -- Read data from file
plot(x1,g1)

gfit = function(xa,c) -- Define function for two peaks
	local x = xa[2]
	return c[1]*math.exp(-c[2]*x) + c[3]*math.exp(-(x-c[4])^2/c[5]^2) +
		c[6]*math.exp(-(x-c[7])^2/c[8]^2) - xa[1]
end

c = {100,.01,90,110,25,80,150,25} -- c[7] between 145 and 160, c[4] between 100 and 120
del,err,nmax = nlstsq({g1,x1},fw,gfit,c) -- Call fitting program

print('RMS error =',err,'#Iter =',nmax)
nd = #x1
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end -- Print fitting coefficients
gcalc,gp1,gp2,gb = {},{},{},{} -- Set up arrays for total and individual peaks
for i=1,nd do gcalc[i] = gfit({0,x1[i]},c) end -- Total fitted function
c1,c3,c6 = c[1],c[3],c[6]
c[3],c[6] = 0,0 -- Set peaks to zero
for i=1,nd do gb[i] = gfit({0,x1[i]},c) end -- Background only
c[1],c[3]  = 0,c3
for i=1,nd do gp1[i] = gfit({0,x1[i]},c) end -- Peak one only
c[3],c[6] = 0,c6
for i=1,nd do gp2[i] = gfit({0,x1[i]},c) end -- Peak two only
plot(x1,g1,gcalc,gp1,gp2,gb) -- plot individual components
write_data('list7_13.dat',x1,g1,gcalc,gp1,gp2,gb) -- Save data

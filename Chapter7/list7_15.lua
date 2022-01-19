-- /* File list7_15.lua */ -- general data fitting with polynomial functions

require"nlstsq"

xd,yd = {},{}
read_data('PhDgrads.dat',xd,yd)
--read_data('gauss3.dat',yd,xd) -- Another data set
-- read_data('list7_1.dat,xd,yd) -- Yet another data set

polyft = function(x,c) -- Define general polynomial function to fit data
	local nc,xx,sum = #c, x[2], 0
	for i=nc,1,-1 do sum = sum*xx + c[i] end
	return sum - x[1]
end

nc = 5 -- Specify number of polynomial coefficients
c = {}; for i=1,nc do c[i] = 0 end -- zero initial guesses
del,err,nmax = nlstsq({yd,xd},fw,polyft,c) -- Call fiting 
print('RMS error =',err,'#Iter =',nmax) -- print results
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
xcalc,ycalc = {},{}
x1,x2,i = 0,60, 1; dx = (x2-x1)/100
for x=x1,x2,dx do 
	xcalc[i],ycalc[i] = x, polyft({0,x},c)
	i = i+1
end
write_data('list7_15.dat',xcalc,ycalc,xd,yd)
plot(xcalc,ycalc);plot(xd,yd)
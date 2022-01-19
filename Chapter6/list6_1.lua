-- /* File list6_1.lua */ -- Linear interpolation 

lintp = function(xd,yd,x,deriv) -- Linear interpolation with first derivative
	local n,y1,y2,x1,x2 = #xd
	while x>xd[il+1] and il<n-1 do il = il+1 end -- Find interval for x
	while x<xd[il] and il>1 do il = il-1 end
	y1,y2,x1,x2 = yd[il],yd[il+1],xd[il],xd[il+1]
	y2 = (y2-y1)/(x2-x1)
	if deriv==1 then return y2 end
	return y1 + y2*(x-x1)
end
setfenv(lintp,{deriv=0,il=1})

lintpf = function(xd,yd) -- Function with built-in data table
	local nd,xdd,ydd = #xd, {}, {} -- So user can't change data
	if xd[nd]<xd[1] then for i=1,nd do xdd[i] = xd[nd+1-i] end
	else for i=1,nd do xdd[i] = xd[i] end end
	for i=1,nd do ydd[i] = yd[i] end
	return function(x,deriv) return lintp(xdd,ydd,x,deriv) end
end
	
xd = {0,1,2,3,4,5}; yd = {3,4,4.5,4,3,2.5}
print('value at x = 2.3 is',lintp(xd,yd,2.3))
print('value at x = 5.5 is',lintp(xd,yd,5.5))

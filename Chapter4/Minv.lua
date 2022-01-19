-- /* File Minv.lua */
-- Inverts a matrix by full Gauss elimination

minv = function (ax,n)
	local i,j,k,jp1,ix,jx,imax,jret
	local ap,row,rowj
	local a = {}
	n = n or #ax
	if n < 2 then -- n=1, special case
		if n<1 then return -1 end
		if ax[1][1] ~= 0.0 then a[1] = {1/ax[1][1]}; return a end
		return -1
	end
	if n==2 then -- n=2, special case
		a[1],a[2] = {},{}
		ap = ax[1][1]*ax[2][2] - ax[1][2]*ax[2][1]
		a[1][1],a[2][2] = ax[2][2]/ap,ax[1][1]/ap
		a[1][2],a[2][1] = -ax[1][2]/ap,-ax[2][1]/ap
		return a
	end
	for i=1,n do -- Make internal copy of matrix ax
		row,rowj = ax[i], {}
		a[i] = rowj
		for j=1,n do -- Save and add columns
			rowj[j],rowj[j+n] = row[j], 0 
		end
		rowj[n+i] = 1
	end
	n2,jret = 2*n, 1
	for i=1,n do -- Make largest element in each row unity
		ap,row = 0.0, a[i]
		for j=1,n do -- Find largest element
			if ap<abs(row[j]) then ap=abs(row[j]) end
		end
		if ap==0.0 then ap=1.0 end
		for j=1,n2 do row[j]=row[j]/ap end -- Divide by largest element
	end
	for j=1,n do -- Find row with maximum value in column j
		jp1,ap = j+1, 0.0
		for i=j,n do
			row=a[i]
			if abs(row[j]) > abs(ap) then
				ap,imax = row[j], i
			end
		end -- ap is now largest remaining column element at row i
		if abs(ap) <= EPS then jret=0 end
		if ap==0.0 then return -2 end
		for k=j,n2 do -- Swap rows if needed and make leading element unity
			row=a[j]
			if imax~=j then a[imax][k],row[k]=row[k],a[imax][k] end
			row[k]=row[k]/ap
		end
		rowj=a[j]
		for ix=1,n do -- Now perform element eliminations for off-diagonal elements
			if ix~=j then 
				row=a[ix]; rowjv = row[j]
				for jx=jp1,n2 do row[jx]=row[jx]-rowjv*rowj[jx] end
			end	
		end
	end
	
	for j=1,n do -- Shift solution to first n columns of a
		row = a[j]
		for k=1,n do row[k],row[k+n] = row[k+n],nil end
	end
	return a,jret
end
		
setfenv(minv,{abs=(Complex or math).abs,EPS=1.e-12})
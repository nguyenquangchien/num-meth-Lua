-- /* File gauss.lua */

gauss = function (a,b,n)
	n = n or #a
	local jret = 1
	if n < 2 then -- Low dimension matrix -- Special case
		if n<1 then return b,-1 end
		row = a[1]
		if row[1] ~= 0.0 then b[1]= b[1]/row[1] return b,1 end
		return b,-1
	end
	for i=1,n do  -- Find largest value in each row and divide row elements by largest value
		row,ap = a[i], 0.0
		for j=1,n do ap = max(ap,abs(row[j])) end -- ap = largest value
		if ap==0.0 then ap=1.0 end
		for j=1,n do row[j]=row[j]/ap end -- Divide by largest value
		b[i]=b[i]/ap
	end
		
	for j=1,n do -- Elimination for each row in matrix
		jp1,ap = j+1, 0.0
		for i=j,n do -- Find largest value in column j, going from row j to row n
			am = a[i][j]
			if abs(am) > abs(ap) then ap,imax = am,i end
		end -- At end ap = largest value in column j to be eliminated and imax has row number
		if abs(ap) <= eps then 
			jret=0 -- Probably singular matrix with no solution
			if ap==0.0 then return b,-2 end -- Singular matrix with no solution
		end
		if imax~=j then -- Swap rows of a and b
			a[imax],a[j],b[imax],b[j] = a[j],a[imax],b[j],b[imax]
		end
		row,b[j] = a[j],b[j]/ap
		for k=j,n do row[k] = row[k]/ap end 
		if j<n then -- Eliminate elements except for last row
			rowj=a[j] -- Selece row j values 
			for ix=jp1,n do -- Step rows from j+1 to row n
				row = a[ix]; rowij = row[j] -- Select row ix values
				for jx=jp1,n do -- Step columns from j+1 to column n
					row[jx]=row[jx]-rowij*rowj[jx] 
				end
				b[ix] = b[ix] - b[j]*rowij
			end
		end
	end
	for j=n-1,1,-1 do -- Back substitute from row n-1 to row 1
		for k=n,j+1,-1 do -- Known values from n to j
			b[j] = b[j] - a[j][k]*b[k]
		end
	end
	return b,jret
end
setfenv(gauss,{abs=(Complex or math).abs,max=math.max,eps=1.e-12})


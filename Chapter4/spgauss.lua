-- /* File spgauss.lua */

spgauss = function (a,b,n) -- Gauss elimination for sparse matrix a[][] and b[]
	local nel = {{},{}} -- arrays for usage statistics
	local jret,jprint,jpold,nct,ug1 = 1,1,1,0,0 -- print, usage parameters
	n = n or #a 
	if type(a[1][1])=='table' then for m=1,n do a[m][m]=a[m][m][1]+a[m][m][2] end end
	if n < 2 then -- Low dimension matrix -- Special case
		if n<1 then return b,-1 end; row = a[1]
		if row[1] ~= 0.0 then b[1]= b[1]/row[1] return b,1 end
		return b,-1
	end
	for i=1,n do  -- Find largest value in each row and divide row elements by largest value
		ap,row=0.0, a[i]
		for j,v in pairs(row) do ap = max(ap,abs(v)) end -- ap = largest value in row i
		if ap==0.0 then ap=1.0 end; ap = 1/ap
		for j,v in pairs(row) do row[j] = v*ap end
		b[i] = (b[i] or 0)*ap -- All elements of b[] exist after this loop
	end
	for j=1,n do -- Elimination for each row in matrix
		jp1,ap = j+1, 0.0
		for i=j,n do -- Find largest value in column j, going from row j to row n
			am = a[i][j] 
			if am and abs(am)>abs(ap) then ap,imax = am,i end
		end -- At end ap = largest value in column j to be eliminated and imax has row number
		if abs(ap) <= eps then -- Probably singular matrix with no solution
			jret=0; if ap==0.0 then return -2 end -- Singular matrix with no solution
		end
		if imax~=j then -- Swap rows of a and b
			a[imax],a[j],b[imax],b[j] = a[j],a[imax],b[j],b[imax] 
		end
		row,b[j] = a[j],b[j]/ap; row[j] = nil -- Normalize new row
		for k,v in pairs(row) do row[k] = v/ap end 
		if j<n then -- Eliminate elements if not last row
			rowj=a[j] -- Selece row j values 
			for ix=jp1,n do -- Step rows from j+1 to row n
				row=a[ix]; rowij = row[j] -- Select row ix values
				if rowij then -- Non nil value
					row[j] = nil
					for jx,v in pairs(rowj) do -- Step through columns of row j
						row[jx] = (row[jx] or 0) - rowij*v
					end
					b[ix] = b[ix] - b[j]*rowij
				end
			end
		end
		if usage==2 then -- Collect usage statistics
			ug1= 0
			for _,_ in pairs(a[j]) do nct = nct+1 end -- Row n
			for jj=j+1,n do -- Rows beyond n
				for _,_ in pairs(a[jj]) do ug1 = ug1+1 end
			end 
			nel[1][j], nel[2][j] = j, nct + ug1
		end
		if nprint then -- Print at each 100 rows
			jprint = floor(j/100)
			if jprint==jpold then
				jpold=jprint+1; print("Completed row ",j," in spgauss");io.flush() 
				if usage==2 then print("Number of matrix elements =",nct+ug1) end
			end
		end
	end
	for j=n-1,1,-1 do -- Back substitute from row n-1 to row 1
		row = a[j]; for k,v in pairs(row) do b[j] = b[j] - v*b[k] end
	end
	if usage==1 then -- Collect usage statistics at end
		nct = 0
		for jj=1,n do
			for _,_ in pairs(a[jj]) do nct = nct+1 end
		end
		nel[1][1],nel[2][1] = n, nct
	end
	return b,jret,nel
end
setfenv(spgauss,{abs=(Complex or math).abs,min=math.min,max=math.max,floor=math.floor,
io=io,print=print,pairs=pairs,type=type,nprint=nil,usage=nil,eps=1.e-12})


Spmat = {}
local spmat_mt = {__index = Spmat, __type = 'SPMAT'}
Spmat.new = function(i, j, dir) -- Sparce matrix function
	dir = dir or 1 -- 1 means x-axis is scanned first
	local mat = {}
	local m = i*j
	mat[0] = {abs(i),abs(j),dir} -- store properties
	while m>0 do mat[m],m = {},m-1 end
	return setmetatable(mat,spmat_mt)
end
setfenv(Spmat.new,{abs=math.abs,setmetatable=setmetatable})

Spmat.print = function(mat)
	local i,j = mat[0][1],mat[0][2]; j = i*j
	local st
	for i = 1,math.min(300,j) do 
		st = string.rep('_',math.min(j,300))
		for k in pairs(mat[i]) do 
			st = string.sub(st,1,k-1)..'+'..string.sub(st,(k+1))
		end
		print(st)
	end
end

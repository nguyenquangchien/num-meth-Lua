-- /* File list4_27.lua */
-- Tests of eigenvalues and eigenvectors

require"Polynomial"
require"Matrix"

n = 10; A = Matrix.new(n)
for i=1,n do
	for j=1,n do
		if j==i-1 then A[i][j] = -1
		elseif j==i then A[i][j] = 2
		elseif j==i+1 then A[i][j] = -1
		else A[i][j] = 0 end
	end
end	
print(A)

ev,rts = Matrix.eigenvectors(A)
print('Eigenvectors \n',ev,'\nEigenvalues')
table.foreach(rts,print)

A = Matrix.new{{300,-200,0},{-200,500,-300},{0,-300,300}}
B = Matrix.new{{4,-1,1},{-1,6,-4},{1,-4,5}}
rts = Matrix.eigenvalues(A,B)
print('Eigenvalues for A and B matrix problem')
table.foreach(rts,print)

a = Matrix.new{{2+j,-2-j,3},{2,j,-3+2*j},{j,4-j,3*j}}
rts = Matrix.eigenvalues(a)
print('Eigenvalues for complex matrix')
table.foreach(rts,print)
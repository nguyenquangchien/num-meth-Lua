-- /* File list8_7.lua */

require"prob"

x,y = {},{}
f = 1.e-3; ds = 'm'
read_data('paper.txt',y)
nd = #y
for i=1,nd do y[i] = f*y[i] end
m,std = stats(y)
print('mean,std = ',engr(m,ds),engr(std,ds))
stderr = std/math.sqrt(nd-1)
 
pbv = {.6827,.9545,.9,.95,.99,.999}
for i=1,#pbv do
	ts = ttable(nd-1,pbv[i])
	dm = ts*stderr
	print(engr((m-dm),ds)..' < mean <'..engr((m+dm),ds)..' at '..100*pbv[i]..' percent confidence level')
end
print(engr((m-stderr),ds)..' < mean <'..engr((m+stderr),ds).. ' for standard error bounds ')
print(engr((m-2*stderr),ds)..' < mean <'..engr((m+2*stderr),ds).. ' for 2 standard error bounds ')

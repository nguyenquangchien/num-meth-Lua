-- /* File list8_8.lua */ Example of prob limits on std

require"prob"

x,y = {},{}
f = 1.e-3
ds = 'm'
read_data('paper.txt',y)
nd = #y
for i=1,nd do y[i] = f*y[i] end
m,std = stats(y)
print('mean,std = '..engr(m,ds)..engr(std,ds))
stdf,fdfc = std*math.sqrt(nd-1), 1/math.sqrt(2*nd)
 
pbv = {.6827,.9545,.9,.95,.99,.999}
for i=1,#pbv do
	tchL,tchH = tchisq(nd-1,pbv[i])
	varL,varH = stdf/math.sqrt(tchL),stdf/math.sqrt(tchH)
	print(engr(varH,ds)..' < std <'..engr(varL,ds)..' at '..100*pbv[i]..' percent confidence level')
end
print(engr(std*(1-fdfc),ds)..' < std <'..engr(std*(1+fdfc),ds)..' for standard error bounds ')
print(engr(std*(1-2*fdfc),ds)..' < std <'..engr(std*(1+2*fdfc),ds)..' for 2 standard error bounds ')

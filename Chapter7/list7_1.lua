-- /* File listing7.1.lua */ -- Some definitions of statistical functions
require"prob"

x7_1 = {0.0,0.2,0.8,1,1.2,1.9,2,2.1,2.95,3}
y7_1 = {0.01,0.22,0.76,1.03,1.18,1.94,2.01,2.08,2.9,2.95}
x7_2 = {1,2,3,4,5,6,7,8,9,10}
y7_2 = {6,9.2,13,14.7,19.7,21.8,22.8,29.1,30.2,32.2}

a1,b1,rsq1 = clinear(x7_1,y7_1) -- Linear least squares routine
print(a1,b1,rsq1)
a2,b2,rsq2 = clinear(x7_2,y7_2)
print(a2,b2,rsq2)
yy7_1 = {}
for i=1,#x7_1 do yy7_1[i] = a1*x7_1[i] + b1 end
write_data('list7_1a.dat',x7_1,y7_1,yy7_1)
yy7_2 = {}
for i=1,#x7_2 do yy7_2[i] = a2*x7_2[i] + b2 end
write_data('list7_1b.dat',x7_2,y7_2,yy7_2)
a3,b3,rsq3 = clinear(y7_2,x7_2)
print(a3,b3,rsq3)
xx7_3 = {}
for i=1,#y7_2 do xx7_3[i] = a3*y7_2[i] + b3 end
write_data('list7_1c.dat',x7_2,y7_2,xx7_3)


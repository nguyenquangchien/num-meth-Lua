-- /* File list5_11.lua */ -- Singular integrals 

require"intg"
getfenv(intg).ERROR=1.e-6

f1 = function(x) return x^-.5 end
iv,n,err = intg(0,4,f1)
ivv = 4
print(n,iv,ivv,err,(iv-ivv)/ivv)

f2 = function(x) return ((x+1)^2)^(-1/3) end
iv1,n1,err1 =  intg(-2,-1,f2)
iv2,n2,err2 = intg(-1,2,f2)
iv = iv1 + iv2
ivv = 3*(3^(1/3)+1 )
print(n1+n2,iv,ivv,err1+err2,(iv-ivv)/ivv)

f3 = function(x) return 1/(4 - x^2)^.5 end
iv,n,err = intg(-2,0,f3) 
ivv = -math.asin(-1)
print(n,iv,ivv,err,(iv-ivv)/ivv)



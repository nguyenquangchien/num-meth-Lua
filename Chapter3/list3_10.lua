 --/* File list3_10.lua */ -- equation with multiple roots
 
require "newton"
getfenv(newton).nprint = 1 -- Print iterative values

function f(x) return (x-2)^m end

m = 10; newton(f,0)
m = 11; newton(f,0)
function dup(min)
 local a = type(min)
 if min == nil then
    a = "." 
 return a..','..a
 else
 return min..','..min
 end
end

function domain(a)
 local b = type(a)
 if a == nil then
    b = "."
 return "-1,1"
 else
 return "-1,-1"
 end
end

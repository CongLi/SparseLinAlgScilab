function [Min, Max] = SEE(b)
   Min = 10000000
   Max = 0
   for i=1:length(b)
     if b(i,1)<Min then
        Min = b(i,1)
     end
     if b(i,1)>Max then
        Max = b(i,1)
     end
   end
endfunction

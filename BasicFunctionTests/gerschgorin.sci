function [Min, Max] = gerschgorin(A, n)
    Max = 0
    Min = max(A)
    for i=1:n
       tmp = 0
	   for j=1:n
	       if j~=i then
	        tmp = tmp + abs(A(i, j))
	       end
	   end
	   if (tmp+A(i,i))>Max then
	       Max = tmp+A(i,i)
    	end
    	if (A(i,i)-tmp)<Min then
            Min = A(i,i)-tmp
    	end
    end
endfunction

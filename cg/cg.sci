function [x, hist] = cg(A, b, Size)
    x = zeros(b)
    r = b
    p = r
    hist(1,:) = [0,(r' * r)]
    for i=1:(2*Size)
        rrp = (r' * r)
        alp = (r' * r) / (p' * A * p)
        r = r - alp * A * p
        x = x + alp * p
        hist((i+1),:) = [i,(r' * r)]
        if ((r' * r) / (b' * b) < 1e-20) then
            break
        end
        bet = (r' * r) / rrp
        p = r + bet * p 
    end
endfunction

function cg_main(filename, A, Size, b, n)
    [x,hist]=cg(A,A*b(1:Size,1), Size);
    for i = 2:n
        [x,tmp]=cg(A,A*b(1:Size,i), Size);
	if hist($,1) < tmp($,1) then
            hist = tmp
	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
            hist = tmp
        end
    end
    str=sprintf("cg_%s_sample=%d.txt",filename,n)
    fprintfMat(str,hist,"%10.7e")
    disp(hist($,1))
endfunction

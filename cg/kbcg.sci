function [x, hist] = kbcg(A, b, k, Size)
    x = zeros(b)
    r = b
    p = r
    hist(1,:) = [0,(r' * r)]
    for i=1:(floor(2*Size/(k+1))+10)
        for j=1:k+2
            if j == 1 then
                S(:,1) = r
		S(:,2) = p
            else
                S(:,2*j-1) = A * S(:,2*j-3)
                S(:,2*j) = A * S(:,2*j-2)
            end
        end
	Q = S' * S
        wr = zeros(2*k+4, 1)    
	wr(1,1) = 1
        wp = zeros(2*k+4, 1)
	wp(2,1) = 1
        wap = zeros(2*k+4, 1)
        wx = zeros(2*k+4, 1)
	for j=1:k+1
	    Gamma = (wr' * Q * wr)
      	    wap(3:2*k+4,1) = wp(1:2*k+2,1)
	    alpha = Gamma / (wp' * Q * wap)
	    wr = wr - alpha * wap
	    wx = wx + alpha * wp
	    Beta = (wr' * Q * wr) / Gamma
	    wp = wr + Beta * wp
        end
	p = S * wp
	r = S * wr
        x = x + S * wx
        hist((i+1),:) = [i*(k+1),(r' * r)]
        if (norm(r) / norm(b) < 1e-10) then
            break
        end
    end
endfunction

function kbcg_main(filename, A, Size, b, n, k)
    [x,hist]=kbcg(A,A*b(1:Size,1),k, Size);
    for i = 2:n
        [x,tmp]=kbcg(A,A*b(1:Size,i),k, Size);
	if hist($,1) < tmp($,1) then
            hist = tmp
	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
            hist = tmp
        end
    end
    str=sprintf("kbcg_%s_k=%d_sample=%d.txt",filename,k,n)
    fprintfMat(str,hist,"%10.7e")
    disp([k,hist($,1)])
endfunction

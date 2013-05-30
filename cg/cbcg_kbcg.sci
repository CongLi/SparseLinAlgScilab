function [x, hist] = cbcg_kbcg(A, b, Min, Max, k, Size)
    rho0 = Max-Min
    rho1 = Max+Min
    x = zeros(b)
    r = b
    p = r
    hist(1,:) = [0,(r' * r)]
    for i=1:(floor(2*Size/(k+1))+10)
        for j=1:k+2
            if j == 1 then
                B(:,1) = r
                B(:,2) = p
            elseif j == 2 then
                B(:,3) = (2/rho0)* A * r - (rho1/rho0) * r
                B(:,4) = (2/rho0)* A * p - (rho1/rho0) * p
            else
                B(:,2*j-1) = (4/rho0) * A * B(:,2*j-3) - (2*rho1/rho0) * B(:,2*j-3) - B(:,2*j-5)
                B(:,2*j)   = (4/rho0) * A * B(:,2*j-2) - (2*rho1/rho0) * B(:,2*j-2) - B(:,2*j-4)
            end
        end
	Q = B' * B
        wr = zeros(2*k+4, 1)   
	wr(1,1) = 1
        wp = zeros(2*k+4, 1)
	wp(2,1) = 1
        wap = zeros(2*k+4, 1)
        wx = zeros(2*k+4, 1)
	for j=1:k+1
            for l=1:2*k+4
                if l==1|l==2 then
                    wap(l, 1) = wp(l+2,1)/2
                elseif l==3|l==4 then
                    wap(l, 1) = wp(l-2,1)+wp(l+2,1)/2
                elseif (l==2*k+3)|(l==2*k+4) then
                    wap(l, 1) = wp(l-2,1)/2
                else 
                    wap(l, 1) = (wp(l-2,1) + wp(l+2,1))/2
                end
            end
            wap = (rho0/2)*wap+(rho1/2)*wp
	    Gamma = (wr' * Q * wr)
	    alpha = Gamma / (wp' * Q * wap)
	    wr = wr - alpha * wap
	    wx = wx + alpha * wp
	    Beta = (wr' * Q * wr) / Gamma
	    wp = wr + Beta * wp
        end
	p = B * wp
	r = B * wr
        x = x + B * wx
        hist((i+1),:) = [i*(k+1),(r' * r)]
//        disp([i*(k+1),r' * r])
        if (norm(r) / norm(b) < 1e-10) then
            break
        end
    end
endfunction

function cbcg_kbcg_main(filename, A, Min, Max, Size, b, n, k)
    [x,hist]=cbcg_kbcg(A,A*b(1:Size,1),Min,Max,k, Size);
    for i = 2:n
        [x,tmp]=cbcg_kbcg(A,A*b(1:Size,i),Min,Max,k, Size);
	if hist($,1) < tmp($,1) then
            hist = tmp
	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
            hist = tmp
        end
    end
    str=sprintf("cbcg_kbcg_%s_k=%d_sample=%d.txt",filename,k,n)
    fprintfMat(str,hist,"%10.7e")
    disp([k,hist($,1)])
endfunction

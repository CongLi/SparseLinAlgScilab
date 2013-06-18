function [x, hist] = cbcg_motoya(A, b, Min, Max, k, Size)
    aa = 2/(Max-Min)
    bb = -(Max+Min)/(Max-Min)
    x = zeros(b)
    r = b
    hist(1,:) = [0,(r' * r)]
    for i=1:(floor(2*Size/(k+1))+10)
        for j=1:k
            if j == 1 then
                S(:,1) = r
            elseif j == 2 then
                S(:,2) = aa * A * r + bb * r
            else
                S(:,j) = 2 * aa * A * S(:,j-1) + 2 * bb * S(:,j-1) - S(:,j-2)
            end
        end
        if i == 1 then
            Q = S
        else
            B = lsq((Q' * A * Q), (Q' * A * S), 1e-15)
            Q = S - Q * B
        end
        a = lsq(Q' * A * Q, Q' * r, 1e-15)
        r = r - A * Q * a
        x = x + Q * a
        hist((i+1),:) = [i*k,(r' * r)]
        if (norm(r) / norm(b) < 1e-10) then
            break
        end
    end
endfunction

// cbcg_1, compare to cbcg_motoya,
// not to use lsq, but try to use inv to find inverse directly
// the pitfall is that sometimes, when Qtranspose_A_Q is singular,
// the inv function fail
function [x, hist] = cbcg_1(A, b, Min, Max, k, Size)
    aa = 2/(Max-Min)
    bb = -(Max+Min)/(Max-Min)
    x = zeros(b)
    r = b
    hist(1,:) = [0,(r' * r)]
    for i=1:(floor(2*Size/(k+1))+10)
        for j=1:k
            if j == 1 then
                S(:,1) = r
            elseif j == 2 then
                S(:,2) = aa * A * r + bb * r
            else
                S(:,j) = 2 * aa * A * S(:,j-1) + 2 * bb * S(:,j-1) - S(:,j-2)
            end
        end
        if i == 1 then
            Q = S
        else
            B = lsq((Q' * A * Q), (Q' * A * S), 1e-15)
            Q = S - Q * B
        end
        a = lsq(Q' * A * Q, Q' * r, 1e-15)
        r = r - A * Q * a
        x = x + Q * a
        hist((i+1),:) = [i*k,(r' * r)]
        if (norm(r) / norm(b) < 1e-10) then
            break
        end
    end
endfunction

function cbcg_main(filename, A, Min, Max, Size, b, n, k)
    [x,hist]=cbcg(A,A*b(1:Size,1),Min,Max,k, Size);
    for i = 2:n
        [x,tmp]=cbcg(A,A*b(1:Size,i),Min,Max,k, Size); 
	if hist($,1) < tmp($,1) then
            hist = tmp
	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
            hist = tmp
        end
    end
    str=sprintf("cbcg_%s_k=%d_sample=%d.txt",filename,k,n)
    fprintfMat(str,hist,"%10.7e")
    disp([k,hist($,1)])
endfunction

stacksize('max')
//cd D:\WorkSpace\SparseLinAlgScilab\cg
//cd /home/sc2012/SparseLinAlgScilab-gh/cg
cd /home/scl/SparseLinAlgScilab-gh/cg
exec('Matrix.sci');

epsilon = 1e-20;
max_iters = 10;
rhs_m = 1;
num_samples = 1;
//b=fscanfMat("/home/skkk/ExperimentsRandom/Random");
//b=rand(5000,rhs_m * num_samples);
b=zeros(5000,rhs_m * num_samples);

//filename="/home/sc2012/MStore/SPD/bcsstk26.mtx";
//filename="/home/sc2012/MStore/SPD/shallow_water2.mtx";
//filename="/home/sc2012/MStore/SPD/nasasrb.mtx";

//filename="/home/scl/MStore/SPD/bcsstk26.mtx";
//filename="/home/scl/MStore/SPD/sts4098.mtx";
filename="/home/scl/MStore/SPD/crystm01.mtx";

//[A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename); // the returned matrix is preconditioed
[A,num_rows,num_cols,entries] = Matrix_nonprecondtioned(filename); // the returned matrix is nonpreconditioed

cbcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);


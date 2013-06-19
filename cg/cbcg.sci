function B = normalize(A)
    [m, n] = size(A)
    for i=1:n
        d = norm(A(:,i))
        B(:,i) = A(:,i) ./ d
    end
endfunction

// 2013/01/09 cbcg version of motoya-san
function [x, hist] = cbcg_motoya_1(A, b, Min, Max, k)
    aa = 2/(Max-Min)
    bb = -(Max+Min)/(Max-Min)
    x = zeros(b)
    r = b
    hist(1,:) = [0,(r' * r)]
    for i=1:1000
        for j=1:k
            if j == 1 then
                S(:,1) = r
            elseif j == 2 then
                S(:,2) = aa * A * r + bb * r
            else
                S(:,j) = 2 * aa * A * S(:,j-1) + 2 * bb * S(:,j-1) - S(:,j-2)
            end
        end
        S = normalize(S)
        if i == 1 then
            Q = S
        else
            B = lsq(Q' * A * Q, Q' * A * S, 1e-15)
            Q = S - Q * B
        end
        Q = normalize(Q)
        a = lsq(Q' * A * Q, Q' * r, 1e-15)
        r = r - A * Q * a
        x = x + Q * a
        hist((i+1),:) = [i*k,(r' * r)]
        if (norm(r) / norm(b) < 1e-10) then
            break
        end
    end
endfunction

function [x, hist] = cbcg_motoya_2(A, b, Min, Max, k, Size)
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
function [x, hist] = cbcg_1(A, b,s_k, max_iters, epsilon)
    exec('gerschgorin.sci');
    b_norm=norm(b);
    [num_row_A, num_col_A] = size(A);
    [eigenValue_min,eigenValue_max]=gerschgorin(A, num_row_A);
    s_alpha = 2.0 / (eigenValue_max - eigenValue_min);
    s_beta  = -(eigenValue_max + eigenValue_min) / (eigenValue_max - eigenValue_min);
    
    x = ones(b);
    r = b - A * x;
    
    hist(1,:) = [0,(r' * r)];
    disp ("iteration starts........");
    
    for i=1:max_iters
        for j=1:s_k
            if j == 1 then
                S(:,1) = r;
            elseif j == 2 then
                S(:,2) = s_alpha * A * r + s_beta * r;
            else
                S(:,j) = 2 * s_alpha * A * S(:,j-1) + 2 * s_beta * S(:,j-1) - S(:,j-2);
            end
        end
        
        if i == 1 then
            Q = S;
        else
            Qtranspose_A_S  = A_Q' * S;
            B = Qtranspose_A_Q_inverse * Qtranspose_A_S;
            Q = S - Q * B
        end
        
        A_Q = A * Q;
        Qtranspose_A_Q = Q' * A_Q;
        Qtranspose_A_Q_inverse = inv(Qtranspose_A_Q);
        
        v_a = Qtranspose_A_Q_inverse * Q' * r;
        x = x + Q * v_a;
        r = r - A_Q * v_a;
        
        s_r_norm = norm(r);
//        disp("iteration",i);
//        disp (s_r_norm);
        hist((i+1),:) = [i*s_k,s_r_norm];
        disp(hist((i+1),:));
        if (s_r_norm / b_norm < epsilon) then
            disp ("converged ... ... ...");
            break;
        end
//        if (s_r_norm < epsilon) then
//            disp ("converged ... ... ...");
//            break;
//        end
    end
endfunction

function cbcg_main(filename, A, b, s_k, max_iters, epsilon)
    [num_rows_b, num_cols_b] = size(b);
    [num_rows_A, num_cols_A] = size(A);
    num_samples = num_cols_b;

    [x,hist]=cbcg_1(A, b(1:num_rows_A,1),s_k, max_iters, epsilon);
    
//    for i = 2:num_samples
//        [x,tmp]=cbcg(A,A*b(1:Size,i),Min,Max,k, Size); 
//	if hist($,1) < tmp($,1) then
//            hist = tmp
//	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
//            hist = tmp
//        end
//    end
//    str=sprintf("cbcg_%s_k=%d_sample=%d.txt",filename,k,n)
//    fprintfMat(str,hist,"%10.7e")
//    disp([k,hist($,1)])
endfunction

stacksize('max')
//cd D:\WorkSpace\SparseLinAlgScilab\cg
//cd /home/sc2012/SparseLinAlgScilab-gh/cg
//cd /home/scl/SparseLinAlgScilab-gh/cg
cd C:\Users\sc2012\Documents\GitHub\SparseLinAlgScilab\cg
exec('Matrix.sci');

epsilon = 1e-15;
max_iters = 150;
s_k = 10; //2, 4, 10
num_samples = 1;
//b=fscanfMat("/home/skkk/ExperimentsRandom/Random");
b=rand(5000, num_samples);
//b=zeros(5000,num_samples);

//filename="/home/sc2012/MStore/SPD/bcsstk26.mtx";
//filename="/home/sc2012/MStore/SPD/shallow_water2.mtx";
//filename="/home/sc2012/MStore/SPD/nasasrb.mtx";
//filename="/home/sc2012/MStore/SPD/crystm01.mtx";

//filename="/home/scl/MStore/SPD/bcsstk26.mtx";
//filename="/home/scl/MStore/SPD/sts4098.mtx";
//filename="/home/scl/MStore/SPD/crystm01.mtx";

//
filename="C:\MStore\SPD\crystm01.mtx";
//filename="C:\MStore\SPD\ex9.mtx";

[A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename); // the returned matrix is preconditioed
//[A,num_rows,num_cols,entries] = Matrix_nonprecondtioned(filename); // the returned matrix is nonpreconditioed


cbcg_main(filename, A, b, s_k, max_iters, epsilon);


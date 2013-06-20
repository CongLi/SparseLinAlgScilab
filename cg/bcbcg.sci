// block chebyshev basis conjugate gradient method implementations
//
function [converge_label, max_difference, min_difference] =  Converge_Checking (R , B , epsilon)
    converge_label = %t;
    max_difference = -%inf;
    min_difference = %inf;
    
    assert_checktrue (size(R) == size(B));
    
    [num_rows, num_cols] = size(R);
    assert_checktrue(num_cols>0);
//    disp (num_cols);
    if num_cols>1 then
        for i = 1: num_cols
            norm_list_R(i) = norm ( R(:,i) );
//          disp (norm_list_R(i));
            norm_list_B(i) = norm ( B(:,i) );
//          disp (norm_list_B(i));
            check_sign = norm_list_R(i) < (norm_list_B(i) * epsilon);
            if ~check_sign then
                converge_label = %f;
                if ( norm_list_R(i) / norm_list_B(i)) > max_difference then
                    max_difference = norm_list_R(i) / norm_list_B(i);
                end
                if ( norm_list_R(i) / norm_list_B(i)) < min_difference then
                    min_difference = norm_list_R(i) / norm_list_B(i);
                end
            end
        end
    else
        r_norm = norm(R);
        check_sign = r_norm / norm(B) < epsilon;
        if ~check_sign then
            converge_label = %f;
            max_difference = r_norm / norm(B);
        end
    end

endfunction
//
function [X, hist] = bcbcg(A, B, rhs_m, s_k, max_iters,epsilon)
//    disp ("size of b_rhs_m is :", size(b_rhs_m));
    exec('gerschgorin.sci');
    [num_row_A, num_col_A] = size(A);
    [eigenValue_min,eigenValue_max]=gerschgorin(A, num_row_A);
    s_alpha = 2.0 / (eigenValue_max - eigenValue_min);
    s_beta  = -(eigenValue_max + eigenValue_min) / (eigenValue_max - eigenValue_min);
    
    X = ones(B);
    R = B - A * X;
    S = ones(num_row_A, s_k * rhs_m);

    [converge_label , max_difference, min_difference] = Converge_Checking (R,B,epsilon);
    hist(1,:) = [0, max_difference, min_difference];
    
    disp ("iteration starts ... ... ... ...");
    for i=1:max_iters
        // generate S of size (num_row_A, s_k * rhs_m)
        for j=1:s_k
            if j == 1 then
                S(:,1:rhs_m) = R;
            elseif j == 2 then
                S(:, (rhs_m + 1) : (2 * rhs_m)) = s_alpha * A * R + s_beta * R;
            else
                S(:,((j-1)*rhs_m + 1):(j*rhs_m)) = 2 * s_alpha * A * S(:,( (j-2) * rhs_m +1 ):( (j-1) * rhs_m) ) + 2 * s_beta * S(:,( (j-2) * rhs_m +1 ):( (j-1) * rhs_m) ) - S(:,( (j-3) * rhs_m +1 ):( (j-2) * rhs_m) );
            end
        end

        if i == 1 then
            // Q size (num_row_A, s_k * s_m)
            Q = S;
        else
            Qtranspose_A_S  = A_Q' * S;
            BETA = Qtranspose_A_Q_inverse * Qtranspose_A_S;
            Q = S - Q * BETA;
        end
        
        A_Q = A * Q;
        // Qtranspose_A_Q size (s_k * s_m, s_k * s_m)
        Qtranspose_A_Q = Q' * A_Q;
        // Qtranspose_A_Q_inverse size (s_k * s_m, s_k * s_m)
        Qtranspose_A_Q_inverse = inv(Qtranspose_A_Q);
        // LAMBDA size(s_k * s_m, s_m)
        LAMBDA = Qtranspose_A_Q_inverse * Q' * R;
        // 
        X = X + Q * LAMBDA; 
        R = R - A_Q * LAMBDA;
        
        [converge_label , max_difference, min_difference] = Converge_Checking (R,B,epsilon);
        hist(i+1,:) = [i*s_k,max_difference, min_difference];
        disp (hist(i+1,:));
        if converge_label then
            disp ("converged ...........");
            break;
        end
    end
endfunction

function [x,hist_bcbcg]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon)
    [num_rows, num_cols] = size(A);
//    disp (size(b(1:num_rows,1:rhs_m)));
    [x,hist_bcbcg]=bcbcg( A , b(1:num_rows,1:rhs_m) , rhs_m, s_k, max_iters,epsilon);
    for i = 2:num_samples
        [x,tmp]=bcbcg( A , b(1:num_rows, ((i-1)*rhs_m + 1):(i* rhs_m)), max_iters,epsilon); 
//      if hist($,1) < tmp($,1) then
//            hist = tmp
//	     elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
//              hist = tmp
//      end
    end
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
max_iters = 350;
rhs_m = 3;
s_k =20;
num_samples = 1;
//b=fscanfMat("/home/skkk/ExperimentsRandom/Random");
b=rand(5000,rhs_m * num_samples);
//b=zeros(5000,rhs_m * num_samples);

//filename="/home/sc2012/MStore/SPD/bcsstk26.mtx";
//filename="/home/sc2012/MStore/SPD/shallow_water2.mtx";
//filename="/home/sc2012/MStore/SPD/nasasrb.mtx";
//filename="/home/sc2012/MStore/SPD/crystm01.mtx";

//filename="/home/scl/MStore/SPD/bcsstk26.mtx";
//filename="/home/scl/MStore/SPD/sts4098.mtx";
//filename="/home/scl/MStore/SPD/crystm01.mtx";

//
//filename="C:\MStore\SPD\crystm01.mtx";
filename="C:\MStore\SPD\bcsstk16.mtx"

[A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename); // the returned matrix is preconditioed
//[A,num_rows,num_cols,entries] = Matrix_nonprecondtioned(filename); // the returned matrix is nonpreconditioed

//[x,hist_bcbcg_m3k1]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//[x,hist_bcbcg_m3k3]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//[x,hist_bcbcg_m3k5]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//[x,hist_bcbcg_m3k8]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//[x,hist_bcbcg_m3k10]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//[x,hist_bcbcg_m3k20]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);


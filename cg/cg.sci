function [x, hist] = cg(A, b, Size, max_iters,epsilon)
    x = ones(b);
    b_norm = norm(b);
    r = b - A * x;
//    b_norm = norm(b);
    p = r;
    rtranspose_r = (r' * r);
    // hist[#iter, algorihm generated r: ||r||/||b||, r=b-Ax: ||r||/||b||]
    hist(1,:) = [0, sqrt(rtranspose_r)/b_norm,sqrt(rtranspose_r)/b_norm];
    A_p = A * p;
    ptranspose_A_p = p' * A_p;
    ptranspose_A_p_inverse = inv(ptranspose_A_p);
    
    for i=1:max_iters
        s_alpha = rtranspose_r * ptranspose_A_p_inverse;
        x = x + s_alpha * p;
        r = r - s_alpha * A_p;
        s_beta = inv(rtranspose_r);
        rtranspose_r = r' * r ;
        r_norm = sqrt(rtranspose_r);
        b_Ax_cg= b - A * x;
        hist((i+1),:) = [i,r_norm / b_norm, norm(b_Ax_cg)/b_norm];
        disp ( hist((i+1),:));
        
        if (r_norm / b_norm < epsilon) then
            break;
        end
//        if (hist((i+1),2) < epsilon) then
//            break;
//        end
        
        s_beta = s_beta * rtranspose_r;
        p = r + s_beta * p;
        
        A_p = A * p;
        ptranspose_A_p = p' * A_p;
        ptranspose_A_p_inverse = inv(ptranspose_A_p);
        
    end
endfunction

function [x,hist_cg]=cg_main (filename, A, num_rows, b, num_samples, max_iters, epsilon)
    
    [x,hist_cg]=cg(A,b(1:num_rows,1), num_rows, max_iters,epsilon);
//    for i = 2:num_samples
//        [x,tmp]=cg(A,b(1:num_rows,i), max_iters,epsilon);
//	if hist($,1) < tmp($,1) then
//            hist = tmp
//	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
//            hist = tmp
//        end
//    end
//    str=sprintf("cg_%s_sample=%d.txt",filename,n)
//    fprintfMat(str,hist,"%10.7e")
//    disp(hist($,1))
endfunction
//
//stacksize('max')
////cd D:\WorkSpace\SparseLinAlgScilab\cg
////cd /home/sc2012/SparseLinAlgScilab-gh/cg
////cd /home/scl/SparseLinAlgScilab-gh/cg
//cd C:\Users\sc2012\Documents\GitHub\SparseLinAlgScilab\cg
//exec('Matrix.sci');
//
//epsilon = 1e-15;
//max_iters = 350;
//num_samples = 1;
////b=fscanfMat("/home/skkk/ExperimentsRandom/Random");
////b=zeros(5000,num_samples);
//b=rand(5000,num_samples);
//
////filename="/home/sc2012/MStore/SPD/bcsstk26.mtx";
////filename="/home/sc2012/MStore/SPD/shallow_water2.mtx";
////filename="/home/sc2012/MStore/SPD/nasasrb.mtx";
////filename="/home/sc2012/MStore/SPD/crystm01.mtx";
//
////filename="/home/scl/MStore/SPD/bcsstk26.mtx";
////filename="/home/scl/MStore/SPD/sts4098.mtx";
////filename="/home/scl/MStore/SPD/shallow_water2.mtx";
////filename="/home/scl/MStore/SPD/crystm01.mtx";
//
//// divergence
////filename="C:\MStore\SPD\msc01050.mtx";
////filename="C:\MStore\SPD\bcsstk21.mtx";
////filename="C:\MStore\SPD\1138_bus.mtx";
////filename="C:\MStore\SPD\nasa4704.mtx";
////filename="C:\MStore\SPD\ex9.mtx";
////
////filename="C:\MStore\SPD\bcsstk15.mtx"; //convergence slow more 1000 iterations
////filename="C:\MStore\SPD\1138_bus.mtx"; //convergence slow more 1000 iterations
////filename="C:\MStore\SPD\sts4098.mtx";  //convergence slow around 700 iterations
////
////filename="C:\MStore\SPD\crystm01.mtx";
////filename="C:\MStore\SPD\mesh2e1.mtx";
////filename="C:\MStore\SPD\bcsstm07.mtx";
////filename="C:\MStore\SPD\Chem97ZtZ.mtx";
//filename="C:\MStore\SPD\bcsstk16.mtx"
//
//
//[A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename); // the returned matrix is preconditioed
////[A,num_rows,num_cols,entries] = Matrix_nonprecondtioned(filename); // the returned matrix is nonpreconditioed
//
//[x,hist_cg]=cg_main(filename, A, num_rows , b, num_samples, max_iters,epsilon);


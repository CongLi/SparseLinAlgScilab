//20130624
stacksize('max');
cd C:\Users\sc2012\Documents\GitHub\SparseLinAlgScilab\cg
exec('Matrix.sci');
exec("cg.sci");
exec("bcg.sci");
exec("cbcg.sci");
exec("bcbcg.sci");
//
MatNameToSave="bcsstk16";
foldername="C:\Users\sc2012\Dropbox\Daily_Doc\20130624\ExpData\";
//
epsilon = 1e-15;
max_iters = 350;
num_samples = 1;
//
// for generating and saving b
//rowSizeMat = 5000;
//rhs_m = 30;
//num_samples = 1;
//b=rand(rowSizeMat,rhs_m * num_samples);
//str_Bmat=sprintf("Bmat_%s_row%d_col%d.txt",MatNameToSave,rowSizeMat,rhs_m);
//fprintfMat(foldername+str_Bmat,b,"%10.20e");
//
clear b;
rowSizeMat = 5000;
rhs_m = 30;
num_samples = 1;
str_Bmat=sprintf("Bmat_%s_row%d_col%d.txt",MatNameToSave,rowSizeMat,rhs_m);
b=fscanfMat(foldername+str_Bmat);
filename="C:\MStore\SPD\bcsstk16.mtx";
[A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename); // the returned matrix is preconditioed
//
////cg experiment
//[x,hist_cg]=cg_main(filename, A, num_rows , b, num_samples, max_iters,epsilon);
//////bcg experiment
//rhs_m = 1;
//[x,hist_bcg_m1]=bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
//rhs_m = 3;
//[x,hist_bcg_m3]=bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
//rhs_m = 5;
//[x,hist_bcg_m5]=bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
//rhs_m = 8;
//[x,hist_bcg_m8]=bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
//rhs_m = 10;
//[x,hist_bcg_m10]=bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
//rhs_m = 20;
//[x,hist_bcg_m20]=bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
////cbcg experiment
//s_k = 3
//[x,hist_cbcg_k3]=cbcg_main(filename, A, b, s_k, max_iters, epsilon);
//s_k = 5
//[x,hist_cbcg_k5]=cbcg_main(filename, A, b, s_k, max_iters, epsilon);
//s_k = 8
//[x,hist_cbcg_k8]=cbcg_main(filename, A, b, s_k, max_iters, epsilon);
//s_k = 10
//[x,hist_cbcg_k10]=cbcg_main(filename, A, b, s_k, max_iters, epsilon);
//s_k = 20
//[x,hist_cbcg_k20]=cbcg_main(filename, A, b, s_k, max_iters, epsilon);
////bcbcg
//rhs_m=3;
//s_k=1;
//[x,hist_bcbcg_m3k1]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//rhs_m=3;
//s_k=3;
//[x,hist_bcbcg_m3k3]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//rhs_m=3;
//s_k=5;
//[x,hist_bcbcg_m3k5]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//rhs_m=3;
//s_k=8;
//[x,hist_bcbcg_m3k8]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//rhs_m=3;
//s_k=10;
//[x,hist_bcbcg_m3k10]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//rhs_m=3;
//s_k=20;
//[x,hist_bcbcg_m3k20]=bcbcg_main(filename, A, b, rhs_m, s_k, num_samples, max_iters,epsilon);
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////WRITE part
//clear foldername;
//foldername="C:\Users\sc2012\Dropbox\Daily_Doc\20130624\ExpData\";
////////////////
//str_cg=sprintf("cg_%s.txt",MatNameToSave);
//////
//str_bcg_m1=sprintf("bcg_m1_%s.txt",MatNameToSave);
//str_bcg_m3=sprintf("bcg_m3_%s.txt",MatNameToSave);
//str_bcg_m5=sprintf("bcg_m5_%s.txt",MatNameToSave);
//str_bcg_m8=sprintf("bcg_m8_%s.txt",MatNameToSave);
//str_bcg_m10=sprintf("bcg_m10_%s.txt",MatNameToSave);
//str_bcg_m20=sprintf("bcg_m20_%s.txt",MatNameToSave);
//////
//str_cbcg_k3=sprintf("cbcg_k3_%s.txt",MatNameToSave);
//str_cbcg_k5=sprintf("cbcg_k5_%s.txt",MatNameToSave);
//str_cbcg_k8=sprintf("cbcg_k8_%s.txt",MatNameToSave);
//str_cbcg_k10=sprintf("cbcg_k10_%s.txt",MatNameToSave);
//str_cbcg_k20=sprintf("cbcg_k20_%s.txt",MatNameToSave);
//////
//str_bcbcg_m3k1=sprintf("bcbcg_m3k1_%s.txt",MatNameToSave);
//str_bcbcg_m3k3=sprintf("bcbcg_m3k3_%s.txt",MatNameToSave);
//str_bcbcg_m3k5=sprintf("bcbcg_m3k5_%s.txt",MatNameToSave);
//str_bcbcg_m3k8=sprintf("bcbcg_m3k8_%s.txt",MatNameToSave);
//str_bcbcg_m3k10=sprintf("bcbcg_m3k10_%s.txt",MatNameToSave);
//str_bcbcg_m3k20=sprintf("bcbcg_m3k20_%s.txt",MatNameToSave);
////////////////
//fprintfMat(foldername+str_cg,hist_cg,"%10.20e");
//////
//fprintfMat(foldername+str_bcg_m1,hist_bcg_m1,"%10.20e");
//fprintfMat(foldername+str_bcg_m3,hist_bcg_m3,"%10.20e");
//fprintfMat(foldername+str_bcg_m5,hist_bcg_m5,"%10.20e");
//fprintfMat(foldername+str_bcg_m8,hist_bcg_m8,"%10.20e");
//fprintfMat(foldername+str_bcg_m10,hist_bcg_m10,"%10.20e");
//fprintfMat(foldername+str_bcg_m20,hist_bcg_m20,"%10.20e");
//////
//fprintfMat(foldername+str_cbcg_k3,hist_cbcg_k3,"%10.20e");
//fprintfMat(foldername+str_cbcg_k5,hist_cbcg_k5,"%10.20e");
//fprintfMat(foldername+str_cbcg_k8,hist_cbcg_k8,"%10.20e");
//fprintfMat(foldername+str_cbcg_k10,hist_cbcg_k10,"%10.20e");
//fprintfMat(foldername+str_cbcg_k20,hist_cbcg_k20,"%10.20e");
//////
//fprintfMat(foldername+str_bcbcg_m3k1,hist_bcbcg_m3k1,"%10.20e");
//fprintfMat(foldername+str_bcbcg_m3k3,hist_bcbcg_m3k3,"%10.20e");
//fprintfMat(foldername+str_bcbcg_m3k5,hist_bcbcg_m3k5,"%10.20e");
//fprintfMat(foldername+str_bcbcg_m3k8,hist_bcbcg_m3k8,"%10.20e");
//fprintfMat(foldername+str_bcbcg_m3k10,hist_bcbcg_m3k10,"%10.20e");
//fprintfMat(foldername+str_bcbcg_m3k20,hist_bcbcg_m3k20,"%10.20e");

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////Reading Part
//hist_cg=fscanfMat(foldername+str_cg);
//hist_cbcg=fscanfMat(foldername+str_cbcg);
//hist_bcg=fscanfMat(foldername+str_bcg);
//hist_bcbcg_m3k1=fscanfMat(foldername+str_bcbcg_m3k1);
//hist_bcbcg_m3k5=fscanfMat(foldername+str_bcbcg_m3k5);
/////////////////
//hist_bcg_m1=fscanfMat(foldername+str_bcg_m1);
//hist_bcg_m3=fscanfMat(foldername+str_bcg_m3);
//hist_bcg_m5=fscanfMat(foldername+str_bcg_m5);
//hist_bcg_m8=fscanfMat(foldername+str_bcg_m8);
//hist_bcg_m10=fscanfMat(foldername+str_bcg_m10);
//hist_bcg_m20=fscanfMat(foldername+str_bcg_m20);
////////////////
//hist_cbcg_k3=fscanfMat(foldername+str_cbcg_k3);
//hist_cbcg_k5=fscanfMat(foldername+str_cbcg_k5);
//hist_cbcg_k8=fscanfMat(foldername+str_cbcg_k8);
//hist_cbcg_k10=fscanfMat(foldername+str_cbcg_k10);
//hist_cbcg_k20=fscanfMat(foldername+str_cbcg_k20);
///////////////
//hist_bcbcg_m3k3=fscanfMat(foldername+str_bcbcg_m3k3);
//hist_bcbcg_m3k5=fscanfMat(foldername+str_bcbcg_m3k5);
//hist_bcbcg_m3k8=fscanfMat(foldername+str_bcbcg_m3k8);
//hist_bcbcg_m3k10=fscanfMat(foldername+str_bcbcg_m3k10);
//hist_bcbcg_m3k20=fscanfMat(foldername+str_bcbcg_m3k20);
//////////////////////////////////////////////////////////////////////////////////////////
////// Drawing part Part
////////////////////////
////// pic 1
//////#inner product computed R
//plot2d("nl",hist_cg(:,1)*2,hist_cg(:,2),style=1);
//plot2d("nl",hist_cbcg_k5(:,1)*3,hist_cbcg_k5(:,2),style=2);
//plot2d("nl",hist_bcg_m3(:,1)*2,hist_bcg_m3(:,2),style=3);
//plot2d("nl",hist_bcbcg_m3k1(:,1)*3,hist_bcbcg_m3k1(:,2),style=4);
//plot2d("nl",hist_bcbcg_m3k5(:,1)*3,hist_bcbcg_m3k5(:,2),style=5);
////////#inner product b-Ax
//plot2d4("nl",hist_cg(:,1)*2,hist_cg(:,3),style=6);
//plot2d4("nl",hist_cbcg_k5(:,1)*3,hist_cbcg_k5(:,3),style=7);
//plot2d4("nl",hist_bcg_m3(:,1)*2,hist_bcg_m3(:,4),style=15); // style 8 cannot see
//plot2d4("nl",hist_bcbcg_m3k1(:,1)*3,hist_bcbcg_m3k1(:,4),style=9);
//plot2d4("nl",hist_bcbcg_m3k5(:,1)*3,hist_bcbcg_m3k5(:,4),style=10);
////////
//xtitle("Converging history","Iterations(#Inner Products)","||r||/||b||");
//set(gca(),"grid",[5 5]);
//legend(['cg';'cbcg,k=5';'bcg,m=3';'bcbcg,m=3,k=1';'bcbcg,m=3,k=5';'cg,r=b-Ax';'cbcg,k=5,r=b-Ax';'bcg,m=3,r=b-Ax';'bcbcg,m=3,k=1,r=b-Ax';'bcbcg,m=3,k=5,r=b-Ax']);
////
//////////////////////////////////////////////////////////
//pic2
////////
//plot2d("nl",hist_bcg_m1(:,1)*2,hist_bcg_m1(:,2),style=1);
//plot2d("nl",hist_bcg_m3(:,1)*2,hist_bcg_m3(:,2),style=2);
//plot2d("nl",hist_bcg_m5(:,1)*2,hist_bcg_m5(:,2),style=3);
//plot2d("nl",hist_bcg_m8(:,1)*2,hist_bcg_m8(:,2),style=4);
//plot2d("nl",hist_bcg_m10(:,1)*2,hist_bcg_m10(:,2),style=5);
//plot2d("nl",hist_bcg_m20(:,1)*2,hist_bcg_m20(:,2),style=6);
//plot2d4("nl",hist_bcg_m1(:,1)*2,hist_bcg_m1(:,4),style=7);
//plot2d4("nl",hist_bcg_m3(:,1)*2,hist_bcg_m3(:,4),style=15);  //style 8 cannot see
//plot2d4("nl",hist_bcg_m5(:,1)*2,hist_bcg_m5(:,4),style=9);
//plot2d4("nl",hist_bcg_m8(:,1)*2,hist_bcg_m8(:,4),style=10);
//plot2d4("nl",hist_bcg_m10(:,1)*2,hist_bcg_m10(:,4),style=11);
//plot2d4("nl",hist_bcg_m20(:,1)*2,hist_bcg_m20(:,4),style=12);
////
//xtitle("Converging history","Iterations(#Inner Products)","||r||/||b||");
//set(gca(),"grid",[5 5]);
//legend(['bcg,m=1';'bcg,m=3';'bcg,m=5';'bcg,m=8';'bcg,m=10';'bcg,m=20';'bcg,m=1,R=B-AX';'bcg,m=3,R=B-AX';'bcg,m=5,R=B-AX';'bcg,m=8,R=B-AX';'bcg,m=10,R=B-AX';'bcg,m=20,R=B-AX']);
//////
//pic3
//// computed r
//plot2d4("nl",hist_cbcg_k3(:,1)*3,hist_cbcg_k3(:,2),style=1);
//plot2d4("nl",hist_cbcg_k5(:,1)*3,hist_cbcg_k5(:,2),style=2);
//plot2d4("nl",hist_cbcg_k8(:,1)*3,hist_cbcg_k8(:,2),style=3);
//plot2d4("nl",hist_cbcg_k10(:,1)*3,hist_cbcg_k10(:,2),style=4);
//plot2d4("nl",hist_cbcg_k20(:,1)*3,hist_cbcg_k20(:,2),style=5);
//////
//plot2d("nl",hist_bcbcg_m3k3(:,1)*3,hist_bcbcg_m3k3(:,2),style=6);
//plot2d("nl",hist_bcbcg_m3k5(:,1)*3,hist_bcbcg_m3k5(:,2),style=7);
//plot2d("nl",hist_bcbcg_m3k8(:,1)*3,hist_bcbcg_m3k8(:,2),style=15); // style 8 cannot see
//plot2d("nl",hist_bcbcg_m3k10(:,1)*3,hist_bcbcg_m3k10(:,2),style=9);
//plot2d("nl",hist_bcbcg_m3k20(:,1)*3,hist_bcbcg_m3k20(:,2),style=10);
//xtitle("Converging history","Iterations(#Inner Products)","||r||/||b||");
//set(gca(),"grid",[5 5]);
//legend(['cbcg,k=3';'cbcg,k=5';'cbcg,k=8';'cbcg,k=10';'cbcg,k=20';'bcbcg,m=3,k=3';'bcbcg,m=3,k=5';'bcbcg,m=3,k=8';'bcbcg,m=3,k=10';'bcbcg,m=3,k=20']);
///////// R=B-AX
//plot2d4("nl",hist_cbcg_k3(:,1)*3,hist_cbcg_k3(:,3),style=11);
//plot2d4("nl",hist_cbcg_k5(:,1)*3,hist_cbcg_k5(:,3),style=12);
//plot2d4("nl",hist_cbcg_k8(:,1)*3,hist_cbcg_k8(:,3),style=13);
//plot2d4("nl",hist_cbcg_k10(:,1)*3,hist_cbcg_k10(:,3),style=14);
//plot2d4("nl",hist_cbcg_k20(:,1)*3,hist_cbcg_k20(:,3),style=25); // 15 is dominated
////
//plot2d("nl",hist_bcbcg_m3k3(:,1)*3,hist_bcbcg_m3k3(:,4),style=16);
//plot2d("nl",hist_bcbcg_m3k5(:,1)*3,hist_bcbcg_m3k5(:,4),style=17);
//plot2d("nl",hist_bcbcg_m3k8(:,1)*3,hist_bcbcg_m3k8(:,4),style=18); 
//plot2d("nl",hist_bcbcg_m3k10(:,1)*3,hist_bcbcg_m3k10(:,4),style=19);
//plot2d("nl",hist_bcbcg_m3k20(:,1)*3,hist_bcbcg_m3k20(:,4),style=20);
//xtitle("Converging history","Iterations(#Inner Products)","||r||/||b||");
//set(gca(),"grid",[5 5]);
//legend(['cbcg,k=3';'cbcg,k=5';'cbcg,k=8';'cbcg,k=10';'cbcg,k=20';'bcbcg,m=3,k=3';'bcbcg,m=3,k=5';'bcbcg,m=3,k=8';'bcbcg,m=3,k=10';'bcbcg,m=3,k=20']);
//legend(['cbcg,k=3';'cbcg,k=5';'cbcg,k=8';'cbcg,k=10';'cbcg,k=20';'bcbcg,m=3,k=3';'bcbcg,m=3,k=5';'bcbcg,m=3,k=8';'bcbcg,m=3,k=10';'bcbcg,m=3,k=20';'cbcg,k=3,R=B-AX';'cbcg,k=5,R=B-AX';'cbcg,k=8,R=B-AX';'cbcg,k=10,R=B-AX';'cbcg,k=20,R=B-AX';'bcbcg,m=3,k=3,R=B-AX';'bcbcg,m=3,k=5,R=B-AX';'bcbcg,m=3,k=8,R=B-AX';'bcbcg,m=3,k=10,R=B-AX';'bcbcg,m=3,k=20,R=B-AX']);


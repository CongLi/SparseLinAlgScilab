//20130624
//    str=sprintf("cg_%s_sample=%d.txt",filename,n)
//    fprintfMat(str,hist,"%10.7e")

//MatNameToSave="bcsstk16";
//foldername="C:\Users\sc2012\Dropbox\Daily_Doc\20130620\";
////WRITE part
//str_cg=sprintf("cg_%s.txt",MatNameToSave);
////
//str_bcg=sprintf("bcg_%s.txt",MatNameToSave);
//str_bcg_m1=sprintf("bcg_m1_%s.txt",MatNameToSave);
//str_bcg_m3=sprintf("bcg_m3_%s.txt",MatNameToSave);
//str_bcg_m5=sprintf("bcg_m5_%s.txt",MatNameToSave);
//str_bcg_m8=sprintf("bcg_m8_%s.txt",MatNameToSave);
//str_bcg_m10=sprintf("bcg_m10_%s.txt",MatNameToSave);
//str_bcg_m20=sprintf("bcg_m20_%s.txt",MatNameToSave);
////
//str_cbcg=sprintf("cbcg_%s.txt",MatNameToSave);
//str_cbcg_k3=sprintf("cbcg_k3_%s.txt",MatNameToSave);
//str_cbcg_k5=sprintf("cbcg_k5_%s.txt",MatNameToSave);
//str_cbcg_k8=sprintf("cbcg_k8_%s.txt",MatNameToSave);
//str_cbcg_k10=sprintf("cbcg_k10_%s.txt",MatNameToSave);
//str_cbcg_k20=sprintf("cbcg_k20_%s.txt",MatNameToSave);
////
//str_bcbcg_m3k1=sprintf("bcbcg_m3k1_%s.txt",MatNameToSave);
//str_bcbcg_m3k3=sprintf("bcbcg_m3k3_%s.txt",MatNameToSave);
//str_bcbcg_m3k5=sprintf("bcbcg_m3k5_%s.txt",MatNameToSave);
//str_bcbcg_m3k8=sprintf("bcbcg_m3k8_%s.txt",MatNameToSave);
//str_bcbcg_m3k10=sprintf("bcbcg_m3k10_%s.txt",MatNameToSave);
//str_bcbcg_m3k20=sprintf("bcbcg_m3k20_%s.txt",MatNameToSave);
//
//depending on the precision, %10.20 can make fscanfMat result almost the 
//same as original data
//fprintfMat(str_cg,hist_cg,"%10.20e");
////
//fprintfMat(str_bcg,hist_bcg,"%10.20e");
//fprintfMat(str_bcg_m1,hist_bcg_m1,"%10.20e");
//fprintfMat(str_bcg_m3,hist_bcg_m3,"%10.20e");
//fprintfMat(str_bcg_m5,hist_bcg_m5,"%10.20e");
//fprintfMat(str_bcg_m8,hist_bcg_m8,"%10.20e");
//fprintfMat(str_bcg_m10,hist_bcg_m10,"%10.20e");
//fprintfMat(str_bcg_m20,hist_bcg_m20,"%10.20e");
////
//fprintfMat(str_cbcg,hist_cbcg,"%10.20e");
//fprintfMat(str_cbcg_k3,hist_cbcg_k3,"%10.20e");
//fprintfMat(str_cbcg_k5,hist_cbcg_k5,"%10.20e");
//fprintfMat(str_cbcg_k8,hist_cbcg_k8,"%10.20e");
//fprintfMat(str_cbcg_k10,hist_cbcg_k10,"%10.20e");
//fprintfMat(str_cbcg_k20,hist_cbcg_k20,"%10.20e");
////
//fprintfMat(str_bcbcg_m3k1,hist_bcbcg_m3k1,"%10.20e");
//fprintfMat(str_bcbcg_m3k3,hist_bcbcg_m3k3,"%10.20e");
//fprintfMat(str_bcbcg_m3k5,hist_bcbcg_m3k5,"%10.20e");
//fprintfMat(str_bcbcg_m3k8,hist_bcbcg_m3k8,"%10.20e");
//fprintfMat(str_bcbcg_m3k10,hist_bcbcg_m3k10,"%10.20e");
//fprintfMat(str_bcbcg_m3k20,hist_bcbcg_m3k20,"%10.20e");
//
// READ part
hist_cg=fscanfMat(foldername+str_cg);
hist_cbcg=fscanfMat(foldername+str_cbcg);
hist_bcg=fscanfMat(foldername+str_bcg);
hist_bcbcg_m3k1=fscanfMat(foldername+str_bcbcg_m3k1);
hist_bcbcg_m3k5=fscanfMat(foldername+str_bcbcg_m3k5);
//


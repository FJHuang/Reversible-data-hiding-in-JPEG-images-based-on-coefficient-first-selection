clear all;clc;
warning off;


[best_psnr1,inc1,runtime1]=optimal_1D('image/80/lena80.jpg',10000,80);
[best_psnr2,inc2,runtime2]=optimal_2D('image/80/lena80.jpg',10000,80);
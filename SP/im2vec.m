% IM2VEC  Reshape 2D image blocks into an array of column vectors
%
%    V=IM2VEC(IM,BLKSIZE[,PADSIZE])
%
%    [V,ROWS,COLS]=IM2VEC(IM,BLKSIZE[,PADSIZE])
%
%    IM is an image to be separated into non-overlapping blocks and
%    reshaped into an MxN array containing N blocks reshaped into Mx1
%    column vectors.  IM2VEC is designed to be the inverse of VEC2IM.
%
%    BLKSIZE is a scalar or 1x2 vector indicating the size of the blocks.
%
%    PADSIZE is a scalar or 1x2 vector indicating the amount of vertical
%    and horizontal space to be skipped between blocks in the image.
%    Default is [0 0].  If PADSIZE is a scalar, the same amount of space
%    is used for both directions.  PADSIZE must be non-negative (blocks
%    must be non-overlapping).
%
%    ROWS indicates the number of rows of blocks found in the image.
%    COLS indicates the number of columns of blocks found in the image.
%
%    See also VEC2IM.

% Phil Sallee 5/03

function [v,rows,cols]=im2vec(im,bsize,padsize)

bsize = bsize + [0 0];

if (nargin < 3)
  padsize = 0;
end
padsize = padsize + [0 0];
if (any(padsize<0))
  error('Pad size must not be negative.');
end
% 拓展矩阵
imsize = size(im);
y=bsize(1)+padsize(1);
x=bsize(2)+padsize(2);
rows = floor((imsize(1)+padsize(1))/y);%rows取小于imsize(1)+padsize(1))/y的值
cols = floor((imsize(2)+padsize(2))/x);%cols取小于imsize(2)+padsize(2))/x的值

t = zeros(y*rows,x*cols);%创建一个空矩阵
imy=y*rows-padsize(1);%计算矩阵的列
imx=x*cols-padsize(2);%计算矩阵的行
t(1:imy,1:imx)=im(1:imy,1:imx);%复制矩阵im
t = reshape(t,y,rows,x,cols);%将t矩阵重塑为y行rows列的x*cols个子矩阵
%将矩阵t按向量[1,3,2,4]的顺序变换，变为y行x列rows*cols的四维矩阵
%再将变换后的矩阵重塑为一个y行x列的子矩阵rows*cols个
t = reshape(permute(t,[1,3,2,4]),y,x,rows*cols);
v = t(1:bsize(1),1:bsize(2),1:rows*cols);
v = reshape(v,y*x,rows*cols);%将v重塑为y*x行，rows*cols列的矩阵

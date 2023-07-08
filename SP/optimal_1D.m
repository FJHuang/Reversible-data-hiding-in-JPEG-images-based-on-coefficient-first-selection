function [best_psnr,inc,runtime]=optimal_1D(filename,messLen,QF)
t1 = clock;
cover=filename(10:length(filename));
% QF = cover(isstrprop(cover,'digit'));
stegoJPEG = 'RDH.jpg';
jpegOrder = [ 1 9 2 3 10 17 25 18   11 4 5 12 19 26 33 41   34 27 20 13 6 7 14 21   28 35 42 49 57 50 43 36,...
             29 22 15 8 16 23 30 37   44 51 58 59 52 45 38 31   24 32 39 46 53 60 61 54    47 40 48 55 62 63 56 64];       
jpgObj = jpeg_read(filename);             %read the coefficients of JPEG image
coverJPEG = filename;
jpgCoef = jpgObj.coef_arrays{1};
[row,col] = size(jpgCoef);
jpgVecCoef = im2vec(jpgCoef,[8,8]);                                                                                                  
jpgEmbCoef = jpgVecCoef(jpegOrder(1:64),:);  % the 1:64 frequency bands are selected for embeding
S = jpgEmbCoef(2:64,:);

QT = jpgObj.quant_tables{1};
QT = im2vec(QT,[8,8]); 
QT = QT(jpegOrder(1:64),:);
QT = QT(2:64);

% 找出最相似的K个子块
k = 700;
Complex = zeros(4096,1);
for i = 1:4096
    tep0 = 0;
    for j = 1:63
        if S(j,i)~=0
            tep0 = tep0+QT(j);
        end
    end
    Complex(i) = tep0;
end


all_lambda = zeros(63,4096);
for i = 1:4096
    sim = abs(Complex(i)-Complex);
    [sort_sim,index_sim] = sort(sim,'ascend');
    select_index = index_sim(1:k);
    select_block = S(:,select_index);
    pz = sum(select_block(:,1:k)==0,2);
    pz = max(1,pz);
    pz = min(pz,k-1);
    lambda  = -2*log(1-pz/k);
    all_lambda(:,i) = lambda;
end


% E的计算
E = zeros(63,4096);
PC = zeros(63,4096);
for i = 1:63
    for j = 1:4096
        pz = CDF(all_lambda(i,j),0.5)-CDF(all_lambda(i,j),-0.5);
        pc = (CDF(all_lambda(i,j),-0.5)-CDF(all_lambda(i,j),-1.5))+(CDF(all_lambda(i,j),1.5)-CDF(all_lambda(i,j),0.5));
        ps = 1 - pz - pc;
        L = pc;
        D = (0.5*pc+ps)*(QT(i)^2);
        E(i,j) = L/D;
        PC(i,j) = pc;
    end
end

% 优先嵌入零系数多的块
Zerosnum = zeros(1,4096);
for i = 1:4095
    Zerosnum(1,i) = sum(S(:,i)==0);
end

[~,index] = sort(Zerosnum,'descend');
S = S(:,index);
E = E(:,index);

E_col = E(:);
[sort_E,~] = sort(E_col,'descend');
sort_E = sort(unique(sort_E),'descend');
[RowE,~]= size(sort_E);


%StartPoint 二分找起始点
left = 1;
right =RowE;

while(left < right)
    mid = floor((left+right)/2);
    e = sort_E(mid);
    SelectPos = (E>e);
    SelectCof = S.*SelectPos;
    EC =sum(SelectCof(:)==1)+sum(SelectCof(:)==-1);
    if EC >=messLen
        right = mid;
    else
        left = mid + 1;
    end
end

StartPoint = left;


S_copy = S;
%Embeddding
msg=randi([0 1],1,messLen);

best_psnr = 0;
inc = 0;
%引入早停机制减少时间复杂度
EarlyStopTime = 0;

for i = StartPoint:RowE
    num = 1;
    E_num = sort_E(i);
    S = S_copy;
    for j = 1:4096
        if num>messLen break; end
        for k  = 1:63
            if num>messLen break; end
            if S(k,j) == 1 && E(k,j) >= E_num
                S(k,j) = S(k,j)+msg(num);
                num = num+1;
            elseif S(k,j) == -1 && E(k,j) >= E_num
                S(k,j) = S(k,j)-msg(num);
                num = num+1;
            elseif  abs(S(k,j))>1 && E(k,j) >= E_num
                S(k,j) = S(k,j)+sign(S(k,j));
            end
        end
    end


    EarlyStopTime = EarlyStopTime + 1;

    S(:,index) = S;
    stegoVecCoef = jpgVecCoef;
    stegoVecCoef(jpegOrder(2:64),:) = S;  

    stegoCoef = vec2im(stegoVecCoef,0, 8, row/8, col/8);

    load(strcat('default_gray_jpeg_obj_', num2str(QF),'.mat'));
    stegoObj = default_gray_jpeg_obj;
    stegoObj.coef_arrays{1} = stegoCoef;
    stegoObj.optimize_coding = 1;

    jpeg_write(stegoObj,stegoJPEG);
    temp = dir(stegoJPEG);
    fileSize = temp.bytes;
    temp2 = dir(coverJPEG);
    infileSize = (temp.bytes-temp2.bytes)*8;
    cover = imread(coverJPEG);
    stego = imread(stegoJPEG);
    diff = double(cover) - double(stego);
    x = sum(sum(diff.*diff));
    PSNR= 10*log10(512*512*255*255/x);

    if(EarlyStopTime > 100)
        break;
    end


    if num>messLen && PSNR>best_psnr
        best_psnr = PSNR;
        inc = infileSize;
        EarlyStopTime = 0;
    end

end
t2 = clock;
runtime = etime(t2,t1);
end


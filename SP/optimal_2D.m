function [best_psnr,inc,runtime]=optimal_2D(filename,messLen,QF) 

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
% k_min = 50;k_max = 400;
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
%     k = k_min + floor((k_max-k_min)*(Complex(i)/63));
    select_index = index_sim(1:k);
    select_block = S(:,select_index);
    pz = sum(select_block(:,1:k)==0,2);
    pz = max(1,pz);
    pz = min(pz,k-1);
    lambda  = -2*log(1-pz/k);
    all_lambda(:,i) = lambda;
end


%系数失真排序
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

%零系数排序
sum0 = zeros(1,4096);
for i = 1:4096
    sum0(i) = sum(S(:,i)==0);
end  
[~,index] = sort(sum0,'descend');

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
    SelectCof = SelectCof(SelectCof~=0);
    num_NACP = floor(numel(SelectCof)/2);
    NACP = reshape(SelectCof(1:2*num_NACP), 2, [])';
    EC = 0;
    for i = 1:numel(NACP(:,1))
        if abs(NACP(i,1))==1 && abs(NACP(i,2))==1
            EC = EC + 1.5;
        elseif abs(NACP(i,1))==1 || abs(NACP(i,2))==1
            EC = EC + 1;
        elseif abs(NACP(i,1))==2 && abs(NACP(i,2))==2
            EC = EC + 1;
        end
    end
%     EC =sum(SelectCof(:)==1)+sum(SelectCof(:)==-1);
    if EC >=messLen
        right = mid;
    else
        left = mid + 1;
    end
end

StartPoint = left;


S_copy = S;
%Embeddding
msg=randi([0 1],1,messLen+1);

best_psnr = 0;
inc = 0;
%引入早退机制减少时间复杂度
EarlyStopTime = 0;

TK = 1024;


for k = StartPoint:RowE
    num = 1;
    E_num = sort_E(k);
    S = S_copy;
    x = 0; y = 0;
    for i = 1:4096
        if num>messLen break; end
        for j = 1:63
            if S(j,i) ~=0 && abs(S(j,i))<TK && E(j,i) >= E_num 
                if x == 0
                    x = S(j,i);
                    x_row = j;
                    x_col = i;
                elseif y == 0
                    y = S(j,i);
                    y_row = j;
                    y_col = i;
                end
            elseif S(j,i) ~=0 && abs(S(j,i))>=TK && E(j,i) >= E_num 
                    S(j,i)  = S(j,i) + sign(S(j,i));
            end
            if x ~=0 && y~=0
                if num>messLen break; end
                % 四个象限嵌入
                           if abs(x) == 1 && abs(y) == 1
                               if msg(num) == 1 && msg(num+1) ==0
                                   S(x_row,x_col) = x+sign(x); num = num+2;
                               elseif msg(num) == 1 && msg(num+1) ==1
                                   S(y_row,y_col) = y+sign(y); num = num+2;
                               else
                                   num = num+1;
                               end
                           elseif abs(x) == 1 && abs(y) >=2
                               if msg(num) == 0
                                   S(y_row,y_col) = y+sign(y);num = num+1;
                               else
                                   S(x_row,x_col) = x+sign(x); 
                                   S(y_row,y_col) = y+sign(y);num = num+1;
                               end
                               
                           elseif  abs(y) ==1 && abs(x) >=2
                               if msg(num) == 0
                                   S(x_row,x_col) = x+sign(x);num = num+1;
                               else
                                   S(x_row,x_col) = x+sign(x); 
                                   S(y_row,y_col) = y+sign(y);num = num+1;
                               end
                               
                           elseif abs(x) == 2 && abs(y) == 2
                               if msg(num) == 0
                                   num = num+1;
                               else
                                   S(x_row,x_col) = x+sign(x);
                                   S(y_row,y_col) = y+sign(y);num = num+1;
                               end
                            elseif abs(x) > 2 || abs(y) > 2
                                   S(x_row,x_col) = x+sign(x);
                                   S(y_row,y_col) = y+sign(y);                      
                           end
                   x = 0; y =0;                   
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

    if(EarlyStopTime > 100  && best_psnr ~=0)
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

if best_psnr == 0
    runtime = 0;
end
end


function [ ] = shibiechengxu3( )
% clear;
% close all;
%运行速度：3.272s，垂直倾斜改为微分
global im
global erzhi_yuzhi
global shibiefangfa
global tezheng
global RegCode
%from kevin
erzhi_yuzhi
shibiefangfa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_s=14;    %strel参数
rectangle_x=20;     %滤波矩形框大小
rectangle_y=100;
rectangle_width=[98,500];%连通区域（车牌）允许大小
rectangle_hight=[25,200];
pr_rato=[2,10];%连通区域（车牌）允许长宽比
supplement_x=5;%从连通区到车牌矩形提取的补充（能否自适应）
supplement_y=15;
%diff_ydata=2;%车牌旋转校正，数据修正，使预拟合曲线光滑，防止跳变（像素点）
size_x=40;%子图标准大小x
size_y=20;%子图标准大小y
shot_height=1080;
shot_width=1920;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=waitbar(0,'请稍等，正在开始识别文件...');

%=================图像读入与预处理=============================
rgb_image=im;
% !figure,imshow(rgb_image),title('原彩色图像');%显示原图像
gray_image=rgb2gray(rgb_image);%原图转灰度图
% !figure,imshow(gray_image),title('原灰度图');%显示灰度图

%获得增强图（去除背景）
s=strel('disk',parameter_s);%!!!!!!!strel参数修改
Bgray_image=imopen(gray_image,s);%开运算操作
% !figure,imshow(Bgray_image),title('背景图像')%显示背景图
Egray_image=imsubtract(gray_image,Bgray_image);%两图相减得增强图
% !figure,imshow(Egray_image),title('增强图像')%显示增强图
waitbar(1/20,h,'请等待，识别中...');

%获得二值化阈值
image_max=double(max(max(Egray_image)));%增强图最大值并输出双精度型
image_min=double(min(min(Egray_image)));%增强图最小值并输出双精度型
level1=(image_max-(image_max-image_min)/3)/255;%最佳阈值修改%获得最佳阈值
level2 = graythresh(gray_image);%OSTU计算阈值
switch erzhi_yuzhi
    case 1
        level=level1;
    case 2
        level=level2;
    case 3
        level=(level1+level2)/2;!!!取平均值
end
bw2=im2bw(rgb_image);%转换图像为二进制图像
bw2=double(bw2);%转换为双精度型
% !figure,imshow(bw2),title('二值图');
waitbar(2/20,h,'请等待，识别中...');

%对二值图像作开闭操作进行滤波(腐蚀膨胀)
grd=edge(bw2,'canny');%用canny算子识别强度图像中的边界
% !figure,imshow(grd);title('图像边缘提取');%输出图像边缘
bg1=imclose(grd,strel('rectangle',[rectangle_x,rectangle_y]));%!!!!!!开闭运算修改%取矩形框的闭运算
% !figure,imshow(bg1);title('图像闭运算[5,19]');%输出闭运算的图像
bg3=imopen(bg1,strel('rectangle',[rectangle_x,rectangle_y]));%取矩形框的开运算
%!figure,imshow(bg3);title('图像开运算[5,19]');%输出开运算的图像
bg2=imopen(bg3,strel('rectangle',[rectangle_y,1]));%取矩形框的开运算
%!figure,imshow(bg2);title('图像开运算[19,1]');%输出开运算的图像
waitbar(3/20,h,'请等待，识别中...');

%对二值图像进行区域提取，并计算区域特征参数。进行区域特征参数比较，提取车牌区域
[L,num] = bwlabel(bg2,8);%标注二进制图像中已连接的部分（L为连通区域矩阵，num为连通区域数目）
Feastats = regionprops(L,'basic');%计算图像区域的系列特征尺寸（'Area'是标量，计算出在图像各个区域中像素总个数。'BoundingBox'是1行ndims(L)*2列的向量，即包含相应区域的最小矩形。BoundingBox 形式为 [ul_corner width]，这里 ul_corner 以 [x y z ...] 的坐标形式给出边界盒子的左上角、boxwidth 以 [x_width y_width ...] 形式指出边界盒子沿着每个维数方向的长度。）
                                  %'Centroid'是1行ndims(L)列的向量，给出每个区域的质心（重心）。
% Area=[Feastats.Area];%区域面积
BoundingBox=[Feastats.BoundingBox];%[x y width height]车牌的框架大小
RGB_image2= label2rgb(L, 'spring', 'k', 'noshuffle'); %标志图像向RGB图像转换
file_name1=strcat('..\车牌标准图\车牌标记图','.jpg');%存储标准子图
imwrite(RGB_image2,file_name1,'jpg') %存储标准子图
% !figure,imshow(RGB_image2);title('图像彩色标记');%输出框架的彩色图像
waitbar(4/20,h,'请等待，识别中...');

%计算连通区域是否符合车牌大小
lx=0;
Getok=zeros(1,num);
for l=1:num
    width=BoundingBox((l-1)*4+3);%框架宽度的计算
    hight=BoundingBox((l-1)*4+4);%框架高度的计算
    if (width>rectangle_width(1) && width<rectangle_width(2) && hight>rectangle_hight(1) && hight<rectangle_hight(2))%框架的宽度和高度的范围
        lx=lx+1;
        Getok(lx)=l;
    end
end
%计算筛选后的连通区域的长宽比！！！！！！！！可考虑仅有1个连通区域减少运算时间
for k= 1:lx
    l=Getok(k);    
    startcol=BoundingBox((l-1)*4+1)-supplement_x;%连通区域左上角坐标x
    startrow=BoundingBox((l-1)*4+2)-supplement_y;%连通区域左上角坐标y
    if startcol<0||startrow<0
        continue
    end
    width=BoundingBox((l-1)*4+3)+supplement_x;%车牌宽
    hight=BoundingBox((l-1)*4+4)+supplement_y;%车牌高
    endcol=startcol+width;
    endrow=startrow+hight;
    if endcol>shot_width||endrow>shot_height
        continue
    end
    rato=width/hight;%计算车牌长宽比
    if rato>pr_rato(1) && rato<pr_rato(2)   
        break;
    end
end
sbw1=bw2(startrow:startrow+hight,startcol:startcol+width+5); %获取车牌二值子图
subcol1=gray_image(startrow:startrow+hight,startcol:startcol+width+5);%获取车牌灰度子图
%grd1=grd(startrow:startrow+hight,startcol:startcol+width+5);%获取车牌边界子图
% !figure,subplot(2,1,1),imshow(subcol1);title('车牌灰度子图');%显示灰度图像
% !subplot(2,1,2),imshow(sbw1);title('车牌二值子图');%显示车牌的二值图
%===============================预处理结束=======================================
waitbar(5/20,h,'请等待，识别中...');

%==============================车牌旋转校正======================================

%==================水平校正旋转=================

binaryImage_horizontal = edge(subcol1,'sobel','horizontal');
% !subplot(4,1,4),imshow(binaryImage_horizontal),title('水平边界图')
theta=0:179;
r=radon(binaryImage_horizontal,theta);
% [m,n]=size(r);%m=角度，n=位置
[~,temp] = find(r>=max(max(r)));%temp记录了倾斜角，最大的倾斜角
angle=90-temp;
subcol = imrotate(subcol1,angle,'bilinear','crop'); %旋转灰度子图
sbw = imrotate(sbw1,angle,'bilinear','crop');%旋转二值子图
% grd3 = edge(subcol,'canny');%重新检测子图边界 已证明相对旋转更好
% grd4 = bwmorph(grd3,'thicken'); %通过在目标外部增加像素加厚目标直到这样做最终使先前未连接目标成为8连通域。得到强化边界
[hight,width]=size(sbw);
%!figure,subplot(4,1,1),imshow(subcol);title('车牌灰度子图(水平校正旋转后)');%输出车牌旋转后的灰度图像标题显示车牌灰度子图
% !subplot(4,1,2),imshow(sbw);title('车牌二值子图（水平校正旋转后）')%输出车牌旋转后的灰度图像
% !subplot(4,1,3),imshow(grd3);title('车牌边界子图（水平校正旋转后）')%输出车牌旋转后的边界子图
% !subplot(4,1,4),imshow(grd4);%输出车牌旋转后的强化边界子图
% !title(['车牌水平校正旋转角: ',num2str(angle),'度'] ,'Color','r');%显示车牌的旋转角度
%===============水平旋转校正结束=============
file_name1=strcat('..\车牌标准图\灰度车牌','.jpg');%存储标准子图
imwrite(subcol,file_name1,'jpg') %存储标准子图
waitbar(7/20,h,'请等待，识别中...');

%=========================车牌去平行框=================================
%旋转车牌后重新计算车牌水平投影，去掉车牌水平边框，获取字符高度
% histcol1=sum(sbw); %计算垂直投影
histrow=sum(sbw'); %计算水平投影
% !figure,subplot(2,1,1),bar(histcol1),title('垂直投影（旋转后）');
% !subplot(2,1,2),bar(histrow),title('水平投影（旋转后）');


% 对水平投影进行峰谷分析
meanrow=mean(histrow);%求水平投影的平均值
minrow=min(histrow);%求水平投影的最小值
levelrow=(meanrow+minrow)/2;%!!!!!!!!!!可修改%求水平投影de波峰波谷阈值
% !figure,subplot(2,1,1),bar(histrow),hold on,plot([1,hight],[levelrow,levelrow]),hold off;title('水平投影（含边框）');%输出水平投影和二值图
% !subplot(2,1,2),imshow(sbw),title('车牌二值子图');

% 计算水平投影的上升点、下降点、波峰宽度、波谷宽度
count1=0;
l=1;
markrow=zeros(1,hight+1);
markrow1=zeros(1,hight+1);
for k=1:hight
    if histrow(k)<=levelrow   %波谷处                          
        count1=count1+1;                                
    else 
        if count1>=1
            markrow(l)=k;%大于阈值的起点（该点已大于阈值)（上升点）
            markrow1(l)=count1;%谷宽度（下降点至下一个上升点）
            l=l+1;
        end
        count1=0;
    end
end
markrow(l)=hight;%最后一段的波谷补充
markrow1(l)=count1;
markrow2=diff(markrow);%峰距离（上升点至下一个上升点）
[~,n1]=size(markrow2);
%l=0;
markrow3=markrow(2:n1+1)-markrow1(2:n1+1);%大于阈值的终点（该点小于阈值）（下降点）
markrow4=markrow3-markrow(1:n1);%峰宽度（上升点至下降点）
markrow5=markrow3-double(uint16(markrow4./2));%峰中心位置（取整，转换格式）
for k=1:n1
    markrow3(k)=markrow(k+1)-markrow1(k+1);%大于阈值的终点（该点小于阈值）（下降点）
    markrow4(k)=markrow3(k)-markrow(k);%峰宽度（上升点至下降点）
    markrow5(k)=markrow3(k)-double(uint16(markrow4(k)/2));%峰中心位置（取整，转换格式）
end

%去水平（上下）边框,获取字符高度
[~,findc]=max(markrow2);%取字符高度（最大峰高度）的位置（第n个峰）
rowtop=markrow(findc);%字符起始高度
rowbot=markrow3(findc);%字符结束高度
% sbw2=sbw(rowtop:rowbot,:);  %二值子图为(rowbot-rowtop+1)行
subcol=subcol(rowtop:rowbot,:);%灰度子图为(rowbot-rowtop+1)行
% maxhight=rowbot-rowtop+1;   %字符高度(rowbot-rowtop+1)
% histcol=sum(sbw2);  %计算垂直投影
% meancol=mean(histcol);%求垂直投影的平均值
% mincol=min(histcol);%求垂直投影的最小值
% levelcol=(meancol+mincol)/4;!!!!!!求垂直投影阈值%求垂直投影的1/4
% figure,subplot(2,1,1),bar(histcol);hold on;plot([1,width],[levelcol,levelcol]);hold off;title('垂直投影（去水平边框后）');%输出车牌的垂直投影图像
%!figure,subplot(2,1,1),imshow(sbw2); title(['车牌字符高度： ',int2str(maxhight)],'Color','r');%输出车牌字符高度
%========================去水平边框结束===========================
waitbar(9/20,h,'请等待，识别中...');  

%==============垂直校正旋转=================
grd4 = edge(subcol,'canny'); 
% grd4 = bwmorph(grd4,'thicken'); 
%subplot(2,1,2),imshow(grd4),title('强化边界去水平边框图')%输出车牌强化边界
theta = -20:0.01:20;
%%进行radon变换来精确计算倾斜角度
% im_g_cp = rgb2gray(grd4);
im_g_cp =grd4;
%%对灰度图进行radon变换
[R,~] = radon(im_g_cp,theta);
%%利用求导和累加，来得到变化幅度最大的角度
diff_R = diff(R');
abs_diff_R = abs(diff_R);
sum_abs_diff_R = sum(abs_diff_R,2);

%%找出最大值以及其下标
[~,indices] = max(sum_abs_diff_R);
acc_angle = theta(1,indices-1);

%%绘制出落差幅度变化趋势
% figure,
% stem(sum_abs_diff_R);
% title('落差幅度变化趋势')
% xlabel('角度')
% ylabel('落差幅度')
% text(indices,Max_diff,'★')

%%依据计算出来的角度进行图像校正
% figure,
% rotated_im = imrotate(subcol,-acc_angle);
% imshow(rotated_im)
angle_x=tand(acc_angle);
H=[1,0,0; -angle_x,1,0;0,0,1];
T=maketform('affine',H);
subcol=imtransform(subcol,T,'bilinear');
% sbw=imtransform(sbw,T,'bilinear');
% figure,subplot(2,1,1),imshow(subcol);title('车牌灰度子图(垂直校正旋转后)');%输出车牌旋转后的灰度图像标题显示车牌灰度子图
% subplot(2,1,2),imshow(sbw);title(['车牌二值子图（垂直校正旋转后）',int2str(angle_x),'°'],'color','r')%输出车牌旋转后的灰度图像
%================垂直校正旋转结束============
waitbar(12/20,h,'请等待，识别中...');  

%====================去垂直边框，字符分割======================
sbw=im2bw(subcol);
histcol=sum(sbw);  %计算垂直投影
meancol=mean(histcol);%求垂直投影的平均值
mincol=min(histcol);%求垂直投影的最小值
levelcol=(meancol+mincol)/4;%!!!!!!求垂直投影阈值%求垂直投影的1/4
levelcol1=(meancol+mincol)/3;%!!!!!!求垂直投影阈值%求垂直投影的1/3
%!figure,bar(histcol);hold on;plot([1,width],[levelcol,levelcol]);hold off;title('垂直投影（去水平边框后）');%输出车牌的垂直投影图像

%计算车牌垂直投影：上升点，谷宽度，峰距离，下降点，峰宽度，峰位置，峰间距离
% meancol=mean(histcol);%求垂直投影的平均值
% mincol=min(histcol);%求垂直投影的平均值
% levelcol=(meancol+mincol)/2;!!!!!!求垂直投影阈值%求垂直投影的1/4
count1=0;
l=1;
for k=1:width
    if histcol(k)<=levelcol 
        count1=count1+1;
    else 
        if count1>=1
            markcol(l)=k; %字符上升点
            markcol1(l)=count1; %谷宽度（下降点至下一个上升点）
            l=l+1;
        end
        count1=0;
    end
end
markcol(l)=width;%最后一段的波谷补充
markcol1(l)=count1;
markcol2=diff(markcol);%峰距离（上升点至下一个上升点）
[~,n1]=size(markcol2);
markcol3=markcol(2:n1+1)-markcol1(2:n1+1);%大于阈值的终点（该点小于阈值）（下降点）
markcol4=markcol3-markcol(1:n1);%峰宽度（上升点至下降点）
for i=1:n1
   [markcol7(i),markcol8(i)]=max(histcol(markcol(i):markcol3(i)));%7-各峰最大值，8各峰最大值位置
   markcol8(i)=markcol8(i)+markcol(i)-1;
end
rem_position=find(markcol7<=levelcol1);
markcol4(rem_position)=[];%去除铆钉
markcol3(rem_position)=[];
markcol2(rem_position)=[];
markcol5=markcol3-double(uint16(markcol4./2));%峰中心位置（取整，转换格式）
markcol6=diff(markcol5); %字符中心距离（字符中心点至下一个字符中心点）
[~,findmax]=max(markcol6); %查找最大值，即为第二字符与第三字符中心距离
markcol6(findmax)=0;
maxwidth=max(markcol6);%查找最大值，即为最大字符宽度
%提取分割字符，去除垂直边框，并变换为22行*14列标准子图
l=1;
[~,n2]=size(sbw);
cleft_most=[];
cright_most=[];
%!figure;
for k=findmax-1:findmax+5%取间隔最大处即铆钉处（第二个字符处）前一个即后五个，去除垂直边框
        cleft=markcol5(k)-maxwidth/2;%字符左边框向量
        cright=markcol5(k)+maxwidth/2-2;%字符右边框向量
        if cleft<1%左边框修正（无垂直边框的情况）
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%右边框修正（无垂直边框的情况）
            cright=n2;
            cleft=n2-maxwidth;
        end
        if l==1
            cleft_most=cleft;
        end
        if l==7
            cright_most=cright;
        end
        SegGray=sbw(:,cleft:cright);
%         subcol1=subcol(:,cleft:cright);
        SegBw2 = imresize(SegGray,[size_x size_y],'bilinear');%变换为40行*20列标准子图 双线性插值        %2013a默认方式'bicubic'（双三次插值）
        [SegBw2]=pre_processing(SegBw2);
%!        subcol2=imresize(subcol1,[size_x size_y],'bilinear');
%!        subplot(3,n1,l),imshow(SegGray);%画每个原分割子图，n1~=7
%!        if l==4%居中显示标题
%!            title(['车牌字符宽度： ',int2str(maxwidth)],'Color','r');
%!        end
%!        subplot(3,n1,n1+l),imshow(SegBw2);title(int2str(l),'Color','r');%画每个标准子图,n1~=7   
%!        subcol2=im2bw(subcol2);
%!        subplot(3,n1,2*n1+l),imshow(subcol2);
        file_name=strcat('..\切割子图\切割子图',int2str(l),'.jpg');%存储标准子图
        imwrite(SegBw2,file_name,'jpg') %存储标准子图
        l=l+1;
end
%handles = guihandles(gcf);
%hiscol_SegGray1=sum(SegGray1);
%!figure,subplot(2,1,1),bar(hiscol_SegGray1),title('去垂直边框垂直投影');
%!subplot(2,1,2),imshow(SegGray1),title('去垂直边框灰度图');
%=========================去垂直边框，字符分割结束====================================
file_name1=strcat('..\车牌标准图\二值车牌','.jpg');%存储标准子图
if isempty(cleft_most)||isempty(cright_most)
    msgbox('识别错误','Author','');
    return
end
SegGray1=sbw(:,cleft_most:cright_most);
imwrite(SegGray1,file_name1,'jpg') %存储标准子图
waitbar(13/20,h,'请等待，识别中...');  

%=================================字符识别============================================
switch shibiefangfa
    case 1 %模板识别
%将计算计算获取的字符图像与样本库进行匹配，自动识别出字符代码。
liccode=char(['0':'9' 'A':'Z' '京辽鲁苏豫粤浙']); %建立自动识别字符代码表  
length_code=length(liccode);
l=1;
[~,n2]=size(sbw);
RegCode=char(zeros(1,14));
%!figure,
for k=findmax-1:findmax+5%字符开始
       cleft=markcol5(k)-maxwidth/2;%字符左边框向量
        cright=markcol5(k)+maxwidth/2-2;%字符右边框向量
        if cleft<1%左边框修正（无垂直边框的情况）
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%右边框修正（无垂直边框的情况）
            cright=n2;
            cleft=n2-maxwidth;
        end
        SegBw1=subcol(:,cleft:cright);
        SegBw2=imresize(SegBw1,[size_x size_y],'bilinear');%变换为40行*20列标准子图
        [SegBw2]=pre_processing(SegBw2);
        SegBw2=im2bw(SegBw2);
%!        subplot(1,n1,k),imshow(SegBw2);
        if l==1                 %第一位汉字识别
            kmin=37;
            kmax=length_code;
        elseif l==2             %第二位 A~Z 字母识别
            kmin=11;
            kmax=36;
        elseif l>=3&&l<=7      %第3-7位 0~9  A~Z字母和数字识别
            kmin=1;
            kmax=36;
        else                    %第五～七位 0~9 数字识别
            kmin=1;
            kmax=10;
        end
        Differences=zeros(size(kmin:kmax));
        
%                 for k2=kmin:kmax
%             fname=strcat('D:\users\Desktop\SRTP\程序\kevin\模板子图\',liccode(k2),'.jpg');
%             SamBw2 = imread(fname);  
%             SamBw2=imresize(SamBw2,[size_x size_y]);%变换为40行*20列标准子图
%             SubBw2 = uint8(SegBw2)-SamBw2;
%             Dmax=sum(sum(SubBw2~=0,2),1);
%             Error(k2)=Dmax;
%                 end
%         Error1=Error(kmin:kmax);%比较误差
%         MinError=min(Error1);%取误差的最小值
%         findc=find(Error1==MinError);%查找最小误差的图像
        
        for k2=kmin:kmax
            file_name=strcat('..\模板子图\',liccode(k2),'.jpg');
%             SamBw2 =imread(file_name);
            SamBw2 = imresize(imread(file_name),[size_x size_y],'bilinear');   
            SamBw2=im2bw(SamBw2);
            SubBw2 = SegBw2-SamBw2;%防止模板子图不符合大小,计算标准子图与模板子图的差别
            Differences(k2-kmin+1)=sum(sum(SubBw2~=0,2),1);%统计差的非零个数，即区别个数
        end
        MinError=min(Differences);%取误差的最小值
        findc=find(Differences==MinError);%查找最小误差的图像
        RegCode(l)=liccode(findc(1)+kmin-1);%最符合的标准子图位置
%         RegCode(l*2-1)=liccode(findc(1)+kmin-1);%最符合的标准子图位置
%         RegCode(l*2)=' ';%输出最小误差图像
        l=l+1;
end
% title (['识别车牌号码:', RegCode],'Color','r');

    case 2%神经网络识别
        [RegCode]=bp_sim();
    case 3%欧式距离识别
%将计算计算获取的字符图像与样本库进行匹配，自动识别出字符代码。
liccode=char(['0':'9' 'A':'Z' '京辽鲁苏豫粤浙']); %建立自动识别字符代码表  
length_code=length(liccode);
l=1;
[~,n2]=size(sbw);
RegCode=char(zeros(1,14));
%!figure,
for k=findmax-1:findmax+5%字符开始
       cleft=markcol5(k)-maxwidth/2;%字符左边框向量
        cright=markcol5(k)+maxwidth/2-2;%字符右边框向量
        if cleft<1%左边框修正（无垂直边框的情况）
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%右边框修正（无垂直边框的情况）
            cright=n2;
            cleft=n2-maxwidth;
        end
        SegBw1=subcol(:,cleft:cright);
        SegBw2=imresize(SegBw1,[size_x size_y],'bilinear');%变换为40行*20列标准子图
        [SegBw2]=pre_processing(SegBw2);
        SegBw2=im2bw(SegBw2);
%!        subplot(1,n1,k),imshow(SegBw2);
        if l==1                 %第一位汉字识别
            kmin=37;
            kmax=length_code;
        elseif l==2             %第二位 A~Z 字母识别
            kmin=11;
            kmax=36;
        elseif l>=3&&l<=7      %第3-7位 0~9  A~Z字母和数字识别
            kmin=1;
            kmax=36;
        else                    %第五～七位 0~9 数字识别
            kmin=1;
            kmax=10;
        end
%         Differences=zeros(size(kmin:kmax));
        ll=1;
        p=zeros(kmax-kmin,1);
        for k2=kmin:kmax
            file_name=strcat('..\模板子图\',liccode(k2),'.jpg');
%             SamBw2 =imread(file_name);
            SamBw2 = imresize(imread(file_name),[size_x size_y],'bilinear');   
            SamBw2=im2bw(SamBw2);
            temp_1=bwdist(SegBw2);
            temp_2=bwdist(SamBw2);
            p(ll)=corr2(temp_2,temp_1);
            ll=ll+1;
        end
        [~,I]=max(abs(p(:)));
        clear p;
%         ll=1;
        RegCode(l)=liccode(I+kmin-1);%最符合的标准子图位置
        l=l+1;
end
        
        
        
end


%储存识别号码
fid=fopen('..\识别车牌号.txt','w');
fprintf(fid,'%s\n',RegCode);
fclose(fid);
%弹出对话框显示识别号码
waitbar(14/20,h,'请等待，识别中...');    
waitbar(1,h,'已完成');
delete(h);
msgbox(['车牌号：',RegCode],'Author','');
% 读取声音
if shibiefangfa~=2 %神经网络不放音
for k=1:7
    file_name=strcat('..\声音文件\',RegCode(k),'.wav');%相对路径
    wavplay(wavread(file_name),44100);
end
end
end


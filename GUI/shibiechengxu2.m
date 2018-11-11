function [ ] = shibiechengxu2( )
% clear;
% close all;
%运行速度：5.418s 水平垂直倾斜改为透视变换
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
supplement_x=10;%从连通区到车牌矩形提取的补充（能否自适应）
supplement_y=15;
%diff_ydata=2;%车牌旋转校正，数据修正，使预拟合曲线光滑，防止跳变（像素点）
size_x=40;%子图标准大小x
size_y=20;%子图标准大小y
shot_height=1080;
shot_width=1920;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hbar=waitbar(0,'请稍等，正在开始识别文件...');

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
waitbar(1/20,hbar,'请等待，识别中...');

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
        level=(level1+level2)/2;%!!!取平均值
end
bw2=im2bw(rgb_image);%转换图像为二进制图像
bw2=double(bw2);%转换为双精度型
% !figure,imshow(bw2),title('二值图');
waitbar(2/20,hbar,'请等待，识别中...');

%对二值图像作开闭操作进行滤波(腐蚀膨胀)
grd=edge(bw2,'canny');%用canny算子识别强度图像中的边界
% !figure,imshow(grd);title('图像边缘提取');%输出图像边缘
bg1=imclose(grd,strel('rectangle',[rectangle_x,rectangle_y]));%!!!!!!开闭运算修改%取矩形框的闭运算
% !figure,imshow(bg1);title('图像闭运算[5,19]');%输出闭运算的图像
bg3=imopen(bg1,strel('rectangle',[rectangle_x,rectangle_y]));%取矩形框的开运算
%!figure,imshow(bg3);title('图像开运算[5,19]');%输出开运算的图像
bg2=imopen(bg3,strel('rectangle',[rectangle_y,1]));%取矩形框的开运算
%!figure,imshow(bg2);title('图像开运算[19,1]');%输出开运算的图像
waitbar(3/20,hbar,'请等待，识别中...');

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
waitbar(4/20,hbar,'请等待，识别中...');

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
endrow=startrow+hight+supplement_y;
endcol=startcol+width+supplement_x;
sbw1=bw2(startrow:endrow,startcol:endcol); %获取车牌二值子图
subcol1=gray_image(startrow:endrow,startcol:endcol);%获取车牌灰度子图
% grd1=grd(startrow:endrow,startcol:endcol);%获取车牌边界子图
% !figure,subplot(2,1,1),imshow(subcol1);title('车牌灰度子图');%显示灰度图像
% !subplot(2,1,2),imshow(sbw1);title('车牌二值子图');%显示车牌的二值图
%===============================预处理结束=======================================
waitbar(5/20,hbar,'请等待，识别中...');

%%
%===============================透视旋转=========================================
file_name=strcat('..\车牌标准图\灰度车牌','.jpg');%存储标准子图
imwrite(subcol1,file_name,'jpg') %存储标准子图
img=subcol1;
[M,N] = size(img);
allow_diff_y=M/10;%y方向允许的下降高度M/10
allow_diff_x=4;%x方向允许的下降高度N/20
bw= im2bw(img);
% figure,imshow(mat2gray(img))
% figure,imshow(bw);
% hold on,plot([N/10,N/10],[1,M]), %减少边角的影响
%         plot([9*N/10,9*N/10],[1,M]),%减少边角的影响
%         plot([1,N],[M/5,M/5]), %减少边角的影响
%         plot([1,N],[4*M/5,4*M/5]), %减少边角的影响
% hold off
%检测上边线
l=1;
 for j=ceil(N/10):2:ceil(9*N/10)
     for i=1:M
         if bw(i,j)==1
             if l==1
             ydata(l)=i;
             diff_data=0;
             else
             diff_data=abs(i-ydata(l-1));
             end
             if diff_data<=allow_diff_y
             xdata(l)=j;%j值，图中横坐标：2
             ydata(l)=i;%i值，图中纵坐标
             l=l+1;
             break
             end
         end
     end
 end
 fresult1=polyfit(xdata,ydata,1);%一次线性拟合
%  z=polyval(fresult1,xdata);
%  figure,
%  imshow(bw),hold on, plot(xdata',ydata','r',xdata,z,'b'),
 %检测下边线
%检测上边线
clear xdata ydata
l=1;
 for j=ceil(N/10):2:ceil(9*N/10)
     for i=M:-1:1
         if bw(i,j)==1
             if l==1
             ydata(l)=i;
             diff_data=0;
             else
             diff_data=abs(i-ydata(l-1));
             end
             if diff_data<=allow_diff_y
             xdata(l)=j;%j值，图中横坐标：2
             ydata(l)=i;%i值，图中纵坐标
             l=l+1;
             break
             end
         end
     end
 end
 fresult2=polyfit(xdata,ydata,1);%一次线性拟合
%  z=polyval(fresult2,xdata);
% hold on, plot(xdata',ydata','r',xdata,z,'b'),
%检测左边线
clear xdata ydata
l=1;
 for i=ceil(M/5):1:ceil(4*M/5)
     for j=1:N
         if bw(i,j)==1
             if l==1
             ydata(l)=i;
             diff_data=0;
             else
             diff_data=abs(i-ydata(l-1));
             end
             if diff_data<=allow_diff_x
             xdata(l)=j;%j值，图中横坐标：2
             ydata(l)=i;%i值，图中纵坐标
             l=l+1;
             break
             end
         end
     end
 end
 fresult3=polyfit(xdata,ydata,1);%一次线性拟合
%  z=polyval(fresult3,xdata);
%  hold on, plot(xdata',ydata','r',xdata,z,'b'),
 %检测右边线
clear xdata ydata
l=1;
 for i=ceil(M/5):1:ceil(4*M/5)
     for j=N:-1:1
         if bw(i,j)==1
             if l==1
             ydata(l)=i;
             diff_data=0;
             else
             diff_data=abs(i-ydata(l-1));
             end
             if diff_data<=allow_diff_x
             xdata(l)=j;%j值，图中横坐标：2
             ydata(l)=i;%i值，图中纵坐标
             l=l+1;
             break
             end
         end
     end
 end
 fresult4=polyfit(xdata,ydata,1);%一次线性拟合
%  z=polyval(fresult4,xdata);
%  hold on, plot(xdata',ydata','r',xdata,z,'b'),hold off
 clear xdata ydata
 %求交点
 dot=zeros(4,2);
 syms x y f1 f2 f3 f4
 f1=poly2sym(fresult1,x);
 f2=poly2sym(fresult2,x);
 f3=poly2sym(fresult3,x);
 f4=poly2sym(fresult4,x);
 %左上
dot(1,1)=round(double(solve(f1==f3)));
dot(1,2)=round(polyval(fresult1,dot(1,1)));
%右上
dot(2,1)=round(double(solve(f1==f4)));
dot(2,2)=round(polyval(fresult1,dot(2,1)));
%左下
dot(3,1)=round(double(solve(f2==f3)));
dot(3,2)=round(polyval(fresult2,dot(3,1)));
%右下
dot(4,1)=round(double(solve(f2==f4)));
dot(4,2)=round(polyval(fresult2,dot(4,1)));
waitbar(7/20,hbar,'请等待，识别中...');
%%
w=round(sqrt((dot(1,1)-dot(2,1))^2+(dot(1,2)-dot(2,2))^2));     %从原四边形获得新矩形宽
h=round(sqrt((dot(1,1)-dot(3,1))^2+(dot(1,2)-dot(3,2))^2));     %从原四边形获得新矩形高

y=[dot(1,1) dot(2,1) dot(3,1) dot(4,1)];        %四个原顶点
x=[dot(1,2) dot(2,2) dot(3,2) dot(4,2)];

%这里是新的顶点，我取的矩形,也可以做成其他的形状
%大可以原图像是矩形，新图像是从dot中取得的点组成的任意四边形.:)
Y=[dot(1,1) dot(1,1) dot(1,1)+h dot(1,1)+h];     
X=[dot(1,2) dot(1,2)+w dot(1,2) dot(1,2)+w];

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';   %变换后的四个顶点，方程右边的值
%联立解方程组，方程的系数
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);             
   0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
   0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
   0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
   0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=A\B;        %用四点求得的方程的解，也是全局变换系数
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
g=fa(7);h=fa(8);

rot=[d e f;
     a b c;
     g h 1];        %公式中第一个数是x,Matlab第一个表示y，所以我矩阵1,2行互换了

pix1=rot*[1 1 1]'/(g*1+h*1+1);  %变换后图像左上点
pix2=rot*[1 N 1]'/(g*1+h*N+1);  %变换后图像右上点
pix3=rot*[M 1 1]'/(g*M+h*1+1);  %变换后图像左下点
pix4=rot*[M N 1]'/(g*M+h*N+1);  %变换后图像右下点

height=round(max([pix1(1) pix2(1) pix3(1) pix4(1)])-min([pix1(1) pix2(1) pix3(1) pix4(1)]));     %变换后图像的高度
width=round(max([pix1(2) pix2(2) pix3(2) pix4(2)])-min([pix1(2) pix2(2) pix3(2) pix4(2)]));      %变换后图像的宽度
% clear subcol sbw

delta_y=round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])));            %取得y方向的负轴超出的偏移量
delta_x=round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));            %取得x方向的负轴超出的偏移量

for i = 1-delta_y:height-delta_y                        %从变换图像中反向寻找原图像的点，以免出现空洞，和旋转放大原理一样
    for j = 1-delta_x:width-delta_x
        pix=rot\[i j 1]';       %求原图像中坐标，因为[YW XW W]=fa*[y x 1],所以这里求的是[YW XW W],W=gy+hx+1;
        pix=[g*pix(1)-1 h*pix(1);g*pix(2) h*pix(2)-1]\[-pix(1) -pix(2)]'; %相当于解[pix(1)*(gy+hx+1) pix(2)*(gy+hx+1)]=[y x],这样一个方程，求y和x，最后pix=[y x];
        
        if pix(1)>=0.5 && pix(2)>=0.5 && pix(1)<=M && pix(2)<=N
            subcol(i+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));     %最邻近插值,也可以用双线性或双立方插值
            sbw(i+delta_y,j+delta_x)=sbw1(round(pix(1)),round(pix(2))); 
        end  
    end
end
% sbw=im2bw(imgn);
% figure,subplot(2,1,1),imshow(subcol);
%        subplot(2,1,2),imshow(sbw);
waitbar(10/20,hbar,'请等待，识别中...');

%%
%=========================车牌去平行框=================================
%旋转车牌后重新计算车牌水平投影，去掉车牌水平边框，获取字符高度
[hight,width]=size(sbw);
% histcol1=sum(sbw); %计算垂直投影
histrow=sum(sbw'); %计算水平投影
% figure,subplot(2,1,1),bar(histcol1),title('垂直投影（旋转后）');
% subplot(2,1,2),bar(histrow),title('水平投影（旋转后）');


% 对水平投影进行峰谷分析
meanrow=mean(histrow);%求水平投影的平均值
minrow=min(histrow);%求水平投影的最小值
levelrow=(meanrow+minrow)/2;%!!!!!!!!!!可修改%求水平投影de波峰波谷阈值
% figure,subplot(2,1,1),bar(histrow),hold on,plot([1,hight],[levelrow,levelrow]),hold off;title('水平投影（含边框）');%输出水平投影和二值图
% subplot(2,1,2),imshow(sbw),title('车牌二值子图');


% 计算水平投影的上升点、下降点、波峰宽度、波谷宽度
count1=0;
l=1;
markrow=zeros(1,hight);
markrow1=zeros(1,hight);
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
l=0;
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
% figure,subplot(2,1,1),imshow(sbw2); title(['车牌字符高度： ',int2str(maxhight)],'Color','r');%输出车牌字符高度
%========================去水平边框结束===========================
waitbar(11/20,hbar,'请等待，识别中...');  

%%
%====================去垂直边框，字符分割======================
sbw=im2bw(subcol);
histcol=sum(sbw);  %计算垂直投影
meancol=mean(histcol);%求垂直投影的平均值
mincol=min(histcol);%求垂直投影的最小值
levelcol=(meancol+mincol)/4;%!!!!!!求垂直投影阈值%求垂直投影的1/4
levelcol1=(meancol+mincol)/3;%!!!!!!求垂直投影阈值%求垂直投影的1/3
% figure,bar(histcol);hold on;plot([1,width],[levelcol,levelcol]);plot([1,width],[levelcol1,levelcol1]);hold off;title('垂直投影（去水平边框后）');%输出车牌的垂直投影图像

%计算车牌垂直投影：上升点，谷宽度，峰距离，下降点，峰宽度，峰位置，峰间距离
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
% figure;
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
%         if l==1
%             cleft_most=cleft;
%         end
%         if l==7
%             cright_most=cright;
%         end
        SegGray=sbw(:,cleft:cright);
%         subcol1=subcol(:,cleft:cright);
        SegBw2 = imresize(SegGray,[size_x size_y],'bilinear');%变换为40行*20列标准子图 双线性插值        %2013a默认方式'bicubic'（双三次插值）
        [SegBw2]=pre_processing(SegBw2);
%         subcol2=imresize(subcol1,[size_x size_y],'bilinear');
%         subplot(3,n1,l),imshow(SegGray);%画每个原分割子图，n1~=7
%         if l==4%居中显示标题
%             title(['车牌字符宽度： ',int2str(maxwidth)],'Color','r');
%         end
%        subplot(3,n1,n1+l),imshow(SegBw2);title(int2str(l),'Color','r');%画每个标准子图,n1~=7   
%         subcol2=im2bw(subcol2);
%         subplot(3,n1,2*n1+l),imshow(subcol2);
        file_name=strcat('..\切割子图\切割子图',int2str(l),'.jpg');%存储标准子图
        imwrite(SegBw2,file_name,'jpg') %存储标准子图
        l=l+1;
end
% SegGray1=sbw(:,cleft_most:cright_most);
% hiscol_SegGray1=sum(SegGray1);
% figure,subplot(2,1,1),bar(hiscol_SegGray1),title('去垂直边框垂直投影');
% subplot(2,1,2),imshow(SegGray1),title('去垂直边框灰度图');
%=========================去垂直边框，字符分割结束====================================
waitbar(13/20,hbar,'请等待，识别中...'); 
%%
%=================================字符识别============================================
switch shibiefangfa
    case 1 %模板识别
%将计算计算获取的字符图像与样本库进行匹配，自动识别出字符代码。
liccode=char(['0':'9' 'A':'H' 'J':'N' 'P':'Z' '京辽鲁苏豫粤浙']); %建立自动识别字符代码表  
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
            kmin=35;
            kmax=length_code;
        elseif l==2             %第二位 A~Z 字母识别
            kmin=11;
            kmax=34;
        elseif l>=3&&l<=7      %第3-7位 0~9  A~Z字母和数字识别
            kmin=1;
            kmax=34;
        else                    %第五～七位 0~9 数字识别
            kmin=1;
            kmax=10;
        end
        Differences=zeros(size(kmin:kmax));     
        for k2=kmin:kmax
            file_name=strcat('..\模板子图\',liccode(k2),'.jpg');
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
            [SamBw2]=pre_processing(SamBw2);
%             SamBw2=im2bw(SamBw2,0.5);
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
waitbar(14/20,hbar,'请等待，识别中...');    
waitbar(1,hbar,'已完成');
delete(hbar);
msgbox(['车牌号：',RegCode],'Author','');
% 读取声音
% if shibiefangfa~=2 %神经网络不放音
% for k=1:7
%     file_name=strcat('..\声音文件\',RegCode(k),'.wav');%相对路径
%     wavplay(wavread(file_name),44100);
% end
% end
end


function [ ] = shibiechengxu2( )
% clear;
% close all;
%�����ٶȣ�5.418s ˮƽ��ֱ��б��Ϊ͸�ӱ任
global im
global erzhi_yuzhi
global shibiefangfa
global tezheng
global RegCode
%from kevin
erzhi_yuzhi
shibiefangfa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_s=14;    %strel����
rectangle_x=20;     %�˲����ο��С
rectangle_y=100;
rectangle_width=[98,500];%��ͨ���򣨳��ƣ������С
rectangle_hight=[25,200];
pr_rato=[2,10];%��ͨ���򣨳��ƣ��������
supplement_x=10;%����ͨ�������ƾ�����ȡ�Ĳ��䣨�ܷ�����Ӧ��
supplement_y=15;
%diff_ydata=2;%������תУ��������������ʹԤ������߹⻬����ֹ���䣨���ص㣩
size_x=40;%��ͼ��׼��Сx
size_y=20;%��ͼ��׼��Сy
shot_height=1080;
shot_width=1920;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hbar=waitbar(0,'���Եȣ����ڿ�ʼʶ���ļ�...');

%=================ͼ�������Ԥ����=============================
rgb_image=im;
% !figure,imshow(rgb_image),title('ԭ��ɫͼ��');%��ʾԭͼ��
gray_image=rgb2gray(rgb_image);%ԭͼת�Ҷ�ͼ
% !figure,imshow(gray_image),title('ԭ�Ҷ�ͼ');%��ʾ�Ҷ�ͼ

%�����ǿͼ��ȥ��������
s=strel('disk',parameter_s);%!!!!!!!strel�����޸�
Bgray_image=imopen(gray_image,s);%���������
% !figure,imshow(Bgray_image),title('����ͼ��')%��ʾ����ͼ
Egray_image=imsubtract(gray_image,Bgray_image);%��ͼ�������ǿͼ
% !figure,imshow(Egray_image),title('��ǿͼ��')%��ʾ��ǿͼ
waitbar(1/20,hbar,'��ȴ���ʶ����...');

%��ö�ֵ����ֵ
image_max=double(max(max(Egray_image)));%��ǿͼ���ֵ�����˫������
image_min=double(min(min(Egray_image)));%��ǿͼ��Сֵ�����˫������
level1=(image_max-(image_max-image_min)/3)/255;%�����ֵ�޸�%��������ֵ
level2 = graythresh(gray_image);%OSTU������ֵ
switch erzhi_yuzhi
    case 1
        level=level1;
    case 2
        level=level2;
    case 3
        level=(level1+level2)/2;%!!!ȡƽ��ֵ
end
bw2=im2bw(rgb_image);%ת��ͼ��Ϊ������ͼ��
bw2=double(bw2);%ת��Ϊ˫������
% !figure,imshow(bw2),title('��ֵͼ');
waitbar(2/20,hbar,'��ȴ���ʶ����...');

%�Զ�ֵͼ�������ղ��������˲�(��ʴ����)
grd=edge(bw2,'canny');%��canny����ʶ��ǿ��ͼ���еı߽�
% !figure,imshow(grd);title('ͼ���Ե��ȡ');%���ͼ���Ե
bg1=imclose(grd,strel('rectangle',[rectangle_x,rectangle_y]));%!!!!!!���������޸�%ȡ���ο�ı�����
% !figure,imshow(bg1);title('ͼ�������[5,19]');%����������ͼ��
bg3=imopen(bg1,strel('rectangle',[rectangle_x,rectangle_y]));%ȡ���ο�Ŀ�����
%!figure,imshow(bg3);title('ͼ������[5,19]');%����������ͼ��
bg2=imopen(bg3,strel('rectangle',[rectangle_y,1]));%ȡ���ο�Ŀ�����
%!figure,imshow(bg2);title('ͼ������[19,1]');%����������ͼ��
waitbar(3/20,hbar,'��ȴ���ʶ����...');

%�Զ�ֵͼ�����������ȡ���������������������������������������Ƚϣ���ȡ��������
[L,num] = bwlabel(bg2,8);%��ע������ͼ���������ӵĲ��֣�LΪ��ͨ�������numΪ��ͨ������Ŀ��
Feastats = regionprops(L,'basic');%����ͼ�������ϵ�������ߴ磨'Area'�Ǳ������������ͼ����������������ܸ�����'BoundingBox'��1��ndims(L)*2�е���������������Ӧ�������С���Ρ�BoundingBox ��ʽΪ [ul_corner width]������ ul_corner �� [x y z ...] ��������ʽ�����߽���ӵ����Ͻǡ�boxwidth �� [x_width y_width ...] ��ʽָ���߽��������ÿ��ά������ĳ��ȡ���
                                  %'Centroid'��1��ndims(L)�е�����������ÿ����������ģ����ģ���
% Area=[Feastats.Area];%�������
BoundingBox=[Feastats.BoundingBox];%[x y width height]���ƵĿ�ܴ�С
RGB_image2= label2rgb(L, 'spring', 'k', 'noshuffle'); %��־ͼ����RGBͼ��ת��
file_name1=strcat('..\���Ʊ�׼ͼ\���Ʊ��ͼ','.jpg');%�洢��׼��ͼ
imwrite(RGB_image2,file_name1,'jpg') %�洢��׼��ͼ
% !figure,imshow(RGB_image2);title('ͼ���ɫ���');%�����ܵĲ�ɫͼ��
waitbar(4/20,hbar,'��ȴ���ʶ����...');

%������ͨ�����Ƿ���ϳ��ƴ�С
lx=0;
Getok=zeros(1,num);
for l=1:num
    width=BoundingBox((l-1)*4+3);%��ܿ�ȵļ���
    hight=BoundingBox((l-1)*4+4);%��ܸ߶ȵļ���
    if (width>rectangle_width(1) && width<rectangle_width(2) && hight>rectangle_hight(1) && hight<rectangle_hight(2))%��ܵĿ�Ⱥ͸߶ȵķ�Χ
        lx=lx+1;
        Getok(lx)=l;
    end
end
%����ɸѡ�����ͨ����ĳ���ȣ����������������ɿ��ǽ���1����ͨ�����������ʱ��
for k= 1:lx
    l=Getok(k);    
    startcol=BoundingBox((l-1)*4+1)-supplement_x;%��ͨ�������Ͻ�����x
    startrow=BoundingBox((l-1)*4+2)-supplement_y;%��ͨ�������Ͻ�����y
    if startcol<0||startrow<0
        continue
    end
    width=BoundingBox((l-1)*4+3)+supplement_x;%���ƿ�
    hight=BoundingBox((l-1)*4+4)+supplement_y;%���Ƹ�
    endcol=startcol+width;
    endrow=startrow+hight;
    if endcol>shot_width||endrow>shot_height
        continue
    end
    rato=width/hight;%���㳵�Ƴ����
    if rato>pr_rato(1) && rato<pr_rato(2)   
        break;
    end
end
endrow=startrow+hight+supplement_y;
endcol=startcol+width+supplement_x;
sbw1=bw2(startrow:endrow,startcol:endcol); %��ȡ���ƶ�ֵ��ͼ
subcol1=gray_image(startrow:endrow,startcol:endcol);%��ȡ���ƻҶ���ͼ
% grd1=grd(startrow:endrow,startcol:endcol);%��ȡ���Ʊ߽���ͼ
% !figure,subplot(2,1,1),imshow(subcol1);title('���ƻҶ���ͼ');%��ʾ�Ҷ�ͼ��
% !subplot(2,1,2),imshow(sbw1);title('���ƶ�ֵ��ͼ');%��ʾ���ƵĶ�ֵͼ
%===============================Ԥ�������=======================================
waitbar(5/20,hbar,'��ȴ���ʶ����...');

%%
%===============================͸����ת=========================================
file_name=strcat('..\���Ʊ�׼ͼ\�Ҷȳ���','.jpg');%�洢��׼��ͼ
imwrite(subcol1,file_name,'jpg') %�洢��׼��ͼ
img=subcol1;
[M,N] = size(img);
allow_diff_y=M/10;%y����������½��߶�M/10
allow_diff_x=4;%x����������½��߶�N/20
bw= im2bw(img);
% figure,imshow(mat2gray(img))
% figure,imshow(bw);
% hold on,plot([N/10,N/10],[1,M]), %���ٱ߽ǵ�Ӱ��
%         plot([9*N/10,9*N/10],[1,M]),%���ٱ߽ǵ�Ӱ��
%         plot([1,N],[M/5,M/5]), %���ٱ߽ǵ�Ӱ��
%         plot([1,N],[4*M/5,4*M/5]), %���ٱ߽ǵ�Ӱ��
% hold off
%����ϱ���
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
             xdata(l)=j;%jֵ��ͼ�к����꣺2
             ydata(l)=i;%iֵ��ͼ��������
             l=l+1;
             break
             end
         end
     end
 end
 fresult1=polyfit(xdata,ydata,1);%һ���������
%  z=polyval(fresult1,xdata);
%  figure,
%  imshow(bw),hold on, plot(xdata',ydata','r',xdata,z,'b'),
 %����±���
%����ϱ���
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
             xdata(l)=j;%jֵ��ͼ�к����꣺2
             ydata(l)=i;%iֵ��ͼ��������
             l=l+1;
             break
             end
         end
     end
 end
 fresult2=polyfit(xdata,ydata,1);%һ���������
%  z=polyval(fresult2,xdata);
% hold on, plot(xdata',ydata','r',xdata,z,'b'),
%��������
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
             xdata(l)=j;%jֵ��ͼ�к����꣺2
             ydata(l)=i;%iֵ��ͼ��������
             l=l+1;
             break
             end
         end
     end
 end
 fresult3=polyfit(xdata,ydata,1);%һ���������
%  z=polyval(fresult3,xdata);
%  hold on, plot(xdata',ydata','r',xdata,z,'b'),
 %����ұ���
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
             xdata(l)=j;%jֵ��ͼ�к����꣺2
             ydata(l)=i;%iֵ��ͼ��������
             l=l+1;
             break
             end
         end
     end
 end
 fresult4=polyfit(xdata,ydata,1);%һ���������
%  z=polyval(fresult4,xdata);
%  hold on, plot(xdata',ydata','r',xdata,z,'b'),hold off
 clear xdata ydata
 %�󽻵�
 dot=zeros(4,2);
 syms x y f1 f2 f3 f4
 f1=poly2sym(fresult1,x);
 f2=poly2sym(fresult2,x);
 f3=poly2sym(fresult3,x);
 f4=poly2sym(fresult4,x);
 %����
dot(1,1)=round(double(solve(f1==f3)));
dot(1,2)=round(polyval(fresult1,dot(1,1)));
%����
dot(2,1)=round(double(solve(f1==f4)));
dot(2,2)=round(polyval(fresult1,dot(2,1)));
%����
dot(3,1)=round(double(solve(f2==f3)));
dot(3,2)=round(polyval(fresult2,dot(3,1)));
%����
dot(4,1)=round(double(solve(f2==f4)));
dot(4,2)=round(polyval(fresult2,dot(4,1)));
waitbar(7/20,hbar,'��ȴ���ʶ����...');
%%
w=round(sqrt((dot(1,1)-dot(2,1))^2+(dot(1,2)-dot(2,2))^2));     %��ԭ�ı��λ���¾��ο�
h=round(sqrt((dot(1,1)-dot(3,1))^2+(dot(1,2)-dot(3,2))^2));     %��ԭ�ı��λ���¾��θ�

y=[dot(1,1) dot(2,1) dot(3,1) dot(4,1)];        %�ĸ�ԭ����
x=[dot(1,2) dot(2,2) dot(3,2) dot(4,2)];

%�������µĶ��㣬��ȡ�ľ���,Ҳ����������������״
%�����ԭͼ���Ǿ��Σ���ͼ���Ǵ�dot��ȡ�õĵ���ɵ������ı���.:)
Y=[dot(1,1) dot(1,1) dot(1,1)+h dot(1,1)+h];     
X=[dot(1,2) dot(1,2)+w dot(1,2) dot(1,2)+w];

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';   %�任����ĸ����㣬�����ұߵ�ֵ
%�����ⷽ���飬���̵�ϵ��
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);             
   0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
   0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
   0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
   0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=A\B;        %���ĵ���õķ��̵Ľ⣬Ҳ��ȫ�ֱ任ϵ��
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
g=fa(7);h=fa(8);

rot=[d e f;
     a b c;
     g h 1];        %��ʽ�е�һ������x,Matlab��һ����ʾy�������Ҿ���1,2�л�����

pix1=rot*[1 1 1]'/(g*1+h*1+1);  %�任��ͼ�����ϵ�
pix2=rot*[1 N 1]'/(g*1+h*N+1);  %�任��ͼ�����ϵ�
pix3=rot*[M 1 1]'/(g*M+h*1+1);  %�任��ͼ�����µ�
pix4=rot*[M N 1]'/(g*M+h*N+1);  %�任��ͼ�����µ�

height=round(max([pix1(1) pix2(1) pix3(1) pix4(1)])-min([pix1(1) pix2(1) pix3(1) pix4(1)]));     %�任��ͼ��ĸ߶�
width=round(max([pix1(2) pix2(2) pix3(2) pix4(2)])-min([pix1(2) pix2(2) pix3(2) pix4(2)]));      %�任��ͼ��Ŀ��
% clear subcol sbw

delta_y=round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])));            %ȡ��y����ĸ��ᳬ����ƫ����
delta_x=round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));            %ȡ��x����ĸ��ᳬ����ƫ����

for i = 1-delta_y:height-delta_y                        %�ӱ任ͼ���з���Ѱ��ԭͼ��ĵ㣬������ֿն�������ת�Ŵ�ԭ��һ��
    for j = 1-delta_x:width-delta_x
        pix=rot\[i j 1]';       %��ԭͼ�������꣬��Ϊ[YW XW W]=fa*[y x 1],�������������[YW XW W],W=gy+hx+1;
        pix=[g*pix(1)-1 h*pix(1);g*pix(2) h*pix(2)-1]\[-pix(1) -pix(2)]'; %�൱�ڽ�[pix(1)*(gy+hx+1) pix(2)*(gy+hx+1)]=[y x],����һ�����̣���y��x�����pix=[y x];
        
        if pix(1)>=0.5 && pix(2)>=0.5 && pix(1)<=M && pix(2)<=N
            subcol(i+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));     %���ڽ���ֵ,Ҳ������˫���Ի�˫������ֵ
            sbw(i+delta_y,j+delta_x)=sbw1(round(pix(1)),round(pix(2))); 
        end  
    end
end
% sbw=im2bw(imgn);
% figure,subplot(2,1,1),imshow(subcol);
%        subplot(2,1,2),imshow(sbw);
waitbar(10/20,hbar,'��ȴ���ʶ����...');

%%
%=========================����ȥƽ�п�=================================
%��ת���ƺ����¼��㳵��ˮƽͶӰ��ȥ������ˮƽ�߿򣬻�ȡ�ַ��߶�
[hight,width]=size(sbw);
% histcol1=sum(sbw); %���㴹ֱͶӰ
histrow=sum(sbw'); %����ˮƽͶӰ
% figure,subplot(2,1,1),bar(histcol1),title('��ֱͶӰ����ת��');
% subplot(2,1,2),bar(histrow),title('ˮƽͶӰ����ת��');


% ��ˮƽͶӰ���з�ȷ���
meanrow=mean(histrow);%��ˮƽͶӰ��ƽ��ֵ
minrow=min(histrow);%��ˮƽͶӰ����Сֵ
levelrow=(meanrow+minrow)/2;%!!!!!!!!!!���޸�%��ˮƽͶӰde���岨����ֵ
% figure,subplot(2,1,1),bar(histrow),hold on,plot([1,hight],[levelrow,levelrow]),hold off;title('ˮƽͶӰ�����߿�');%���ˮƽͶӰ�Ͷ�ֵͼ
% subplot(2,1,2),imshow(sbw),title('���ƶ�ֵ��ͼ');


% ����ˮƽͶӰ�������㡢�½��㡢�����ȡ����ȿ��
count1=0;
l=1;
markrow=zeros(1,hight);
markrow1=zeros(1,hight);
for k=1:hight
    if histrow(k)<=levelrow   %���ȴ�                          
        count1=count1+1;                                
    else 
        if count1>=1
            markrow(l)=k;%������ֵ����㣨�õ��Ѵ�����ֵ)�������㣩
            markrow1(l)=count1;%�ȿ�ȣ��½�������һ�������㣩
            l=l+1;
        end
        count1=0;
    end
end
markrow(l)=hight;%���һ�εĲ��Ȳ���
markrow1(l)=count1;
markrow2=diff(markrow);%����루����������һ�������㣩
[~,n1]=size(markrow2);
l=0;
markrow3=markrow(2:n1+1)-markrow1(2:n1+1);%������ֵ���յ㣨�õ�С����ֵ�����½��㣩
markrow4=markrow3-markrow(1:n1);%���ȣ����������½��㣩
markrow5=markrow3-double(uint16(markrow4./2));%������λ�ã�ȡ����ת����ʽ��
for k=1:n1
    markrow3(k)=markrow(k+1)-markrow1(k+1);%������ֵ���յ㣨�õ�С����ֵ�����½��㣩
    markrow4(k)=markrow3(k)-markrow(k);%���ȣ����������½��㣩
    markrow5(k)=markrow3(k)-double(uint16(markrow4(k)/2));%������λ�ã�ȡ����ת����ʽ��
end

%ȥˮƽ�����£��߿�,��ȡ�ַ��߶�
[~,findc]=max(markrow2);%ȡ�ַ��߶ȣ�����߶ȣ���λ�ã���n���壩
rowtop=markrow(findc);%�ַ���ʼ�߶�
rowbot=markrow3(findc);%�ַ������߶�
% sbw2=sbw(rowtop:rowbot,:);  %��ֵ��ͼΪ(rowbot-rowtop+1)��
subcol=subcol(rowtop:rowbot,:);%�Ҷ���ͼΪ(rowbot-rowtop+1)��
% maxhight=rowbot-rowtop+1;   %�ַ��߶�(rowbot-rowtop+1)
% histcol=sum(sbw2);  %���㴹ֱͶӰ
% meancol=mean(histcol);%��ֱͶӰ��ƽ��ֵ
% mincol=min(histcol);%��ֱͶӰ����Сֵ
% levelcol=(meancol+mincol)/4;!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/4
% figure,subplot(2,1,1),bar(histcol);hold on;plot([1,width],[levelcol,levelcol]);hold off;title('��ֱͶӰ��ȥˮƽ�߿��');%������ƵĴ�ֱͶӰͼ��
% figure,subplot(2,1,1),imshow(sbw2); title(['�����ַ��߶ȣ� ',int2str(maxhight)],'Color','r');%��������ַ��߶�
%========================ȥˮƽ�߿����===========================
waitbar(11/20,hbar,'��ȴ���ʶ����...');  

%%
%====================ȥ��ֱ�߿��ַ��ָ�======================
sbw=im2bw(subcol);
histcol=sum(sbw);  %���㴹ֱͶӰ
meancol=mean(histcol);%��ֱͶӰ��ƽ��ֵ
mincol=min(histcol);%��ֱͶӰ����Сֵ
levelcol=(meancol+mincol)/4;%!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/4
levelcol1=(meancol+mincol)/3;%!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/3
% figure,bar(histcol);hold on;plot([1,width],[levelcol,levelcol]);plot([1,width],[levelcol1,levelcol1]);hold off;title('��ֱͶӰ��ȥˮƽ�߿��');%������ƵĴ�ֱͶӰͼ��

%���㳵�ƴ�ֱͶӰ�������㣬�ȿ�ȣ�����룬�½��㣬���ȣ���λ�ã�������
count1=0;
l=1;
for k=1:width
    if histcol(k)<=levelcol 
        count1=count1+1;
    else 
        if count1>=1
            markcol(l)=k; %�ַ�������
            markcol1(l)=count1; %�ȿ�ȣ��½�������һ�������㣩
            l=l+1;
        end
        count1=0;
    end
end
markcol(l)=width;%���һ�εĲ��Ȳ���
markcol1(l)=count1;
markcol2=diff(markcol);%����루����������һ�������㣩
[~,n1]=size(markcol2);
markcol3=markcol(2:n1+1)-markcol1(2:n1+1);%������ֵ���յ㣨�õ�С����ֵ�����½��㣩
markcol4=markcol3-markcol(1:n1);%���ȣ����������½��㣩
for i=1:n1
   [markcol7(i),markcol8(i)]=max(histcol(markcol(i):markcol3(i)));%7-�������ֵ��8�������ֵλ��
   markcol8(i)=markcol8(i)+markcol(i)-1;
end
rem_position=find(markcol7<=levelcol1);
markcol4(rem_position)=[];%ȥ��í��
markcol3(rem_position)=[];
markcol2(rem_position)=[];
markcol5=markcol3-double(uint16(markcol4./2));%������λ�ã�ȡ����ת����ʽ��
markcol6=diff(markcol5); %�ַ����ľ��루�ַ����ĵ�����һ���ַ����ĵ㣩
[~,findmax]=max(markcol6); %�������ֵ����Ϊ�ڶ��ַ�������ַ����ľ���
markcol6(findmax)=0;
maxwidth=max(markcol6);%�������ֵ����Ϊ����ַ����
%��ȡ�ָ��ַ���ȥ����ֱ�߿򣬲��任Ϊ22��*14�б�׼��ͼ
l=1;
[~,n2]=size(sbw);
% figure;
for k=findmax-1:findmax+5%ȡ�����󴦼�í�������ڶ����ַ�����ǰһ�����������ȥ����ֱ�߿�
        cleft=markcol5(k)-maxwidth/2;%�ַ���߿�����
        cright=markcol5(k)+maxwidth/2-2;%�ַ��ұ߿�����
        if cleft<1%��߿��������޴�ֱ�߿�������
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%�ұ߿��������޴�ֱ�߿�������
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
        SegBw2 = imresize(SegGray,[size_x size_y],'bilinear');%�任Ϊ40��*20�б�׼��ͼ ˫���Բ�ֵ        %2013aĬ�Ϸ�ʽ'bicubic'��˫���β�ֵ��
        [SegBw2]=pre_processing(SegBw2);
%         subcol2=imresize(subcol1,[size_x size_y],'bilinear');
%         subplot(3,n1,l),imshow(SegGray);%��ÿ��ԭ�ָ���ͼ��n1~=7
%         if l==4%������ʾ����
%             title(['�����ַ���ȣ� ',int2str(maxwidth)],'Color','r');
%         end
%        subplot(3,n1,n1+l),imshow(SegBw2);title(int2str(l),'Color','r');%��ÿ����׼��ͼ,n1~=7   
%         subcol2=im2bw(subcol2);
%         subplot(3,n1,2*n1+l),imshow(subcol2);
        file_name=strcat('..\�и���ͼ\�и���ͼ',int2str(l),'.jpg');%�洢��׼��ͼ
        imwrite(SegBw2,file_name,'jpg') %�洢��׼��ͼ
        l=l+1;
end
% SegGray1=sbw(:,cleft_most:cright_most);
% hiscol_SegGray1=sum(SegGray1);
% figure,subplot(2,1,1),bar(hiscol_SegGray1),title('ȥ��ֱ�߿�ֱͶӰ');
% subplot(2,1,2),imshow(SegGray1),title('ȥ��ֱ�߿�Ҷ�ͼ');
%=========================ȥ��ֱ�߿��ַ��ָ����====================================
waitbar(13/20,hbar,'��ȴ���ʶ����...'); 
%%
%=================================�ַ�ʶ��============================================
switch shibiefangfa
    case 1 %ģ��ʶ��
%����������ȡ���ַ�ͼ�������������ƥ�䣬�Զ�ʶ����ַ����롣
liccode=char(['0':'9' 'A':'H' 'J':'N' 'P':'Z' '����³��ԥ����']); %�����Զ�ʶ���ַ������  
length_code=length(liccode);
l=1;
[~,n2]=size(sbw);
RegCode=char(zeros(1,14));
%!figure,
for k=findmax-1:findmax+5%�ַ���ʼ
       cleft=markcol5(k)-maxwidth/2;%�ַ���߿�����
        cright=markcol5(k)+maxwidth/2-2;%�ַ��ұ߿�����
        if cleft<1%��߿��������޴�ֱ�߿�������
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%�ұ߿��������޴�ֱ�߿�������
            cright=n2;
            cleft=n2-maxwidth;
        end
        SegBw1=subcol(:,cleft:cright);
        SegBw2=imresize(SegBw1,[size_x size_y],'bilinear');%�任Ϊ40��*20�б�׼��ͼ
        [SegBw2]=pre_processing(SegBw2);
        SegBw2=im2bw(SegBw2);
%!        subplot(1,n1,k),imshow(SegBw2);
        if l==1                 %��һλ����ʶ��
            kmin=35;
            kmax=length_code;
        elseif l==2             %�ڶ�λ A~Z ��ĸʶ��
            kmin=11;
            kmax=34;
        elseif l>=3&&l<=7      %��3-7λ 0~9  A~Z��ĸ������ʶ��
            kmin=1;
            kmax=34;
        else                    %���填��λ 0~9 ����ʶ��
            kmin=1;
            kmax=10;
        end
        Differences=zeros(size(kmin:kmax));     
        for k2=kmin:kmax
            file_name=strcat('..\ģ����ͼ\',liccode(k2),'.jpg');
            SamBw2 = imresize(imread(file_name),[size_x size_y],'bilinear');
            SamBw2=im2bw(SamBw2);
            SubBw2 = SegBw2-SamBw2;%��ֹģ����ͼ�����ϴ�С,�����׼��ͼ��ģ����ͼ�Ĳ��
            Differences(k2-kmin+1)=sum(sum(SubBw2~=0,2),1);%ͳ�Ʋ�ķ�����������������
        end
        MinError=min(Differences);%ȡ������Сֵ
        findc=find(Differences==MinError);%������С����ͼ��
        RegCode(l)=liccode(findc(1)+kmin-1);%����ϵı�׼��ͼλ��
%         RegCode(l*2-1)=liccode(findc(1)+kmin-1);%����ϵı�׼��ͼλ��
%         RegCode(l*2)=' ';%�����С���ͼ��
        l=l+1;
end
% title (['ʶ���ƺ���:', RegCode],'Color','r');

    case 2%������ʶ��
        [RegCode]=bp_sim();
    case 3%ŷʽ����ʶ��
%����������ȡ���ַ�ͼ�������������ƥ�䣬�Զ�ʶ����ַ����롣
liccode=char(['0':'9' 'A':'Z' '����³��ԥ����']); %�����Զ�ʶ���ַ������  
length_code=length(liccode);
l=1;
[~,n2]=size(sbw);
RegCode=char(zeros(1,14));
%!figure,
for k=findmax-1:findmax+5%�ַ���ʼ
       cleft=markcol5(k)-maxwidth/2;%�ַ���߿�����
        cright=markcol5(k)+maxwidth/2-2;%�ַ��ұ߿�����
        if cleft<1%��߿��������޴�ֱ�߿�������
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%�ұ߿��������޴�ֱ�߿�������
            cright=n2;
            cleft=n2-maxwidth;
        end
        SegBw1=subcol(:,cleft:cright);
        SegBw2=imresize(SegBw1,[size_x size_y],'bilinear');%�任Ϊ40��*20�б�׼��ͼ
        [SegBw2]=pre_processing(SegBw2);
        SegBw2=im2bw(SegBw2);
%!        subplot(1,n1,k),imshow(SegBw2);
        if l==1                 %��һλ����ʶ��
            kmin=37;
            kmax=length_code;
        elseif l==2             %�ڶ�λ A~Z ��ĸʶ��
            kmin=11;
            kmax=36;
        elseif l>=3&&l<=7      %��3-7λ 0~9  A~Z��ĸ������ʶ��
            kmin=1;
            kmax=36;
        else                    %���填��λ 0~9 ����ʶ��
            kmin=1;
            kmax=10;
        end
%         Differences=zeros(size(kmin:kmax));
        ll=1;
        p=zeros(kmax-kmin,1);
        for k2=kmin:kmax
            file_name=strcat('..\ģ����ͼ\',liccode(k2),'.jpg');
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
        RegCode(l)=liccode(I+kmin-1);%����ϵı�׼��ͼλ��
        l=l+1;
end
        
        
        
end


%����ʶ�����
fid=fopen('..\ʶ���ƺ�.txt','w');
fprintf(fid,'%s\n',RegCode);
fclose(fid);
%�����Ի�����ʾʶ�����
waitbar(14/20,hbar,'��ȴ���ʶ����...');    
waitbar(1,hbar,'�����');
delete(hbar);
msgbox(['���ƺţ�',RegCode],'Author','');
% ��ȡ����
% if shibiefangfa~=2 %�����粻����
% for k=1:7
%     file_name=strcat('..\�����ļ�\',RegCode(k),'.wav');%���·��
%     wavplay(wavread(file_name),44100);
% end
% end
end


function [ ] = blue_recog( )
% clear;
% close all;
%1.732s Ĭ��
global im
global erzhi_yuzhi
global shibiefangfa
global tezheng
global RegCode
global shuipingxuanzhuan
global time_weizhi
%from kevin
%Ԥ���� ���� ��ɫʶ��
%ˮƽ�任 hough
%�����ַ�ϸ��ȥ��ģʽ
%����ȥí������
%����ˮƽ�߿����
%������RGB+Sͼ����Ϊ������
erzhi_yuzhi
shibiefangfa
shuipingxuanzhuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot_height=2048;  %ԭͼ��С
shot_width=1536;   %ԭͼ��С
parameter_s=14;    %strel����
rectangle_x=20;     %�˲����ο��С
rectangle_y=40;
rectangle_width=[shot_height*0.10,shot_height*0.45];%��ͨ���򣨳��ƣ������С��Χ
rectangle_hight=[shot_width*0.10,shot_width*0.3];
pr_rato=[2,10];%��ͨ���򣨳��ƣ��������
supplement_x=10;%����ͨ�������ƾ�����ȡ�Ĳ��䣨�ܷ�����Ӧ��
supplement_y=15;
% diff_ydata=2;%������תУ��������������ʹԤ������߹⻬����ֹ���䣨���ص㣩
size_x=40;%��ͼ��׼��Сx
size_y=20;%��ͼ��׼��Сy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=waitbar(0,'���Եȣ����ڿ�ʼʶ���ļ�...');

%=================ͼ�������Ԥ����=============================
rgb_image=im;
gray_image=rgb2gray(rgb_image);%ԭͼת�Ҷ�ͼ
[shot_width,shot_height]=size(gray_image);  %ԭͼ��С
ratio=shot_width/shot_height;
rgb_image=imresize(rgb_image,[1024*ratio,1024]);
% !figure,imshow(rgb_image),title('ԭ��ɫͼ��');%��ʾԭͼ��
gray_image=rgb2gray(rgb_image);%ԭͼת�Ҷ�ͼ
% !figure,imshow(gray_image),title('ԭ�Ҷ�ͼ');%��ʾ�Ҷ�ͼ
[shot_width,shot_height]=size(gray_image);  %ԭͼ��С
rectangle_width=[shot_height*0.10,shot_height*0.45];%��ͨ���򣨳��ƣ������С��Χ
rectangle_hight=[shot_width*0.010,shot_width*0.3];
R=rgb_image(:,:,1);
G=rgb_image(:,:,2);
B=rgb_image(:,:,3);
RW=roicolor(R,10,70);%ɸѡ��ɫ����
GW=roicolor(G,30,140);%ɸѡ��ɫ����
BW=roicolor(B,90,240);%ɸѡ��ɫ����
TOTLE_BW=RW+GW+BW;
TOTLE_BW=floor(TOTLE_BW/3);
hsv=rgb2hsv(rgb_image);
S=hsv(:,:,2);
% yuzhi_S=imhist(S,100);
SW=roicolor(S,0.55,1);%ɸѡ���Ͷȷ�������
SW=floor((SW+TOTLE_BW)/2);
bw2=SW;









% %�����ǿͼ��ȥ��������
% s=strel('disk',parameter_s);%!!!!!!!strel�����޸�
% Bgray_image=imopen(gray_image,s);%���������
% % !figure,imshow(Bgray_image),title('����ͼ��')%��ʾ����ͼ
% Egray_image=imsubtract(gray_image,Bgray_image);%��ͼ�������ǿͼ
% % !figure,imshow(Egray_image),title('��ǿͼ��')%��ʾ��ǿͼ
% waitbar(1/20,h,'��ȴ���ʶ����...');
% 
% %��ö�ֵ����ֵ
% image_max=double(max(max(Egray_image)));%��ǿͼ���ֵ�����˫������
% image_min=double(min(min(Egray_image)));%��ǿͼ��Сֵ�����˫������
% level1=(image_max-(image_max-image_min)/3)/255;%�����ֵ�޸�%��������ֵ
% level2 = graythresh(gray_image);%OSTU������ֵ
% switch erzhi_yuzhi
%     case 1
%         level=level1;
%     case 2
%         level=level2;
%     case 3
%         level=(level1+level2)/2;!!!ȡƽ��ֵ
% end
% bw2=im2bw(Egray_image,level);%ת��ͼ��Ϊ������ͼ��
% bw2=double(bw2);%ת��Ϊ˫������
% % !figure,imshow(bw2),title('��ֵͼ');
waitbar(2/20,h,'��ȴ���ʶ����...');

%�Զ�ֵͼ�������ղ��������˲�(��ʴ����)
grd=edge(bw2,'canny');%��canny����ʶ��ǿ��ͼ���еı߽�
%figure,imshow(grd);title('ͼ���Ե��ȡ');%���ͼ���Ե
bg1=imclose(grd,strel('rectangle',[rectangle_x,rectangle_y]));!!!!!!���������޸�%ȡ���ο�ı�����
%figure,imshow(bg1);title(['ͼ�������[',num2str(rectangle_x),',',num2str(rectangle_y),']']);%����������ͼ��
bg3=imopen(bg1,strel('rectangle',[rectangle_x,rectangle_y]));%ȡ���ο�Ŀ�����
%figure,imshow(bg3);title(['ͼ������[',num2str(rectangle_x),',',num2str(rectangle_y),']']);%����������ͼ��
bg2=imopen(bg3,strel('rectangle',[rectangle_x,1]));%ȡ���ο�Ŀ�����
%figure,imshow(bg2);title(['ͼ������[',num2str(rectangle_y),',1]']);%����������ͼ��
waitbar(3/20,h,'��ȴ���ʶ����...');

%�Զ�ֵͼ�����������ȡ���������������������������������������Ƚϣ���ȡ��������
[L,num] = bwlabel(bg2,8);%��ע������ͼ���������ӵĲ��֣�LΪ��ͨ�������numΪ��ͨ������Ŀ��
Feastats = regionprops(L,'basic');%����ͼ�������ϵ�������ߴ磨'Area'�Ǳ������������ͼ����������������ܸ�����'BoundingBox'��1��ndims(L)*2�е���������������Ӧ�������С���Ρ�BoundingBox ��ʽΪ [ul_corner width]������ ul_corner �� [x y z ...] ��������ʽ�����߽���ӵ����Ͻǡ�boxwidth �� [x_width y_width ...] ��ʽָ���߽��������ÿ��ά������ĳ��ȡ���
                                  %'Centroid'��1��ndims(L)�е�����������ÿ����������ģ����ģ���
Area=[Feastats.Area];%�������
BoundingBox=[Feastats.BoundingBox];%[x y width height]���ƵĿ�ܴ�С
RGB_image2= label2rgb(L, 'spring', 'k', 'noshuffle'); %��־ͼ����RGBͼ��ת��
file_name1=strcat('..\���Ʊ�׼ͼ\���Ʊ��ͼ','.jpg');%�洢��׼��ͼ
imwrite(RGB_image2,file_name1,'jpg') %�洢��׼��ͼ
% !figure,imshow(RGB_image2);title('ͼ���ɫ���');%�����ܵĲ�ɫͼ��
waitbar(4/20,h,'��ȴ���ʶ����...');

%������ɫͳ�Ƴ���ɸѡ
l=1;
ll=1;
num1=num;
% figure
while (l<=num1)
    startcol=floor(BoundingBox((l-1)*4+1));%y
    startrow=floor(BoundingBox((l-1)*4+2));%x
    if startcol==0
        startcol=1;
    end
    if startrow==0
        startrow=1;
    end
    width=BoundingBox((l-1)*4+3);%x��
    hight=BoundingBox((l-1)*4+4);%y��
    bw_sub=TOTLE_BW(startrow:startrow+hight,startcol:startcol+width);
%     subplot(num,1,ll),imshow(bw_sub);
    sum(sum(bw_sub))
    if sum(sum(bw_sub))<=0.01*Area(l)%��ɫ��������Ƿ�����������1%
        BoundingBox((l-1)*4+1:(l-1)*4+4)=[];
        l=l-1;
        num1=num1-1;
    end
    l=l+1;
    ll=ll+1;
end


if num1~=1%����ʣһ���������������ж�
%������ͨ�����Ƿ���ϳ��ƴ�С
lx=0;
Getok=zeros(1,num1);
for l=1:num1
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
else 
    startcol=BoundingBox(1)-supplement_x;%��ͨ�������Ͻ�����x
    startrow=BoundingBox(2)-supplement_y;%��ͨ�������Ͻ�����y
    width=BoundingBox(3)+supplement_x;%���ƿ�
    hight=BoundingBox(4)+supplement_y;%���Ƹ�
end

endrow=startrow+hight+supplement_y;
endcol=startcol+width+supplement_x;
R_sub=R(startrow:endrow,startcol:endcol);%ȡ��ͼ��ͼ
G_sub=G(startrow:endrow,startcol:endcol);%ȡ��ͼ��ͼ
B_sub=B(startrow:endrow,startcol:endcol);%ȡ��ͼ��ͼ
rgb_sub1=rgb_image(startrow:endrow,startcol:endcol);%ȡ��ͼ��ͼ
subcol1=gray_image(startrow:endrow,startcol:endcol);%��ȡ���ƻҶ���ͼ
sbw1=im2bw(subcol1);%��ȡ���ƶ�ֵ��ͼ
grd1=edge(sbw1,'canny');%��ȡ���Ʊ߽���ͼ
% figure,subplot(3,1,1),imshow(subcol1);title('���ƻҶ���ͼ');%��ʾ�Ҷ�ͼ��
% subplot(3,1,2),imshow(sbw1);title('���ƶ�ֵ��ͼ');%��ʾ���ƵĶ�ֵͼ
% subplot(3,1,3),imshow(grd1);title('���Ʊ߽���ͼ');%��ʾ���Ƶı߽�
%===============================Ԥ�������=======================================
waitbar(5/20,h,'��ȴ���ʶ����...');
%%
%==============================������תУ��======================================
if shuipingxuanzhuan==1
    %==================ˮƽУ����ת=================

binaryImage_horizontal = edge(subcol1,'sobel','horizontal');
% subplot(4,1,4),imshow(binaryImage_horizontal),title('ˮƽ�߽�ͼ')
theta=0:179;
r=radon(binaryImage_horizontal,theta);
[m,n]=size(r);%m=�Ƕȣ�n=λ��
[I,temp] = find(r>=max(max(r)));%temp��¼����б�ǣ�������б��
angle=90-temp;

R_sub = imrotate(R_sub,angle); %��ת��ɫ��ͼ
G_sub = imrotate(G_sub,angle); %��ת��ɫ��ͼ
B_sub = imrotate(B_sub,angle); %��ת��ɫ��ͼ
rgb_sub(:,:,1)=R_sub;rgb_sub(:,:,2)=G_sub;rgb_sub(:,:,3)=B_sub;%�ϳɲ�ɫ��ͼ
subcol = imrotate(subcol1,angle,'bilinear','crop'); %��ת�Ҷ���ͼ
sbw = imrotate(sbw1,angle,'bilinear','crop');%��ת��ֵ��ͼ
grd3 = edge(subcol,'canny');%���¼����ͼ�߽� ��֤�������ת����
grd4 = bwmorph(grd3,'thicken'); %ͨ����Ŀ���ⲿ�������ؼӺ�Ŀ��ֱ������������ʹ��ǰδ����Ŀ���Ϊ8��ͨ�򡣵õ�ǿ���߽�
[hight,width]=size(sbw);
% figure,subplot(5,1,1),imshow(subcol);title('���ƻҶ���ͼ(ˮƽУ����ת��)');%���������ת��ĻҶ�ͼ�������ʾ���ƻҶ���ͼ
% subplot(5,1,2),imshow(sbw);title('���ƶ�ֵ��ͼ��ˮƽУ����ת��')%���������ת��ĻҶ�ͼ��
% subplot(5,1,3),imshow(grd3);title('���Ʊ߽���ͼ��ˮƽУ����ת��')%���������ת��ı߽���ͼ
% subplot(5,1,4),imshow(grd4);%���������ת���ǿ���߽���ͼ
% subplot(5,1,5),imshow(rgb_sub);%���������ת���ǿ���߽���ͼ
% title(['����ˮƽУ����ת��: ',num2str(angle),'��'] ,'Color','r');%��ʾ���Ƶ���ת�Ƕ�


%===============ˮƽ��תУ������=============
elseif shuipingxuanzhuan==2
%===============================͸����ת=========================================
% file_name=strcat('.\���Ʊ�׼ͼ\�Ҷȳ���','.jpg');%�洢��׼��ͼ
% img= imread(file_name);
img=subcol1;
tiaobian=10;
[M,N] = size(img);
allow_diff_y=M/40;%y����������½��߶�M/10
allow_diff_x=N/40;%x����������½��߶�N/20
bw= im2bw(img);
% figure,imshow(mat2gray(img))
% figure,imshow(bw);
hold on,plot([N/10,N/10],[1,M]), %���ٱ߽ǵ�Ӱ��
        plot([9*N/10,9*N/10],[1,M]),%���ٱ߽ǵ�Ӱ��
        plot([1,N],[M/5,M/5]), %���ٱ߽ǵ�Ӱ��
        plot([1,N],[4*M/5,4*M/5]), %���ٱ߽ǵ�Ӱ��
hold off
%����ϱ���
l=1;
 for j=ceil(N/10):tiaobian:ceil(9*N/10)
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
 z=polyval(fresult1,xdata);
 1
 figure,
 imshow(bw),hold on, plot(xdata',ydata','r',xdata,z,'b'),
 %����±���
%����ϱ���
clear xdata ydata
l=1;
 for j=ceil(N/10):tiaobian:ceil(9*N/10)
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
 z=polyval(fresult2,xdata);
hold on, plot(xdata',ydata','r',xdata,z,'b'),
2
%��������
clear xdata ydata
l=1;
 for i=ceil(M/5):tiaobian/2:ceil(4*M/5)
     for j=1:N
         if bw(i,j)==1
             if l==1
             xdata(l)=j;
             diff_data=0;
             else
             diff_data=abs(j-xdata(l-1));
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
 z=polyval(fresult3,xdata);
 hold on, plot(xdata',ydata','r',xdata,z,'b'),
 3
 %����ұ���
clear xdata ydata
l=1;
 for i=ceil(M/5):tiaobian/2:ceil(4*M/5)
     for j=N:-1:1
         if bw(i,j)==1
             if l==1
             xdata(l)=j;
             diff_data=0;
             else
             diff_data=abs(j-xdata(l-1));
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
 z=polyval(fresult4,xdata);
 hold on, plot(xdata',ydata','r',xdata,z,'b'),hold off
 4
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
imgn=zeros(height,width);

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
figure,subplot(2,1,1),imshow(subcol),title('͸�ӱ任�Ҷ�ͼ');
       subplot(2,1,2),imshow(sbw),title('͸�ӱ任��ֵͼ');
else
    rgb_sub=rgb_sub1;
    sbw=sbw1;
    subcol=subcol1;
end

% binaryImage_horizontal = edge(subcol1,'sobel','horizontal');
% % !subplot(4,1,4),imshow(binaryImage_horizontal),title('ˮƽ�߽�ͼ')
% theta=0:179;
% r=radon(binaryImage_horizontal,theta);
% % [m,n]=size(r);%m=�Ƕȣ�n=λ��
% [~,temp] = find(r>=max(max(r)));%temp��¼����б�ǣ�������б��
% angle=90-temp;
% subcol = imrotate(subcol1,angle,'bilinear','crop'); %��ת�Ҷ���ͼ
% sbw = imrotate(sbw1,angle,'bilinear','crop');%��ת��ֵ��ͼ
% % grd3 = edge(subcol,'canny');%���¼����ͼ�߽� ��֤�������ת����
% % grd4 = bwmorph(grd3,'thicken'); %ͨ����Ŀ���ⲿ�������ؼӺ�Ŀ��ֱ������������ʹ��ǰδ����Ŀ���Ϊ8��ͨ�򡣵õ�ǿ���߽�
% [hight,width]=size(sbw);
%!figure,subplot(4,1,1),imshow(subcol);title('���ƻҶ���ͼ(ˮƽУ����ת��)');%���������ת��ĻҶ�ͼ�������ʾ���ƻҶ���ͼ
% !subplot(4,1,2),imshow(sbw);title('���ƶ�ֵ��ͼ��ˮƽУ����ת��')%���������ת��ĻҶ�ͼ��
% !subplot(4,1,3),imshow(grd3);title('���Ʊ߽���ͼ��ˮƽУ����ת��')%���������ת��ı߽���ͼ
% !subplot(4,1,4),imshow(grd4);%���������ת���ǿ���߽���ͼ
% !title(['����ˮƽУ����ת��: ',num2str(angle),'��'] ,'Color','r');%��ʾ���Ƶ���ת�Ƕ�
%===============ˮƽ��תУ������=============
file_name1=strcat('..\���Ʊ�׼ͼ\�Ҷȳ���','.jpg');%�洢��׼��ͼ
imwrite(subcol,file_name1,'jpg') %�洢��׼��ͼ
waitbar(7/20,h,'��ȴ���ʶ����...');

%=========================����ȥƽ�п�=================================
%��ת���ƺ����¼��㳵��ˮƽͶӰ��ȥ������ˮƽ�߿򣬻�ȡ�ַ��߶�
[hight,width]=size(sbw);
histcol1=sum(sbw); %���㴹ֱͶӰ
histrow=sum(sbw'); %����ˮƽͶӰ
% !figure,subplot(2,1,1),bar(histcol1),title('��ֱͶӰ����ת��');
% !subplot(2,1,2),bar(histrow),title('ˮƽͶӰ����ת��');


% ��ˮƽͶӰ���з�ȷ���
meanrow=mean(histrow);%��ˮƽͶӰ��ƽ��ֵ
minrow=min(histrow);%��ˮƽͶӰ����Сֵ
levelrow=(meanrow+minrow)/1.82;%!!!!!!!!!!���޸�%��ˮƽͶӰde���岨����ֵ
% !figure,subplot(2,1,1),bar(histrow),hold on,plot([1,hight],[levelrow,levelrow]),hold off;title('ˮƽͶӰ�����߿�');%���ˮƽͶӰ�Ͷ�ֵͼ
% !subplot(2,1,2),imshow(sbw),title('���ƶ�ֵ��ͼ');

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
%l=0;
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
%!figure,subplot(2,1,1),imshow(sbw2); title(['�����ַ��߶ȣ� ',int2str(maxhight)],'Color','r');%��������ַ��߶�
%========================ȥˮƽ�߿����===========================
waitbar(9/20,h,'��ȴ���ʶ����...');  
%%
%==============��ֱУ����ת=================
% grd4 = edge(subcol,'canny'); 
% grd4 = bwmorph(grd4,'thicken'); 
% % !subplot(2,1,2),imshow(grd4),title('ǿ���߽�ȥˮƽ�߿�ͼ')%�������ǿ���߽�
% theta = -90:89;
% [R,~] = radon(grd4,theta);%R=�Ƕȣ�xp=λ��
% [R1,~] = max(R);
% theta_max = 90;
% while (theta_max > 20 || theta_max< -20)
%     [~,theta_max] = max(R1);                      
%     R1(theta_max) = 0; 
% %     theta_max = theta_max - 91;
% end
% angle_x=tan(theta_max);
% H=[1,0,0; angle_x,1,0;0,0,1];
% T=maketform('affine',H);
% subcol=imtransform(subcol,T,'bilinear');
% % sbw=imtransform(sbw,T,'bilinear');
% % !figure,subplot(2,1,1),imshow(subcol);title('���ƻҶ���ͼ(��ֱУ����ת��)');%���������ת��ĻҶ�ͼ�������ʾ���ƻҶ���ͼ
% % !subplot(2,1,2),imshow(sbw);title(['���ƶ�ֵ��ͼ����ֱУ����ת��',int2str(angle_x),'��'],'color','r')%���������ת��ĻҶ�ͼ��
%================��ֱУ����ת����============
waitbar(11/20,h,'��ȴ���ʶ����...');  

%====================ȥ��ֱ�߿��ַ��ָ�======================
sbw=im2bw(subcol);
histcol=sum(sbw);  %���㴹ֱͶӰ
meancol=mean(histcol);%��ֱͶӰ��ƽ��ֵ
mincol=min(histcol);%��ֱͶӰ����Сֵ
levelcol=(meancol+mincol)/4;%!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/4
levelcol1=(meancol+mincol)/2.75;%!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/2.75
levelcol2=meancol;
%!figure,bar(histcol);hold on;plot([1,width],[levelcol,levelcol]);hold off;title('��ֱͶӰ��ȥˮƽ�߿��');%������ƵĴ�ֱͶӰͼ��

%���㳵�ƴ�ֱͶӰ�������㣬�ȿ�ȣ�����룬�½��㣬���ȣ���λ�ã�������
% meancol=mean(histcol);%��ֱͶӰ��ƽ��ֵ
% mincol=min(histcol);%��ֱͶӰ��ƽ��ֵ
% levelcol=(meancol+mincol)/2;!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/4
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
rem_position=find(markcol7<=levelcol2);
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
% cleft_most=[];
% cright_most=[];
%!figure;
for k=findmax-1:findmax+5%ȡ�����󴦼�í�������ڶ����ַ�����ǰһ�����������ȥ����ֱ�߿�
        cleft=markcol5(k)-maxwidth/2;%�ַ���߿�����
        cright=markcol5(k)+maxwidth/2;%�ַ��ұ߿�����
        if cleft<1%��߿��������޴�ֱ�߿�������
            cleft=1;
            cright=maxwidth;
        end
        if cright>n2%�ұ߿��������޴�ֱ�߿�������
            cright=n2;
            cleft=n2-maxwidth;
        end
        if l==1
            cleft_most=cleft;
        end
        if l==7
            cright_most=cright;
        end
        SegGray=subcol(:,cleft:cright);
        [SegBw2]=pre_processing(SegGray);
        [L,num] = bwlabel(SegBw2,8);%��ע������ͼ���������ӵĲ��֣�LΪ��ͨ�������numΪ��ͨ������Ŀ��
        Feastats = regionprops(L,'basic');%����ͼ�������ϵ�������ߴ磨'Area'�Ǳ������������ͼ����������������ܸ�����'BoundingBox'��1��ndims(L)*2�е���������������Ӧ�������С���Ρ�BoundingBox ��ʽΪ [ul_corner width]������ ul_corner �� [x y z ...] ��������ʽ�����߽���ӵ����Ͻǡ�boxwidth �� [x_width y_width ...] ��ʽָ���߽��������ÿ��ά������ĳ��ȡ���
                                         %'Centroid'��1��ndims(L)�е�����������ÿ����������ģ����ģ���
        Area=[Feastats.Area];%�������
        most_area=Area;
        BoundingBox=round([Feastats.BoundingBox]);%[x y width height]���ƵĿ�ܴ�С
%         RGB_image2= label2rgb(L, 'spring', 'k', 'noshuffle'); %��־ͼ����RGBͼ��ת��
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %ȥ���߽紦΢С��ͨ��
         if k==findmax-1   %���ݺ�����Ӣ���������𣬺�������˵������С������Ӣ����ͨ������ϴ󣬿ɹ��������
             refer_area=40;
         else
             refer_area=110;
         end
         BoundingBox1=BoundingBox;
         BoundingBox1(4:4:4+4*(num-1))=BoundingBox1(2:4:2+4*(num-1))+BoundingBox1(4:4:4+4*(num-1))-1; %���½Ǿ���x
         BoundingBox1(3:4:3+4*(num-1))=BoundingBox1(1:4:1+4*(num-1))+BoundingBox1(3:4:3+4*(num-1))-1;  %���½Ǿ���y
         BoundingBox2=BoundingBox1;
         count=num;
         i=1;
         while i<=count
             if BoundingBox2(4*(i-1)+1)<=1||BoundingBox2(4*(i-1)+2)<=1||BoundingBox2(4*(i-1)+4)>=40||BoundingBox2(4*(i-1)+3)>=20
                 if Area(i)<refer_area
                     BoundingBox2(4*(i-1)+1:4*(i-1)+4)=[];%������40��λ�ڱ߽����ͨ��
                     Area(i)=[];
                     count=count-1;
                     i=i-1;
                 end
             end
             i=i+1;
         end
         if count==0%��ֹʧ��ɾ������
             i=find(most_area==max(most_area));
             BoundingBox2(1:4)=BoundingBox1(1+4*(i-1):4+4*(i-1));
             count=1;
         end
             Bounding_y=min(BoundingBox2(1:4:1+4*(count-1)));
             Bounding_x=min(BoundingBox2(2:4:2+4*(count-1)));
             Bounding_yy=max(BoundingBox2(3:4:3+4*(count-1)));
             Bounding_xx=max(BoundingBox2(4:4:4+4*(count-1)));
          SegBw2=SegBw2(Bounding_x:Bounding_xx,Bounding_y:Bounding_yy);
          [SegBw2]=pre_processing(SegBw2);
%           subplot(5,n1,4*n1+l),imshow(SegBw2); 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subcol1=subcol(:,cleft:cright);
%         SegBw2 = imresize(SegGray,[size_x size_y],'bilinear');%�任Ϊ40��*20�б�׼��ͼ ˫���Բ�ֵ        %2013aĬ�Ϸ�ʽ'bicubic'��˫���β�ֵ��
%         [SegBw2]=pre_processing(SegBw2);
%!        subcol2=imresize(subcol1,[size_x size_y],'bilinear');
%!        subplot(3,n1,l),imshow(SegGray);%��ÿ��ԭ�ָ���ͼ��n1~=7
%!        if l==4%������ʾ����
%!            title(['�����ַ���ȣ� ',int2str(maxwidth)],'Color','r');
%!        end
%!        subplot(3,n1,n1+l),imshow(SegBw2);title(int2str(l),'Color','r');%��ÿ����׼��ͼ,n1~=7   
%!        subcol2=im2bw(subcol2);
%!        subplot(3,n1,2*n1+l),imshow(subcol2);
        file_name=strcat('..\�и���ͼ\�и���ͼ',int2str(l),'.jpg');%�洢��׼��ͼ
        imwrite(SegBw2,file_name,'jpg') %�洢��׼��ͼ
        l=l+1;
end
%handles = guihandles(gcf);
%hiscol_SegGray1=sum(SegGray1);
%!figure,subplot(2,1,1),bar(hiscol_SegGray1),title('ȥ��ֱ�߿�ֱͶӰ');
%!subplot(2,1,2),imshow(SegGray1),title('ȥ��ֱ�߿�Ҷ�ͼ');
%=========================ȥ��ֱ�߿��ַ��ָ����====================================
file_name1=strcat('..\���Ʊ�׼ͼ\��ֵ����','.jpg');%�洢��׼��ͼ
if isempty(cleft_most)||isempty(cright_most)
    msgbox('ʶ�����','Author','');
    return
end
SegGray1=sbw(:,cleft_most:cright_most);
imwrite(SegGray1,file_name1,'jpg') %�洢��׼��ͼ
waitbar(13/20,h,'��ȴ���ʶ����...');  
%%
%=================================�ַ�ʶ��============================================
switch shibiefangfa
    case 1 %ģ��ʶ��
%����������ȡ���ַ�ͼ�������������ƥ�䣬�Զ�ʶ����ַ����롣
liccode=char(['0':'9' 'A':'H' 'J':'N' 'P':'Z' '����³��ԥ�����ڶ���������']); %�����Զ�ʶ���ַ������  
length_code=length(liccode);
l=1;
[~,n2]=size(sbw);
RegCode=char(zeros(1,14));
%!figure,
for k=findmax-1:findmax+5%�ַ���ʼ
       file_name=strcat('..\�и���ͼ\�и���ͼ',int2str(l),'.jpg');%�洢��׼��ͼ
       SegBw2=imread(file_name); %�洢��׼��ͼ
       SegBw2=im2bw(SegBw2);
%        subplot(1,n1,k),imshow(SegBw2);
        if l==1                 %��һλ����ʶ��
            kmin=35;
            kmax=length_code;
        elseif l==2             %�ڶ�λ A~Z ��ĸʶ��
            kmin=11;
            kmax=34;
        elseif l>=3&&l<=7      %��3-7λ 0~9  A~Z��ĸ������ʶ��
            kmin=1;
            kmax=34;
        else                    
            kmin=1;
            kmax=10;
        end
        Differences=zeros(size(kmin:kmax));
        
%                 for k2=kmin:kmax
%             fname=strcat('D:\users\Desktop\SRTP\����\kevin\ģ����ͼ\',liccode(k2),'.jpg');
%             SamBw2 = imread(fname);  
%             SamBw2=imresize(SamBw2,[size_x size_y]);%�任Ϊ40��*20�б�׼��ͼ
%             SubBw2 = uint8(SegBw2)-SamBw2;
%             Dmax=sum(sum(SubBw2~=0,2),1);
%             Error(k2)=Dmax;
%                 end
%         Error1=Error(kmin:kmax);%�Ƚ����
%         MinError=min(Error1);%ȡ������Сֵ
%         findc=find(Error1==MinError);%������С����ͼ��
        
        for k2=kmin:kmax
            file_name=strcat('..\ģ����ͼ\',liccode(k2),'.jpg');
%             SamBw2 =imread(file_name);
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
liccode=char(['0':'9' 'A':'H' 'J':'N' 'P':'Z' '����³��ԥ�����ڶ���������']); %�����Զ�ʶ���ַ������  
length_code=length(liccode);
l=1;
[~,n2]=size(sbw);
RegCode=char(zeros(1,14));
%!figure,
for k=findmax-1:findmax+5%�ַ���ʼ
      file_name=strcat('..\�и���ͼ\�и���ͼ',int2str(l),'.jpg');%�洢��׼��ͼ
       SegBw2=imread(file_name); %�洢��׼��ͼ
       SegBw2=im2bw(SegBw2);
%        subplot(1,n1,k),imshow(SegBw2);
        if l==1                 %��һλ����ʶ��
            kmin=35;
            kmax=length_code;
        elseif l==2             %�ڶ�λ A~Z ��ĸʶ��
            kmin=11;
            kmax=34;
        elseif l>=3&&l<=7      %��3-7λ 0~9  A~Z��ĸ������ʶ��
            kmin=1;
            kmax=34;
        else                    
            kmin=1;
            kmax=10;
        end
%         Differences=zeros(size(kmin:kmax));
%         Differences=zeros(size(kmin:kmax));
        ll=1;
        p=zeros(kmax-kmin,1);
        for k2=kmin:kmax
            file_name=strcat('..\ģ����ͼ\',liccode(k2),'.jpg');
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
        RegCode(l)=liccode(I+kmin-1);%����ϵı�׼��ͼλ��

        l=l+1;
end
        
        
        
end
RegCode
%�����Ի�����ʾʶ�����
waitbar(14/20,h,'��ȴ���ʶ����...');    
waitbar(1,h,'�����');
delete(h);
msgbox(['���ƺţ�',RegCode],'Author','');
% ��ȡ����
if shibiefangfa~=2 %�����粻����
for k=1:7
    file_name=strcat('..\�����ļ�\',RegCode(k),'.wav');%���·��
    wavplay(wavread(file_name),44100);
end
end

%����ʶ�����
RegCode(1)=liccode(time_weizhi+1);%%%%%λ���޸ķ�
fid=fopen('..\ʶ���ƺ�.txt','w');
fprintf(fid,'%s\n',RegCode);
fclose(fid);

end


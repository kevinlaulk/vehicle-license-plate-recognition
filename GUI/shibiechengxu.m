function [ ] = shibiechengxu( )
% clear;
% close all;
%1.732s Ĭ��
global im
global erzhi_yuzhi
global shibiefangfa
global tezheng
global RegCode
%from keviglobal im
erzhi_yuzhi
shibiefangfa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_s=14;    %strel����
rectangle_x=5;     %�˲����ο��С
rectangle_y=50;
rectangle_width=[98,400];%��ͨ���򣨳��ƣ������С
rectangle_hight=[25,100];
pr_rato=[2,10];%��ͨ���򣨳��ƣ��������
supplement_x=5;%����ͨ�������ƾ�����ȡ�Ĳ��䣨�ܷ�����Ӧ��
supplement_y=15;
%diff_ydata=2;%������תУ��������������ʹԤ������߹⻬����ֹ���䣨���ص㣩
size_x=40;%��ͼ��׼��Сx
size_y=20;%��ͼ��׼��Сy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=waitbar(0,'���Եȣ����ڿ�ʼʶ���ļ�...');

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
waitbar(1/20,h,'��ȴ���ʶ����...');

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
        level=(level1+level2)/2;!!!ȡƽ��ֵ
end
bw2=im2bw(Egray_image,level);%ת��ͼ��Ϊ������ͼ��
bw2=double(bw2);%ת��Ϊ˫������
% !figure,imshow(bw2),title('��ֵͼ');
waitbar(2/20,h,'��ȴ���ʶ����...');

%�Զ�ֵͼ�������ղ��������˲�(��ʴ����)
grd=edge(bw2,'canny');%��canny����ʶ��ǿ��ͼ���еı߽�
% !figure,imshow(grd);title('ͼ���Ե��ȡ');%���ͼ���Ե
bg1=imclose(grd,strel('rectangle',[rectangle_x,rectangle_y]));%!!!!!!���������޸�%ȡ���ο�ı�����
% !figure,imshow(bg1);title('ͼ�������[5,19]');%����������ͼ��
bg3=imopen(bg1,strel('rectangle',[rectangle_x,rectangle_y]));%ȡ���ο�Ŀ�����
%!figure,imshow(bg3);title('ͼ������[5,19]');%����������ͼ��
bg2=imopen(bg3,strel('rectangle',[rectangle_y,1]));%ȡ���ο�Ŀ�����
%!figure,imshow(bg2);title('ͼ������[19,1]');%����������ͼ��
waitbar(3/20,h,'��ȴ���ʶ����...');

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
waitbar(4/20,h,'��ȴ���ʶ����...');

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
    rato=width/hight;%���㳵�Ƴ����
    if rato>pr_rato(1) && rato<pr_rato(2)   
        break;
    end
end
sbw1=bw2(startrow:startrow+hight,startcol:startcol+width+5); %��ȡ���ƶ�ֵ��ͼ
subcol1=gray_image(startrow:startrow+hight,startcol:startcol+width+5);%��ȡ���ƻҶ���ͼ
%grd1=grd(startrow:startrow+hight,startcol:startcol+width+5);%��ȡ���Ʊ߽���ͼ
% !figure,subplot(2,1,1),imshow(subcol1);title('���ƻҶ���ͼ');%��ʾ�Ҷ�ͼ��
% !subplot(2,1,2),imshow(sbw1);title('���ƶ�ֵ��ͼ');%��ʾ���ƵĶ�ֵͼ
%===============================Ԥ�������=======================================
waitbar(5/20,h,'��ȴ���ʶ����...');

%==============================������תУ��======================================
%���㳵��ˮƽͶӰ������ˮƽͶӰ���з�ȷ���
% histcol1=sum(sbw1);      %���㴹ֱͶӰ
% histrow=sum(sbw1');      %����ˮƽͶӰ
% !figure,subplot(2,1,1),bar(histcol1),title('��ֱͶӰ�����߿�');%�����ֱͶӰ��ˮƽͶӰ
% !subplot(2,1,2),bar(histrow),title('ˮƽͶӰ�����߿�');
% !figure,subplot(2,1,1),bar(histrow),title('ˮƽͶӰ�����߿�');%���ˮƽͶӰ�Ͷ�ֵͼ
% !subplot(2,1,2),imshow(sbw1),title('���ƶ�ֵ��ͼ');
% %��ˮƽͶӰ���з�ȷ���
% meanrow=mean(histrow);%��ˮƽͶӰ��ƽ��ֵ
% minrow=min(histrow);%��ˮƽͶӰ����Сֵ
% levelrow=(meanrow+minrow)/2;!!!!!!!!!!���޸�%��ˮƽͶӰde���岨����ֵ
% waitbar(6/20,h,'��ȴ���ʶ����...');
% 
% %����ˮƽͶӰ�������㡢�½��㡢�����ȡ����ȿ��
% count1=0;
% l=1;
% % markrow=zeros(1,hight+1);
% % markrow1=zeros(1,hight+1);
% for k=1:hight
%     if histrow(k)<=levelrow   %���ȴ�                          
%         count1=count1+1;                                
%     else 
%         if count1>=1
%             markrow(l)=k;%������ֵ����㣨�õ��Ѵ�����ֵ)�������㣩
%             markrow1(l)=count1;%�ȿ�ȣ��½�������һ�������㣩
%             l=l+1;
%         end
%         count1=0;
%     end
% end
% markrow(l)=hight;%���һ�εĲ��Ȳ���
% markrow1(l)=count1;
% markrow2=diff(markrow);%����루����������һ�������㣩
% [~,n1]=size(markrow2);
% l=0;
% markrow3=markrow(2:n1+1)-markrow1(2:n1+1);%������ֵ���յ㣨�õ�С����ֵ�����½��㣩
% markrow4=markrow3-markrow(1:n1);%���ȣ����������½��㣩
% markrow5=markrow3-double(uint16(markrow4./2));%������λ�ã�ȡ����ת����ʽ��
% % for k=1:n1
% %     markrow3(k)=markrow(k+1)-markrow1(k+1);%������ֵ���յ㣨�õ�С����ֵ�����½��㣩
% %     markrow4(k)=markrow3(k)-markrow(k);%���ȣ����������½��㣩
% %     markrow5(k)=markrow3(k)-double(uint16(markrow4(k)/2));%������λ�ã�ȡ����ת����ʽ��
% % end
% waitbar(7/20,h,'��ȴ���ʶ����...');

% %���㳵����ת�Ƕ�
% %(1)�����������½����ҵ�һ��Ϊ1�ĵ�
% [m2,n2]=size(sbw1);%sbw1��ͼ���С
% [~,n1]=size(markrow4);%markrow4�Ĵ�С
% maxw=max(markrow4);%�����Ϊ�ַ�
% if markrow4(1) ~= maxw%�����ϱ߼���ϱ�
%     ysite=1;
%     k1=1;
%     for l=1:n2
%     for k=1:markrow3(ysite)%�Ӷ�������һ�����½���ɨ��
%         if sbw1(k,l)==1
%             if k1==1
%             xdata(k1)=l;%���ϱ߽�x����
%             ydata(k1)=k;%���ϱ߽�y����
%             k1=k1+1;
%             break;
%             elseif abs(ydata(k1-1)-k)<diff_ydata%����������ʹ��Ϲ⻬����ֹ����
%             xdata(k1)=l;%���ϱ߽�x����
%             ydata(k1)=k;%���ϱ߽�y����
%             k1=k1+1;
%             break;
%             end
%         end
%     end
%     end
% else  %���ϱ߼���±�
%     ysite=n1;
%     if markrow4(n1) ==0
%         if markrow4(n1-1) ==maxw
%            ysite= 0; %���±�
%        else
%            ysite= n1-1;
%        end
%     end
%     if ysite ~=0
%         k1=1;
%         for l=1:n2
%             k=m2;
%             while k>=markrow(ysite) %�ӵױ������һ�����������ɨ��
%                 if sbw1(k,l)==1
%                     if k1==1
%                     xdata(k1)=l;%���±߽�x����
%                     ydata(k1)=k;%���±߽�y����
%                     k1=k1+1;
%                     break;
%                     elseif abs(ydata(k1-1)-k)<diff_ydata%����������ʹ��Ϲ⻬����ֹ����
%                     xdata(k1)=l;%���±߽�x����
%                     ydata(k1)=k;%���±߽�y����
%                     k1=k1+1;   
%                     break;
%                     end
%                 end
%                 k=k-1;
%             end
%         end
%     end
% end       
% waitbar(8/20,h,'��ȴ���ʶ����...');
% 
% %������ϣ�������x�н�
% fresult=polyfit(xdata',ydata',1);%һ���������
% angle=atand(fresult(1))*2.5;!!!!!!���ƽ�������Ľ�%�����Ƕȣ��ȣ�
% 
% %��ת����ͼ��
%==================ˮƽУ����ת=================

binaryImage_horizontal = edge(subcol1,'sobel','horizontal');
% !subplot(4,1,4),imshow(binaryImage_horizontal),title('ˮƽ�߽�ͼ')
theta=0:179;
r=radon(binaryImage_horizontal,theta);
% [m,n]=size(r);%m=�Ƕȣ�n=λ��
[~,temp] = find(r>=max(max(r)));%temp��¼����б�ǣ�������б��
angle=90-temp;
subcol = imrotate(subcol1,angle,'bilinear','crop'); %��ת�Ҷ���ͼ
sbw = imrotate(sbw1,angle,'bilinear','crop');%��ת��ֵ��ͼ
% grd3 = edge(subcol,'canny');%���¼����ͼ�߽� ��֤�������ת����
% grd4 = bwmorph(grd3,'thicken'); %ͨ����Ŀ���ⲿ�������ؼӺ�Ŀ��ֱ������������ʹ��ǰδ����Ŀ���Ϊ8��ͨ�򡣵õ�ǿ���߽�
[hight,width]=size(sbw);
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
% histcol1=sum(sbw); %���㴹ֱͶӰ
histrow=sum(sbw'); %����ˮƽͶӰ
% !figure,subplot(2,1,1),bar(histcol1),title('��ֱͶӰ����ת��');
% !subplot(2,1,2),bar(histrow),title('ˮƽͶӰ����ת��');


% ��ˮƽͶӰ���з�ȷ���
meanrow=mean(histrow);%��ˮƽͶӰ��ƽ��ֵ
minrow=min(histrow);%��ˮƽͶӰ����Сֵ
levelrow=(meanrow+minrow)/2;%!!!!!!!!!!���޸�%��ˮƽͶӰde���岨����ֵ
% !figure,subplot(2,1,1),bar(histrow),hold on,plot([1,hight],[levelrow,levelrow]),hold off;title('ˮƽͶӰ�����߿�');%���ˮƽͶӰ�Ͷ�ֵͼ
% !subplot(2,1,2),imshow(sbw),title('���ƶ�ֵ��ͼ');

% ����ˮƽͶӰ�������㡢�½��㡢�����ȡ����ȿ��
count1=0;
l=1;
markrow=zeros(1,hight+1);
markrow1=zeros(1,hight+1);
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

%==============��ֱУ����ת=================
grd4 = edge(subcol,'canny'); 
grd4 = bwmorph(grd4,'thicken'); 
% !subplot(2,1,2),imshow(grd4),title('ǿ���߽�ȥˮƽ�߿�ͼ')%�������ǿ���߽�
theta = -90:89;
[R,~] = radon(grd4,theta);%R=�Ƕȣ�xp=λ��
[R1,~] = max(R);
theta_max = 90;
while (theta_max > 20 || theta_max< -20)
    [~,theta_max] = max(R1);                      
    R1(theta_max) = 0; 
%     theta_max = theta_max - 91;
end
angle_x=tan(theta_max);
H=[1,0,0; angle_x,1,0;0,0,1];
T=maketform('affine',H);
subcol=imtransform(subcol,T,'bilinear');
% sbw=imtransform(sbw,T,'bilinear');
% !figure,subplot(2,1,1),imshow(subcol);title('���ƻҶ���ͼ(��ֱУ����ת��)');%���������ת��ĻҶ�ͼ�������ʾ���ƻҶ���ͼ
% !subplot(2,1,2),imshow(sbw);title(['���ƶ�ֵ��ͼ����ֱУ����ת��',int2str(angle_x),'��'],'color','r')%���������ת��ĻҶ�ͼ��
%================��ֱУ����ת����============
waitbar(11/20,h,'��ȴ���ʶ����...');  

%====================ȥ��ֱ�߿��ַ��ָ�======================
sbw=im2bw(subcol);
histcol=sum(sbw);  %���㴹ֱͶӰ
meancol=mean(histcol);%��ֱͶӰ��ƽ��ֵ
mincol=min(histcol);%��ֱͶӰ����Сֵ
levelcol=(meancol+mincol)/4;%!!!!!!��ֱͶӰ��ֵ%��ֱͶӰ��1/4
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
markcol5=markcol3-double(uint16(markcol4./2));%������λ�ã�ȡ����ת����ʽ��
markcol6=diff(markcol5); %�ַ����ľ��루�ַ����ĵ�����һ���ַ����ĵ㣩
[~,findmax]=max(markcol6); %�������ֵ����Ϊ�ڶ��ַ�������ַ����ľ���
markcol6(findmax)=0;
maxwidth=max(markcol6);%�������ֵ����Ϊ����ַ����
%��ȡ�ָ��ַ���ȥ����ֱ�߿򣬲��任Ϊ22��*14�б�׼��ͼ
l=1;
[~,n2]=size(sbw);
cleft_most=[];
cright_most=[];
%!figure;
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
        if l==1
            cleft_most=cleft;
        end
        if l==7
            cright_most=cright;
        end
        SegGray=sbw(:,cleft:cright);
%         subcol1=subcol(:,cleft:cright);
        SegBw2 = imresize(SegGray,[size_x size_y],'bilinear');%�任Ϊ40��*20�б�׼��ͼ ˫���Բ�ֵ        %2013aĬ�Ϸ�ʽ'bicubic'��˫���β�ֵ��
        [SegBw2]=pre_processing(SegBw2);
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

%=================================�ַ�ʶ��============================================
switch shibiefangfa
    case 1 %ģ��ʶ��
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
        elseif l>=3&&l<=5      %��������λ 0~9  A~Z��ĸ������ʶ��
            kmin=1;
            kmax=36;
        else                    %���填��λ 0~9 ����ʶ��
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
        elseif l>=3&&l<=5      %��������λ 0~9  A~Z��ĸ������ʶ��
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


%����ʶ�����
fid=fopen('..\ʶ���ƺ�.txt','w');
fprintf(fid,'%s\n',RegCode);
fclose(fid);
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
end


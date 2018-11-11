function [RegCode]=bp_sim()
global tezheng
%=========����=======================
load bp_net net %��������
N=5;    %����
M=36;     %������������0-9����ĸA-Z
bp_liccode=char(['0':'9' 'A':'Z']);%���� ��ĸ�����б�
size_x =32;
size_y =16;
for l=2:7
%====��ȡ����/������ȡ====
file_name=strcat('..\�и���ͼ\�и���ͼ',int2str(l),'.jpg');%�洢��׼��ͼ
sim_sample=imresize(imread(file_name),[size_x size_y],'bilinear');%��С����
sample=im2bw(sim_sample);%��ֵ��,ģ��ԭͼ
% figure,imshow(sample);title('�и�ԭͼ')
        switch tezheng
        case 1
        %====ͳ��====
        histcol=sum(sample); %���㴹ֱͶӰ
        histrow=sum(sample'); %����ˮƽͶӰ
        p=[histcol,histrow];
        case 2
             %===�ṹ������ȡ====
        sample=bwmorph(sample,'clean');
        sample=bwmorph(sample,'fill');
%        sample=bwmorph(sample,'remove');
%         sample=bwmorph(sample,'thicken');
%         s=strel('disk',2);!!!!!!!strel�����޸�
%         sample=imopen(sample,s);%���������
%         figure, imshow(sample);title('����ͼ')
        sample=bwmorph(sample,'skel',inf);
%         sample=bwmorph(sample,'thin',2);
%         sample=bwmorph(sample,'thin');
%         sample=bwmorph(sample,'thin');
%         sample=bwmorph(sample,'thin');
%         figure, imshow(sample);title('ϸ��ͼ')
        %ֱ��������ȡ
        histcol=sum(sample); %���㴹ֱͶӰ
        histrow=sum(sample'); %����ˮƽͶӰ
        %��ֱͶӰ����
        pos=find(histcol>=10);
        [~,col_size]=size(pos);
        left_col=0;mid_col=0;right_col=0;
        for col=1:col_size
        if pos(col)<5
            left_col=1;
        elseif pos(col)<10
            mid_col=1;
        else
            right_col=1;
        end
        end
        %ˮƽͶӰ����
        pos=find(histrow>=5);
        [~,row_size]=size(pos);
        left_row=0;mid_row=0;right_row=0;
        for row=1:row_size
        if pos(row)<10
            left_row=1;
        elseif pos(row)<20
            mid_row=1;
        else
            right_row=1;
        end
        end
        %��������ȡ ���������
        sample_8=zeros(size_x+2,size_y+2);
        sample_8(2:size_x+1,2:size_y+1)=sample;
        rec_1=zeros(3,2);rec_2=zeros(3,2);rec_3=zeros(3,2);rec_4=zeros(3,2);%Ԥ֧�ڴ�
        count_1=0;count_2=0;count_3=0;count_4=0;
        for i=2:size_x+1
            for j=2:size_y+1
                if sample_8(i,j)==1
                    sample_8(i,j)=sample_8(i-1,j-1)+sample_8(i-1,j)+sample_8(i-1,j+1)+sample_8(i,j-1)+sample_8(i,j+1)+sample_8(i+1,j-1)+sample_8(i+1,j)+sample_8(i+1,j+1);
                end
                switch sample_8(i,j)
                case 1%�׵��ۼ���Ϊ1�����˵�
                   if i<=11&&j<=9
                       rec_1(1,1)=rec_1(1,1)+1;
                   elseif i<=11&&j<=17
                       rec_1(1,2)=rec_1(1,2)+1;
                   elseif i<=21&&j<=9
                       rec_1(2,1)=rec_1(2,1)+1;
                   elseif i<=21&&j<=17
                       rec_1(2,2)=rec_1(2,2)+1;
                   elseif j<=9
                       rec_1(3,1)=rec_1(3,1)+1;
                   else
                       rec_1(3,2)=rec_1(3,2)+1;
                   end
                   count_1=count_1+1;
                case 2%�׵��ۼ���Ϊ2�����յ�
                   if i<=11&&j<=9
                       rec_2(1,1)=rec_2(1,1)+1;
                   elseif i<=11&&j<=17
                       rec_2(1,2)=rec_2(1,2)+1;
                   elseif i<=21&&j<=9
                       rec_2(2,1)=rec_2(2,1)+1;
                   elseif i<=21&&j<=17
                       rec_2(2,2)=rec_2(2,2)+1;
                   elseif j<=9
                       rec_2(3,1)=rec_2(3,1)+1;
                   else
                       rec_2(3,2)=rec_2(3,2)+1;
                   end
                   count_2=count_2+1;
                case 3%�׵��ۼ���Ϊ3���������
                   if i<=11&&j<=9
                       rec_3(1,1)=rec_3(1,1)+1;
                   elseif i<=11&&j<=17
                       rec_3(1,2)=rec_3(1,2)+1;
                   elseif i<=21&&j<=9
                       rec_3(2,1)=rec_3(2,1)+1;
                   elseif i<=21&&j<=17
                       rec_3(2,2)=rec_3(2,2)+1;
                   elseif j<=9
                       rec_3(3,1)=rec_3(3,1)+1;
                   else
                       rec_3(3,2)=rec_3(3,2)+1;
                   end
                   count_3=count_3+1;
                   case 4%�׵��ۼ���Ϊ4�����Ĳ��
                   if i<=11&&j<=9
                       rec_4(1,1)=rec_4(1,1)+1;
                   elseif i<=11&&j<=17
                       rec_4(1,2)=rec_4(1,2)+1;
                   elseif i<=21&&j<=9
                       rec_4(2,1)=rec_4(2,1)+1;
                   elseif i<=21&&j<=17
                       rec_4(2,2)=rec_4(2,2)+1;
                   elseif j<=9
                       rec_4(3,1)=rec_4(3,1)+1;
                   else
                       rec_4(3,2)=rec_4(3,2)+1;
                   end
                   count_4=count_4+1;
                   otherwise
                        continue;
                end
            end
        end
        %��ͨ�պ�����������
%         count_bihe=count_1-count_3*3-count_4*4;
%         if count_bihe<0
%             count_bihe=1;
%         else 
%             count_bihe=0;
%         end
        [L,num] = bwlabel(sample,8);
        count_bihe=num;%��ͨ������
        Feastats = regionprops(L,'basic');%����ͼ�������ϵ�������ߴ磨'Area'�Ǳ������������ͼ����������������ܸ�����'BoundingBox'��1��ndims(L)*2�е���������������Ӧ�������С���Ρ�BoundingBox ��ʽΪ [ul_corner width]������ ul_corner �� [x y z ...] ��������ʽ�����߽���ӵ����Ͻǡ�boxwidth �� [x_width y_width ...] ��ʽָ���߽��������ÿ��ά������ĳ��ȡ���
                                  %'Centroid'��1��ndims(L)�е�����������ÿ����������ģ����ģ���
        Area=[Feastats.Area];%��ͨ���������
        %��ɫ���ط���ͳ��
        white_count=zeros(1,7);
        white_count(1,7)=sum(sum(sample));%�ܸ���
        white_count(1,1)=sum(sum(sample(1:10,1:8)));%����
        white_count(1,2)=sum(sum(sample(1:10,9:16)));%����
        white_count(1,3)=sum(sum(sample(10:20,1:8)));
        white_count(1,4)=sum(sum(sample(10:20,9:16)));
        white_count(1,5)=sum(sum(sample(20:32,1:8)));
        white_count(1,6)=sum(sum(sample(20:32,9:16)));
        %�Խ��߰�ɫ���ظ���ͳ��
        white_co1=0;white_co2=0;
        for j=1:16
                white_co1=white_co1+sample(2*j,j);
        end
        for j=1:16
                white_co2=white_co2+sample(32-2*j+1,j);
        end
        end
        %���߰�ɫ���ظ���ͳ��
        white_co3=sum(sample(:,8));
        white_co4=sum(sample(16,:));
        %=========������ȡ����=========
        p=[left_col,mid_col,right_col,left_row,mid_row,right_row,reshape(rec_1,[1,6]),reshape(rec_2,[1,6]),reshape(rec_3,[1,6]),reshape(rec_4,[1,6]),count_bihe,Area,...
                    white_count,white_co1,white_co2,white_co3,white_co4];

p=p';
%=========������ȡ����==========
A=sim(net,p);
[~,position] = max(A);
RegCode(l)=bp_liccode(position);
%!subplot(1,6,l-1),imshow(sample)
end
RegCode
%title(['ʶ���ƺ���:', RegCode],'Color','r')
end
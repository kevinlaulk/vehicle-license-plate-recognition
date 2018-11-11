function [RegCode]=bp_sim()
global tezheng
%=========仿真=======================
load bp_net net %下载网络
N=5;    %组数
M=36;     %样本数：数字0-9，字母A-Z
bp_liccode=char(['0':'9' 'A':'Z']);%数字 字母索引列表
size_x =32;
size_y =16;
for l=2:7
%====读取样本/特征提取====
file_name=strcat('..\切割子图\切割子图',int2str(l),'.jpg');%存储标准子图
sim_sample=imresize(imread(file_name),[size_x size_y],'bilinear');%大小重塑
sample=im2bw(sim_sample);%二值化,模板原图
% figure,imshow(sample);title('切割原图')
        switch tezheng
        case 1
        %====统计====
        histcol=sum(sample); %计算垂直投影
        histrow=sum(sample'); %计算水平投影
        p=[histcol,histrow];
        case 2
             %===结构特征提取====
        sample=bwmorph(sample,'clean');
        sample=bwmorph(sample,'fill');
%        sample=bwmorph(sample,'remove');
%         sample=bwmorph(sample,'thicken');
%         s=strel('disk',2);!!!!!!!strel参数修改
%         sample=imopen(sample,s);%开运算操作
%         figure, imshow(sample);title('背景图')
        sample=bwmorph(sample,'skel',inf);
%         sample=bwmorph(sample,'thin',2);
%         sample=bwmorph(sample,'thin');
%         sample=bwmorph(sample,'thin');
%         sample=bwmorph(sample,'thin');
%         figure, imshow(sample);title('细化图')
        %直线特征提取
        histcol=sum(sample); %计算垂直投影
        histrow=sum(sample'); %计算水平投影
        %垂直投影分类
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
        %水平投影分类
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
        %特征点提取 八邻域计算
        sample_8=zeros(size_x+2,size_y+2);
        sample_8(2:size_x+1,2:size_y+1)=sample;
        rec_1=zeros(3,2);rec_2=zeros(3,2);rec_3=zeros(3,2);rec_4=zeros(3,2);%预支内存
        count_1=0;count_2=0;count_3=0;count_4=0;
        for i=2:size_x+1
            for j=2:size_y+1
                if sample_8(i,j)==1
                    sample_8(i,j)=sample_8(i-1,j-1)+sample_8(i-1,j)+sample_8(i-1,j+1)+sample_8(i,j-1)+sample_8(i,j+1)+sample_8(i+1,j-1)+sample_8(i+1,j)+sample_8(i+1,j+1);
                end
                switch sample_8(i,j)
                case 1%白点累计量为1，即端点
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
                case 2%白点累计量为2，即拐点
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
                case 3%白点累计量为3，即三叉点
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
                   case 4%白点累计量为4，即四叉点
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
        %连通闭合曲线数计算
%         count_bihe=count_1-count_3*3-count_4*4;
%         if count_bihe<0
%             count_bihe=1;
%         else 
%             count_bihe=0;
%         end
        [L,num] = bwlabel(sample,8);
        count_bihe=num;%连通区个数
        Feastats = regionprops(L,'basic');%计算图像区域的系列特征尺寸（'Area'是标量，计算出在图像各个区域中像素总个数。'BoundingBox'是1行ndims(L)*2列的向量，即包含相应区域的最小矩形。BoundingBox 形式为 [ul_corner width]，这里 ul_corner 以 [x y z ...] 的坐标形式给出边界盒子的左上角、boxwidth 以 [x_width y_width ...] 形式指出边界盒子沿着每个维数方向的长度。）
                                  %'Centroid'是1行ndims(L)列的向量，给出每个区域的质心（重心）。
        Area=[Feastats.Area];%连通区区域面积
        %白色像素分区统计
        white_count=zeros(1,7);
        white_count(1,7)=sum(sum(sample));%总个数
        white_count(1,1)=sum(sum(sample(1:10,1:8)));%左上
        white_count(1,2)=sum(sum(sample(1:10,9:16)));%右上
        white_count(1,3)=sum(sum(sample(10:20,1:8)));
        white_count(1,4)=sum(sum(sample(10:20,9:16)));
        white_count(1,5)=sum(sum(sample(20:32,1:8)));
        white_count(1,6)=sum(sum(sample(20:32,9:16)));
        %对角线白色像素个数统计
        white_co1=0;white_co2=0;
        for j=1:16
                white_co1=white_co1+sample(2*j,j);
        end
        for j=1:16
                white_co2=white_co2+sample(32-2*j+1,j);
        end
        end
        %中线白色像素个数统计
        white_co3=sum(sample(:,8));
        white_co4=sum(sample(16,:));
        %=========特征提取结束=========
        p=[left_col,mid_col,right_col,left_row,mid_row,right_row,reshape(rec_1,[1,6]),reshape(rec_2,[1,6]),reshape(rec_3,[1,6]),reshape(rec_4,[1,6]),count_bihe,Area,...
                    white_count,white_co1,white_co2,white_co3,white_co4];

p=p';
%=========特征提取结束==========
A=sim(net,p);
[~,position] = max(A);
RegCode(l)=bp_liccode(position);
%!subplot(1,6,l-1),imshow(sample)
end
RegCode
%title(['识别车牌号码:', RegCode],'Color','r')
end
function [SegBw2]=pre_processing(SegBw2)
        size_x=40;
        size_y=20;
        SegBw1=im2bw(SegBw2);
        [i,j]=find(SegBw1==1);%图像预处理
        imin=min(i);
        imax=max(i);
        jmin=min(j);
        jmax=max(j);
        SegBw1=SegBw1(imin:imax,jmin:jmax);
        SegBw2=imresize(SegBw1,[size_x,size_y],'bilinear');%变换为40行*20列标准子图
        SegBw2=im2bw(SegBw2);
end

function [SegBw2]=pre_processing(SegBw2)
        size_x=40;
        size_y=20;
        SegBw1=im2bw(SegBw2);
        [i,j]=find(SegBw1==1);%ͼ��Ԥ����
        imin=min(i);
        imax=max(i);
        jmin=min(j);
        jmax=max(j);
        SegBw1=SegBw1(imin:imax,jmin:jmax);
        SegBw2=imresize(SegBw1,[size_x,size_y],'bilinear');%�任Ϊ40��*20�б�׼��ͼ
        SegBw2=im2bw(SegBw2);
end

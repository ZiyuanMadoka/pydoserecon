% MATLAB script
% Rib_detection.m
% detect landmarks and retrive measurements from abdominal DRR in mhd format.
% Copyright (c) 2019 Ziyuan Wang

function LMs = Rib_detection(Initdir, filename, ii, scale,ref_size)
    global vloc20;
    global I20;
    paren = @(x,varargin) x(varargin{:}); % indexing function-return results
    [img, ] = read_mhd(strcat(Initdir,filename));  
    DRR = img.data';
    temp_spacing = img.spacing; % pixel size in two dimensions
    margin_mm = @(length) round(length/temp_spacing(1));
    
    % preprocessing to remove the window and apply averaging filter 3*3 %imshow(I2);
    [img_start,I2] = precrop_filt(DRR,3);  
    [I_height, I_width] = size(I2);  % height and width of the image
    
    % initial estimation of vertebral column
    h_Sig = sum(I2(1:fix(I_height/2),:),1);  % horizontal_signal: sum along column, half of the columns
    [rib_last1,rib_width1] = find_comp(h_Sig,1,0.3);  % select the rib region based on threshold, 0.3 between max and min
    rib_start1 = rib_last1-rib_width1+1;
    h_dev = abs(diff(h_Sig));
    central_start = int16(max(rib_start1 + rib_width1*0.25, I_width*0.3));
    central_end = int16(min(rib_last1 - rib_width1*0.25,I_width*0.7));
    signal =smooth(smooth(h_dev(central_start:central_end))); % take the middle half of the h_Sig
    [peaks,loc]=findpeaks(signal,'MinPeakWidth',margin_mm(2),'MinPeakProminence',max(signal)*0.3,'NPeaks',2,'MinPeakDistance',margin_mm(12)); 
    ver_col_start=central_start+loc(1)-margin_mm(2);  
    ver_col_end = central_start+loc(2)+margin_mm(2); 
    
    %initial estimation of intervertebral disks
    v_verte = sum(I2(:,ver_col_start-margin_mm(4):ver_col_end+margin_mm(4)),2); 
    [vpeak1, vloc1]=findpeaks(-smooth(v_verte),'MinPeakDistance',I_height/12, 'NPeaks',9);  
    % minpeakdistance should be coorparated with pixel size of the image, it was 30*0.391 = 11.73 mm 
    vloc1 = vloc1(vloc1>10);
    if length(vloc1)<8  
        [vpeak1, vloc1]=findpeaks(diff(smooth(v_verte)),'MinPeakDistance',I_height/12, 'NPeaks',9);  % was 9, if the image was cut
        % minpeakdistance should be coorparated with pixel size of the image, it was 30*0.391 = 11.73 mm 
    end
    vloc1 = vloc1(vloc1>10); % exclude the peak location <10
   
    % estimatino of vertebral column of each bone
    length_x = min(length(vloc1),8);  % maximum 8 x points
    x1 =zeros(1,length_x);x2 = x1;vlocm = x1;
    loc_seg0 = [ver_col_start, ver_col_end];
    for i = 0:length_x-1
        xrange = ver_col_start-margin_mm(20):ver_col_end+margin_mm(20); % -25 and +25 used to be 50 * 0.391 = 19.55 mm
        minpeakdistance=margin_mm(16);
        if i==0 
            v_sig_seg = sum(I2(1:vloc1(i+1),:),1); % vector_signal: sum along column, half column
            vlocm(i+1) = (1+vloc1(i+1))/2;
        else
            if i > length(vloc1)-3
                xrange = ver_col_start-margin_mm(25):ver_col_end+margin_mm(25);
            end
            v_sig_seg = sum(I2(vloc1(i):vloc1(i+1),:),1); % vector_signal: sum along column, half column
            vlocm(i+1) = (vloc1(i)+vloc1(i+1))/2;
            minpeakdistance=max(margin_mm(16),loc_seg0(2)-loc_seg0(1)-margin_mm(10));
        end        
        v_dev_seg = abs(diff(v_sig_seg)); %.* (v_sig_seg(1:length(v_sig_seg)-1)/max(v_sig_seg)).^3;
        minpeakprominence = max(smooth(smooth((v_dev_seg(xrange)))))*0.15; %25 was 50 * 0.391 = 19.55 mm, was *0.2
        [peaks,loc_seg]=findpeaks(smooth(smooth(v_dev_seg(xrange))),'MinPeakWidth',margin_mm(3),'MinPeakProminence',minpeakprominence,'NPeaks',3,'MinPeakDistance',minpeakdistance);
        idxleft = loc_seg< length(xrange)/2;   
        if(length(loc_seg(idxleft)) >= 1 && length(loc_seg(~idxleft)) >= 1)  
            [~,Ileft]=sort(peaks(idxleft),'descend');
            [~,Iright] = sort(peaks(~idxleft),'descend');
            loc_seg0 = [paren(loc_seg(idxleft),Ileft(1)),paren(loc_seg(~idxleft),Iright(1))];
            x1(i+1) = loc_seg0(1)+xrange(1);
            x2(i+1) = loc_seg0(2)+xrange(1);
        else
            x1(i+1) = nan;
            x2(i+1) = nan;
            loc_seg = loc_seg0;
        end
    end
    if (length(x1) == 8) && (abs(x1(end)-mean(x1(6:7)))<margin_mm(10)) && (abs(x2(end)-mean(x2(6:7)))<margin_mm(10)) 
        x1 = round(smooth(x1))';x2 = round(smooth(x2))';
    else
        x1 = round(smooth(x1(1:7)))'; x2 = round(smooth(x2(1:7)))';
    end
    
    % estimate the tilt angle of the vertebral column
    Fit_c3 = polyfit(vlocm(1:length(x1)),(x1+x2)/2,1);
    colli = atan(-Fit_c3(1))/pi*180;
    if round(colli) == 0
        collimator_angle = colli;
    else
        collimator_angle = round(colli/abs(colli)*round(abs(colli)));
    end
    
    %second round estimation for vlocs
    v_verte_2 = v_verte; %copy the size of h_verte
    if length(x1) < 8 
        for kk =(length(x1)+1):8
            x1(kk) = x1(end);
            x2(kk) = x2(end);
        end
    end
    x1_ex =[x1(1),x1];
    x2_ex = [x2(2),x2];
    for i = 1: length(vlocm)+1
        if i==1
            block = 1:vlocm(1);
        elseif i == length(vlocm)+1
            block = round(vlocm(end)+1):I_height;
        else
            block = round(vlocm(i-1)+1):vlocm(i);
        end
        v_verte_2(block) =sum(I2(block,x1_ex(i)-margin_mm(4):x2_ex(i)+margin_mm(4)),2); % -5 and +5  
    end
    MPP=(max(-smooth(v_verte_2))-min(-smooth(v_verte_2)))*0.001;
    [~, vloc2]=findpeaks(-smooth(v_verte_2),'MinPeakDistance',I_height/12, 'MinPeakProminence',MPP, 'NPeaks',9);  
    vloc2 = vloc2(vloc2>margin_mm(9));
    if ii>1 && vloc2(end)>(I_height-margin_mm(6))
        Isize0 = size(I20);
        Isize = size(I2);
        vloc2(end) = vloc20(end)*Isize(1)/Isize0(1);
    end
    vloc2 = vloc2(vloc2<=(I_height-margin_mm(6))); 

    if ii>1 && length(vloc2)<8
        Isize0 = size(I20);
        Isize = size(I2);
        vloc2 = vloc20*Isize(1)/Isize0(1);
    end
    vlocm2= vlocm;
    y1=zeros(size(x1));
    y2 = y1;
    for i=1:length(vloc2)-1
        if i==1
            vlocm2(i) = (1+vloc2(i))/2;
        else
            vlocm2(i) = (vloc2(i-1)+vloc2(i))/2;
        end
        Dy = (x2(i)-x1(i))/2*tan(Fit_c3(1)); % delta y to correct for bending
        y1(i)=vlocm2(i)+Dy;
        y2(i)=vlocm2(i)-Dy;
    end
    
    % rotate the image based on tilt angle
    Irot = imrotate(I2,collimator_angle);
    I3 = Irot(crop_rot(Irot),:);
    [I3_height,I3_width]=size(I3);
    [~,L1_rot] = rotate((x1(3)+x2(3))/2-I_width/2, vloc2(3)-I_height/2,-collimator_angle);
    v3_Sig = sum(I3(1:fix(I3_height/2),:),1);  % vector_signal: sum along column, half column
    threshold = 0.3;
    [rib_left_rot,rib_width] = find_comp2(v3_Sig,1,threshold);
    rib_last = rotate(rib_left_rot-I3_width/2,L1_rot,collimator_angle)+I_width/2;
    rib_start= rotate(rib_left_rot-rib_width+1-I3_width/2,L1_rot,collimator_angle)+I_width/2;
    Fit_v =[-Fit_c3(1),(1/Fit_c3(1)+Fit_c3(1))*polyval(Fit_c3,vloc2(3))-Fit_c3(2)/Fit_c3(1)];
    Rib_right = [rib_start,Fit_v(1)*rib_start+Fit_v(2)];
    Rib_left = [rib_last,Fit_v(1)*rib_last+Fit_v(2)];
    L1 = [polyval(Fit_c3,vloc2(3)),vloc2(3)];
    
    rib_width_right = sqrt((Rib_right(1)-L1(1))^2+(Rib_right(2)-L1(2))^2);
    rib_width_left = sqrt((Rib_left(1)-L1(1))^2+(Rib_left(2)-L1(2))^2);
    vloc20 = vloc2;
    I20 = I2;
    
    % calculate landmarks and covert back to physical size size
    scale_factor = ref_size./(img.size.*img.spacing); % scaled for reference
  
    % change middle of T11 to middle of L4
    T11L4=sqrt(((polyval(Fit_c3,vlocm2(2))-polyval(Fit_c3,vlocm2(7)))*scale_factor(1))^2+((vlocm2(7)-vlocm2(2))*scale_factor(2))^2)*img.spacing(1);
    % assign values to a structure to store landmarks
    LMs = struct;
    LMs.Rib_width = scale_factor(1)*rib_width*img.spacing(1);
    LMs.T11L4 = T11L4;
    LMs.collimator = mod(-collimator_angle,360);
    LMs.Th10_bottom_xyz = coordi_CT(polyval(Fit_c3,vloc2(1)),vloc2(1),img,img_start);
    LMs.L4_bottom_xyz = coordi_CT(polyval(Fit_c3,vloc2(7)),vloc2(7),img,img_start);
    LMs.L1_bottom_xyz = coordi_CT(polyval(Fit_c3,vloc2(4)),vloc2(4),img,img_start);
    LMs.Th12_right_cor_xyz = coordi_CT(x1(3),y1(3),img,img_start);
    LMs.Th12_left_cor_xyz = coordi_CT(x2(3),y2(3),img,img_start);
    LMs.L2_right_cor_xyz = coordi_CT(x1(5),y1(5),img,img_start);
    LMs.L2_left_cor_xyz = coordi_CT(x2(5),y2(5),img,img_start);
    LMs.Rib_right_xyz = coordi_CT(Rib_right(1),Rib_right(2),img,img_start);
    LMs.Rib_left_xyz = coordi_CT(Rib_left(1),Rib_left(2),img,img_start);
    LMs.Verte_width = scale_factor(1)*mean(x2-x1)/cosd(collimator_angle);
    LMs.Rib_width_left = scale_factor(1)*rib_width_left * img.spacing(1);
    LMs.Rib_width_right = scale_factor(1)*rib_width_right * img.spacing(1);
end



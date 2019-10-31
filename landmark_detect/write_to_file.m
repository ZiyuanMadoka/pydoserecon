% MATLAB function
% write_to_file.m
% Copyright (c) 2019 Ziyuan Wang
% write the measurements saved in struct LMs to a txt file saved as outputfile

function write_to_file(outputfile, LMs)
    fid = fopen(outputfile,'w');
    str_f = '%% measurements from DRR, size(mm) ,angle (degree), system coordinate system as DRR \r\n';
    formatSpec = ["Rib_width %4.2f \r\n","Length_T11L4 %4.2f \r\n", "Collimator_angle %d \r\n","Th10_bottom_xyz %4.2f %4.2f %4.2f \r\n",...
        "L4_bottom_xyz %4.2f %4.2f %4.2f \r\n","L1_bottom_xyz %4.2f %4.2f %4.2f \r\n",...
        "Th12_right_cor_xyz %4.2f %4.2f %4.2f \r\n", "Th12_left_cor_xyz %4.2f %4.2f %4.2f \r\n",...
        "L2_right_cor_xyz %4.2f %4.2f %4.2f \r\n","L2_left_cor_xyz %4.2f %4.2f %4.2f \r\n",...
        "Rib_left_xyz %4.2f %4.2f %4.2f \r\n", "Rib_right_xyz %4.2f %4.2f %4.2f \r\n",...
        "Verte_width %4.2f \r\n","Rib_width_left %4.2f \r\n","Rib_width_right %4.2f \r\n"];
    fprintf(fid,str_f);
    %fprintf(fid,'%% measurements from DRR, scaled size 100 * 100 pixel,angle (degree), system coordinate system as DRR \r\n'); 
    %fprintf(fid,'%% the widht of the widest part of the rib \r\n');
    fprintf(fid,formatSpec(1),LMs.Rib_width);
    fprintf(fid,formatSpec(2),LMs.T11L4);
    %fprintf(fid,'%% considering bending, physical length, including Th11 and L5 \r\n');
    fprintf(fid,formatSpec(3),LMs.collimator);
    fprintf(fid,formatSpec(4),LMs.Th10_bottom_xyz);
    %fprintf(fid,'%% middle_bottem of TH10 \r\n');
    fprintf(fid,formatSpec(5),LMs.L4_bottom_xyz);
    fprintf(fid,formatSpec(6),LMs.L1_bottom_xyz);
    %fprintf(fid,'%% middle_bottem of L5 \r\n');
    fprintf(fid,formatSpec(7),LMs.Th12_right_cor_xyz);
    %fprintf(fid,'%% Right_middle of boundary of L5 \r\n');
    fprintf(fid,formatSpec(8),LMs.Th12_left_cor_xyz);
    %fprintf(fid,'%% Left_middle of boundary of Th12 \r\n');
    fprintf(fid,formatSpec(9),LMs.L2_right_cor_xyz);
    %fprintf(fid,'%% Right_middle of boundary of L5 \r\n');
    fprintf(fid,formatSpec(10),LMs.L2_left_cor_xyz);
    %fprintf(fid,'%% Rib left boundary at L1 \r\n');
    fprintf(fid,formatSpec(11),LMs.Rib_left_xyz);
    %fprintf(fid,'%% Right right boundary at L1 \r\n');
    fprintf(fid,formatSpec(12),LMs.Rib_right_xyz);
    %fprintf(fid,'%% Vertebral width \r\n');
    fprintf(fid,formatSpec(13),LMs.Verte_width);
    %fprintf(fid,'%% rib width left part \r\n');
    fprintf(fid,formatSpec(14),LMs.Rib_width_left);
    %fprintf(fid,'%% rib width right part \r\n');
    fprintf(fid,formatSpec(15),LMs.Rib_width_right);
    fclose(fid);
end
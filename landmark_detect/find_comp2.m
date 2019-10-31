function[idx_end, max_com] = find_comp2(DRR,dim,threshold)
    sum_D = sum(DRR,dim);
    if (dim==2)
        low_bg = min(min(DRR));
    else
        low_bg = min(sum_D);
    end
    sum_D = smooth(sum_D);
    idx_bg = sum_D> low_bg+(max(sum_D)-min(sum_D))*threshold; % rows, threshold:[0,1]
    idx_bg2= sum_D> low_bg+(max(sum_D)-min(sum_D))*0.1;
    component =1; pre = 1; max_com = 0;idx_end = length(idx_bg);
    idx_pass = find(idx_bg);
    %max_com = idx_end-idx_pass(1)+1;
    for kk = idx_pass(1):idx_pass(end)
        if(idx_bg2(kk)==0 || kk == idx_pass(end))
            pre = 0;
            if(component >max_com)
                max_com = component;
                idx_end = kk;
            end
            component =0;
        elseif(idx_bg2(kk)==1 && pre==1 )
            component=component+1;
        else
            component = 1;
            pre = 1;
        end
    end
end
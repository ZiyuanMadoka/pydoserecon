function[idx_end, max_com] = find_comp(DRR,dim,threshold)
    sum_D = sum(DRR,dim);
    if (dim==2)
        low_bg = min(min(DRR));
    else
        low_bg = min(sum_D);
    end
    sum_D = smooth(sum_D);
    idx_bg = sum_D> low_bg+(max(sum_D)-min(sum_D))*threshold; % rows, threshold:[0,1]
    component =1; pre = 1; max_com = 0;idx_end = length(idx_bg);
    for kk = 1:length(idx_bg)
        if(idx_bg(kk)==0 || kk == length(idx_bg))
            pre = 0;
            if(component >max_com)
                max_com = component;
                idx_end = kk;
            end
            component =0;
        elseif(idx_bg(kk)==1 && pre==1 )
            component=component+1;
        else
            component = 1;
            pre = 1;
        end
    end
end
% Its function that calculates the half power beam width (HPBW)
% a is a matrix that the first column is the degress
% the second column is the matching power in elevation plane
% the third column is the matching power in aznimuth plane
% the powers is supposed no to be normalized
function [hpbw_el,hpbw_az]=hpbw_sc(a)

deg=a(:,1);
pow_el=a(:,2);
pow_az=a(:,3);

pow_el=pow_el-max(pow_el);
pow_az=pow_az-max(pow_az);

ind_el=find(pow_el>-3.2);
ind_az=find(pow_az>-3.2);

angles_el=deg(ind_el);
angles_az=deg(ind_az);



max_deg_up_ind_el=max(find(angles_el<90));
max_deg_down_ind_el=min(find(angles_el>270 & angles_el<360));

max_deg_up_ind_az=max(find(angles_az<90));
max_deg_down_ind_az=min(find(angles_az>270 & angles_az<360));




max_deg_up_el=angles_el(max_deg_up_ind_el);
max_deg_down_el=360-angles_el(max_deg_down_ind_el);

hpbw_el=max_deg_up_el+max_deg_down_el;


max_deg_up_az=angles_az(max_deg_up_ind_az);
max_deg_down_az=360-angles_az(max_deg_down_ind_az);

hpbw_az=max_deg_up_az+max_deg_down_az;

end

function [my,bins] = window_mean(x,y,bins)
my   = zeros(1,length(bins)-1);
for it = 2:length(bins)
    loc     =  y > bins(it-1) &  y < bins(it);
    my(it)  =  nanmean(x(loc));
end
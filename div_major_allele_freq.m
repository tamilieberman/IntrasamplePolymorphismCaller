function [majorAF, majorNT,minorNT, minorAF] = div_major_allele_freq(cnts)

%Tami Lieberman, 2012, Kishony lab

c=cnts(1:4,:,:)+cnts(5:8,:,:);

[sorted, sortedpositions] = sort(c,1);
maxcount = sorted(end,:,:);
minorcount = sorted(end-1,:,:);

majorAF=double(maxcount)./sum(c,1);
minorAF=double(minorcount)./sum(c,1);

majorNT = squeeze(sortedpositions(end,:,:));
minorNT = squeeze(sortedpositions(end-1,:,:));

majorAF=squeeze(majorAF);
minorAF=squeeze(minorAF);

return


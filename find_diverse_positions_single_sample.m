function p = find_diverse_positions_single_sample(infile, parameterfile, reportfile, MAF_control)

%Tami Lieberman, 2015, Kishony lab

%Names of parameters should be self explanatory when combined with
%description of data matrix in pileup_to_diversity_matrix.m

%All tests are strictly greater or less than the provided threshold


load(infile);
load(parameterfile);

if nargin <4
    MAF_control=ones(1,numel(size(data,2)));
end
    

%get coverage histogram
%coveragethresholds contains cdf cutoffs .01:.01:1.0.
%only counts positions where at least 1 read aligned
%purpose of this data structure is to allow downstream processes to
%remove positions with excess coverage in a way that accounts for
%variation in coverage between samples
coveragethresholds=zeros(100,1);
cov=sum(double(data(1:8,:)));
cov(cov<1)=[];
if isempty(cov) | sum(cov)/length(data)< 0.5
    error('Not enough coverage')
end
cov=sort(cov);
cutoffs=.01:.01:1;
for j=1:numel(cutoffs);
    coveragethresholds(j)=cov(floor(cutoffs(j)*numel(cov)));
end
    


num_reqs=11;

NPositions=size(data,2);

%unpack parameters
minorfreqthreshold=params.minorfreqthreshold;
minreads_perstrand=params.minreads_perstrand;
maxreads_perstrand=coveragethresholds(params.maxreads_perstrand_percentile);
min_bq=params.min_bq;
min_mq=params.min_mq;
min_td=params.min_td;
max_td=params.max_td;
max_sbp=params.max_sbp;
max_bqp=params.max_bqp;
%max_mqp=p.max_mqp;
max_tdp=params.max_tdp;
max_percent_indels=(params.max_percent_indels)/100;
max_percent_ends=(params.max_percent_ends)/100;
min_control_MAF=params.min_control_MAF;
minreads_perstrand_per_allele=params.minreads_perstrand_per_allele;


[~, majorNT, minorNT] = div_major_allele_freq(data);

positionsv=(1:NPositions)';

n1=majorNT';
n2=minorNT';

minorfreq=double((data(sub2ind(size(data),n2,positionsv))+data(sub2ind(size(data),n2+4,positionsv))))'./sum(data(1:8,:));
readsf=sum(data(1:4,:));
readsr=sum(data(5:8,:));
f2 = data(sub2ind(size(data),n2,positionsv)); %minor allele counts on forward strand
r2 = data(sub2ind(size(data),n2+4,positionsv)); %minor allele counts on reverse strand
majorbqf=data(sub2ind(size(data),n1+8,positionsv));
majorbqr=data(sub2ind(size(data),n1+12,positionsv));
minorbqf=data(sub2ind(size(data),n2+8,positionsv));
minorbqr=data(sub2ind(size(data),n2+12,positionsv));
majormqf=data(sub2ind(size(data),n1+16,positionsv));
majormqr=data(sub2ind(size(data),n1+20,positionsv));
minormqf=data(sub2ind(size(data),n2+16,positionsv));
minormqr=data(sub2ind(size(data),n2+20,positionsv));
majortdf=(data(sub2ind(size(data),n1+24,positionsv)));
minortdf=(data(sub2ind(size(data),n2+24,positionsv)));
majortdr=(data(sub2ind(size(data),n1+28,positionsv)));
minortdr=(data(sub2ind(size(data),n2+28,positionsv)));
percent_indels=double(data(end,:))./(sum(data(1:8,:))+double(data(end,:)));
percent_ends=double(data(end-1,:))./sum(data(1:8,:));

SBp=data(end-6,:);
BQp=data(end-5,:);
%MQp=data(end-4,:);
TDFp=data(end-3,:);
TDRp=data(end-2,:);



%Find true/false of meeting thresholds
Tminor = minorfreq > minorfreqthreshold;
Treads= (readsf > minreads_perstrand) & (readsr > minreads_perstrand) &...
    ((readsf +readsr) < maxreads_perstrand) & (f2' > minreads_perstrand_per_allele) ...
    & (r2' > minreads_perstrand_per_allele);
Tbq= ((majorbqf > min_bq) & (minorbqf > min_bq) & (majorbqr > min_bq) & (minorbqr > min_bq))';
Tmq = ((majormqf > min_mq) & (minormqf > min_mq) & (majormqr > min_mq) & (minormqr > min_mq))';
Ttd = ((majortdf > min_td) & (majortdf < max_td) & (majortdr < max_td) & (majortdr > min_td)...
    & (minortdf > min_td) & (minortdf < max_td) & (minortdr > min_td) & (minortdr < max_td))';
Tid = percent_indels < max_percent_indels;
Te = percent_ends < max_percent_ends;


TSBp = SBp < max_sbp;
TBQp = BQp < max_bqp;
%TMQp = MQp < max_mqp;
TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);

if numel(MAF_control > 1)
    control = MAF_control > min_control_MAF;
else
    control = ones(Npositions,1);
end



%Report how many positions met each requirement
f=fopen(reportfile ,'w');
fprintf(f,'MinorAlleleFreq: %g  \n',sum(Tminor)) ;
fprintf(f,'Cov: %g  \n',sum(Treads)) ;
fprintf(f,'minBQ: %g  \n',sum(Tbq)) ;
fprintf(f,'minMQ: %g  \n',sum(Tmq)) ;
fprintf(f,'SBp: %g  \n',sum(TSBp)) ;
fprintf(f,'BQp: %g  \n',sum(TBQp)) ;
%fprintf(f,'MQp: %g  \n',sum(TMQp)) ;
fprintf(f,'TDp: %g  \n',sum(TTDp)) ;
fprintf(f,'maxIndels: %g  \n',sum(Tid)) ;
fprintf(f,'maxEndsOfReads: %g  \n',sum(Te)) ;
fprintf(f,'acceptableTD: %g  \n',sum(Ttd)) ;


%Records positions that met all requirements
allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + Te + TSBp + TBQp +TTDp+ control; % TMQp 
fprintf(f,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;
p(allreqs==num_reqs)=1;
p=find(p);

%Report how many positions are removed because of Isogenic control
good=Tminor + Treads + Tbq + Tmq + Ttd + Tid + Te + TSBp + TBQp +TTDp;
removedbycontrol= good==(num_reqs-1) &  ~control;
fprintf(f,'Positions only removed because of isogenic control: %g  \n',sum(removedbycontrol));
fprintf(f,'\n');

   

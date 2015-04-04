function maf=maf_from_matrix(filename)

%Tami Lieberman, 2015, Kishony lab

load(filename)
[maf, ~, ~] = div_major_allele_freq(data);

end

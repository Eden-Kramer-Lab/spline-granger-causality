function [ out ] = pval2grstat( in, Type)
% PVAL2GRSTAT = computes the associated p-value to the
% Grendander-Rosenblatt test statistic and vice versa.
%
% INPUTS:
% . in   = p-value or gr statisitc
% . type = indicates whether the value of 'in' 'pvalue' or 'grstatistic'
%
% OUTPUTS:
% . out  = computes the p-value if input is 'grstatistic' or the
%          grstatistic if input is 'pvalue'
%
% Example: pval2grstat(val, 'pvalue') converts val from a pvalue to its
% associate GR statistic



pval  = 1 -[0.999 0.99 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 ...
    0.45 0.4 0.35 0.3 0.18524 0.10267 0.04136 0.00916 ]';

grstat = [3.4808 2.8070 2.2414 1.96 1.7805 1.6449 1.5341 1.4395 1.3562 ...
    1.2812 1.2126 1.149 1.0892 1.0322 0.9774 0.9238 0.8 0.7 0.6 0.5]';


pvalii=[min(pval):0.001:max(pval)]';
grstatii = interp1q(pval,grstat,pvalii) ;

nelectrodes = length(in);

if strcmp(Type, 'pvalue') %% if input type is pvalue, output type is grstatistic
    fprintf('Returning gr statistic:')
    [~,i] = min(abs(pvalii-in));
    out = grstatii(i);
    
elseif strcmp(Type, 'grstatistic') %% if input type is grstatistic, output type pvalue
    fprintf('Returning p value:')
    if in > 3.6
        fprintf('GR statistic > 3.6; lots of evidence to reject')
        out = 0.0;
    elseif in < 0.4
        fprintf('GR statistic < 0.4; not enough evidence to reject')
        out = 1.0;
    else
        [~,i] = min(abs(grstatii-in));
        out = pvalii(i);
    end
else
    fprintf('No type specified')
    out = NaN;
end

end


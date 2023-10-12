function p = bootstrapping(a,times)

% a is [1x(participants)]
% bootstrapping across pariticipants

b = [];
for t = 1:times
    bootstrap_sample = randsample(a,length(a),true); % random sapmle with replacements
    b(t) = mean(bootstrap_sample);
end  
b_shift = b - mean(b);
%p = (sum(b_shift>abs(mean(a)))+sum(b_shift<-abs(mean(a))))/times;
%p = 2*min([sum(b_shift>abs(mean(a))),sum(b_shift<-abs(mean(a)))])/times;
p = 2*min([sum(b_shift<mean(a)),sum(b_shift>mean(a))])/times;

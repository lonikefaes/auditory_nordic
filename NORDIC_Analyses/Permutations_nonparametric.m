function [p] = Permutations_nonparametric(dataVar, NumSub, contrasts)

for i = 1:size(contrasts,1)

    contr = contrasts(i,:);

    tmp  = diff(dataVar([contr],:));

    ObsStat = tmp*ones(10,1);

    Permutations = zeros(2^NumSub-2,NumSub);
    count1       = 1;
    for ind = 1:NumSub-1
        % cycling from 1 to n-1 rather than 0 to n
        C = nchoosek(1:NumSub,ind);
        P = ones(size(C,1),NumSub);
        for indP = 1:size(C,1)
            P(indP,C(indP,:))=-1;
        end
        Permutations(count1 + (0:size(P,1)-1),:) = P;
        count1  = count1 + size(P,1);
    end

    Permstats = tmp*Permutations';                                          % Histogram centered around 0 for all permutations.
    p(i) = (sum(abs(Permstats)>= abs(ObsStat) ) + 1)/(numel(Permstats)+2);  % p-value

end

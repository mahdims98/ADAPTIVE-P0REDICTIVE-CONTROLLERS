function P = remove_first_zeros_threshold(P, threshold)  
    index = -1;
    for i=1:length(P)
        if abs(P(i)) > threshold
            index = i;
            break
        end
    end
    if index ~=-1
        P = P(index:end);
    end
end
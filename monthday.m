function [M,D] = monthday(doy)

    First_of_month = [1,32,60,91,121,152,182,213,244,274,305,335,366];
    
    M = find(First_of_month > doy,1)-1;
    D = doy-First_of_month(M)+1;
end
function [subCorMask, subEle] = find_subCorMask(initialSci,subSci,initialRegion)
    %   find correction region in subSci
    %   inputs:
    %   -- initialSci: large science plane that includes subSci
    %   -- subSci: sub-science plane
    %   -- initialRegion: regions that includes correction mask in original
    %   science plane
    %   outputs:
    %   -- subCorMask: correction mask in subSci plane
    %   -- subEle: indices of correction mask in subSci plane
    
    numEles = length(initialRegion.CorEle);
    subCorMask = zeros(size(subSci.xxlD));
    for iq = 1:numEles
        ipx = find(subSci.xxlD == initialSci.xxlD(initialRegion.CorEle(iq)));
        ipy = find(subSci.yylD == initialSci.yylD(initialRegion.CorEle(iq)));
        ip = intersect(ipx,ipy);
        subEle(iq) = ip;
        subCorMask(ip) = 1;
    end

    subEle = transpose(subEle);

end


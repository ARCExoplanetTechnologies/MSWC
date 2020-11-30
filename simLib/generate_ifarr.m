function [ifArr] = generate_ifarr(Ndm,pupil,sigma)


InfFuncSigma = (pupil.D/Ndm*sigma);
infdx = zeros(Ndm,length(pupil.x));

for q = 1:Ndm
    x_cent = q-Ndm/2-1/2;
    infdx(q,:) = exp(-4*log(2)*((pupil.x-pupil.D*x_cent/(Ndm-1)).^2)./(InfFuncSigma)^2);
end
InfFunc = exp(-4*log(2)*(pupil.xx.^2 + pupil.yy.^2)./(InfFuncSigma)^2);

ifArr.if = InfFunc;
ifArr.ifSigma = InfFuncSigma;
ifArr.ifdx = infdx;

end


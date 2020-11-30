function regions = initialize_regions(imageOut, params)

regionType = params.regionType; 
ScoreMaskReg = params.ScoreMaskReg;
IWZScoreReg = params.IWZScoreReg;
OWZScoreReg = params.OWZScoreReg;
angularSize = params.CorAngularSize; % angular extent in deg. of dark hole if annular
shiftReg = params.shiftReg;
CorMaskReg = params.CorMaskReg;


disp('Initializing Control Regions')

% Create the masks for scoring and correcting
if strcmp(regionType, 'rectangular')
    %CorMask = (USLAMD>=CorMaskReg(1) & USLAMD<=CorMaskReg(2) & VSLAMD>=CorMaskReg(3) & VSLAMD<=CorMaskReg(4));
    CorMask = (imageOut.xxlD >= CorMaskReg(1) & imageOut.xxlD <= CorMaskReg(2) & ...
        imageOut.yylD >= CorMaskReg(3) & imageOut.yylD <= CorMaskReg(4));
    
    ScoreMask = (imageOut.xxlD >= ScoreMaskReg(1) & imageOut.xxlD <= ScoreMaskReg(2) & ...
        imageOut.yylD >= ScoreMaskReg(3) & imageOut.yylD <= ScoreMaskReg(4));
    
    IWZscore = (imageOut.xxlD >= IWZScoreReg(1) & imageOut.xxlD <= IWZScoreReg(2) & ...
        imageOut.yylD >= IWZScoreReg(3) & imageOut.yylD <= IWZScoreReg(4));

    OWZscore = (imageOut.xxlD >= OWZScoreReg(1) & imageOut.xxlD <= OWZScoreReg(2) & ...
        imageOut.yylD >= OWZScoreReg(3) & imageOut.yylD <= OWZScoreReg(4));
    
    pts_per_side = 2;
    Score_perimeter_x = [linspace(ScoreMaskReg(1),ScoreMaskReg(2),pts_per_side) linspace(ScoreMaskReg(2),ScoreMaskReg(2),pts_per_side) linspace(ScoreMaskReg(2),ScoreMaskReg(1),pts_per_side) linspace(ScoreMaskReg(1),ScoreMaskReg(1),pts_per_side)];
    Score_perimeter_y = [linspace(ScoreMaskReg(3),ScoreMaskReg(3),pts_per_side) linspace(ScoreMaskReg(3),ScoreMaskReg(4),pts_per_side) linspace(ScoreMaskReg(4),ScoreMaskReg(4),pts_per_side) linspace(ScoreMaskReg(4),ScoreMaskReg(3),pts_per_side)];

    
    IWZ_perimeter_x = [linspace(IWZScoreReg(1),IWZScoreReg(2),pts_per_side) linspace(IWZScoreReg(2),IWZScoreReg(2),pts_per_side) linspace(IWZScoreReg(2),IWZScoreReg(1),pts_per_side) linspace(IWZScoreReg(1),IWZScoreReg(1),pts_per_side)];
    IWZ_perimeter_y = [linspace(IWZScoreReg(3),IWZScoreReg(3),pts_per_side) linspace(IWZScoreReg(3),IWZScoreReg(4),pts_per_side) linspace(IWZScoreReg(4),IWZScoreReg(4),pts_per_side) linspace(IWZScoreReg(4),IWZScoreReg(3),pts_per_side)];

    OWZ_perimeter_x = [linspace(OWZScoreReg(1),OWZScoreReg(2),pts_per_side) linspace(OWZScoreReg(2),OWZScoreReg(2),pts_per_side) linspace(OWZScoreReg(2),OWZScoreReg(1),pts_per_side) linspace(OWZScoreReg(1),OWZScoreReg(1),pts_per_side)];
    OWZ_perimeter_y = [linspace(OWZScoreReg(3),OWZScoreReg(3),pts_per_side) linspace(OWZScoreReg(3),OWZScoreReg(4),pts_per_side) linspace(OWZScoreReg(4),OWZScoreReg(4),pts_per_side) linspace(OWZScoreReg(4),OWZScoreReg(3),pts_per_side)];
   
    Cor_perimeter_x = [linspace(CorMaskReg(1),CorMaskReg(2),pts_per_side) linspace(CorMaskReg(2),CorMaskReg(2),pts_per_side) linspace(CorMaskReg(2),CorMaskReg(1),pts_per_side) linspace(CorMaskReg(1),CorMaskReg(1),pts_per_side)];
    Cor_perimeter_y = [linspace(CorMaskReg(3),CorMaskReg(3),pts_per_side) linspace(CorMaskReg(3),CorMaskReg(4),pts_per_side) linspace(CorMaskReg(4),CorMaskReg(4),pts_per_side) linspace(CorMaskReg(4),CorMaskReg(3),pts_per_side)];

elseif strcmp(regionType, 'shiftAnnular')
    
    [tLamD, rLamD] = cart2pol(imageOut.xxlD, imageOut.yylD);
      
    CorMask = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= CorMaskReg(1) & rLamD < CorMaskReg(2) & imageOut.xxlD >= shiftReg);
    
    ScoreMask = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= ScoreMaskReg(1) & rLamD < ScoreMaskReg(2) & imageOut.xxlD >= shiftReg);

    IWZscore = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= IWZScoreReg(1) & rLamD < IWZScoreReg(2) & imageOut.xxlD >= shiftReg);

    OWZscore = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= OWZScoreReg(1) & rLamD < OWZScoreReg(2) & imageOut.xxlD >= shiftReg);

    th_center = 0;
    th_range = deg2rad(angularSize);
    IWA = ScoreMaskReg(1);
    OWA = ScoreMaskReg(2);
    
    th_range_top = th_range - 2*asin(shiftReg/OWA);
    th_range_bottom = th_range - 2*asin(shiftReg/IWA);
    
    pts_per_side = 2000;

    %DZ_perimeter_th = [(th_center - (th_range/2)):0.01:(th_center + (th_range/2)) (th_center + (th_range/2)):-0.01:(th_center - (th_range/2)) th_center - (th_range/2)];
%     DZ_perimeter_th_bottom =  [(th_center - (th_range_bottom/2)):0.01:(th_center + (th_range_bottom/2))];
%     DZ_perimeter_th_top = [(th_center + (th_range_top/2)):-0.01:(th_center - (th_range_top/2))];
    DZ_perimeter_th_bottom =  linspace((th_center - (th_range_bottom/2)), (th_center + (th_range_bottom/2)), pts_per_side);
    DZ_perimeter_th_top = linspace((th_center + (th_range_top/2)), (th_center - (th_range_top/2)), pts_per_side);
    
    DZ_perimeter_th = [DZ_perimeter_th_bottom DZ_perimeter_th_top th_center - (th_range_bottom/2)];
    DZ_perimeter_r = [IWA.*ones(1,(length(DZ_perimeter_th_bottom))) OWA.*ones(1, (length(DZ_perimeter_th_top))) IWA];
    
    %DZ_perimeter_th = [(th_center - (th_range_bottom/2)):0.01:(th_center + (th_range_bottom/2)) (th_center + (th_range_top/2)):-0.01:(th_center - (th_range_top/2))  th_center - (th_range_bottom/2)];
    %DZ_[erimeter_th = [(th_center - th_range_bottom/2)
    %DZ_perimeter_r = [IWA.*ones(1,floor(length(DZ_perimeter_th)/2)) OWA.*ones(1, floor(length(DZ_perimeter_th)/2)) IWA];
    %DZ_perimeter_th = [DZ_perimeter_th th_center - (th_range_bottom/2)];

    Score_perimeter_x = DZ_perimeter_r.*cos(DZ_perimeter_th);
    Score_perimeter_y = DZ_perimeter_r.*sin(DZ_perimeter_th);
    
   	th_center = 0;
    th_range = deg2rad(angularSize);
    IWA = CorMaskReg(1);
    OWA = CorMaskReg(2);
    
    th_range_top = th_range - 2*asin(shiftReg/OWA);
    th_range_bottom = th_range - 2*asin(shiftReg/IWA);
    
    DZ_perimeter_th_bottom =  linspace((th_center - (th_range_bottom/2)), (th_center + (th_range_bottom/2)), pts_per_side);
    DZ_perimeter_th_top = linspace((th_center + (th_range_top/2)), (th_center - (th_range_top/2)), pts_per_side);
    %DZ_perimeter_th_top = [(th_center + (th_range_top/2)):-0.001:(th_center - (th_range_top/2))];
    
    DZ_perimeter_th = [DZ_perimeter_th_bottom DZ_perimeter_th_top th_center - (th_range_bottom/2)];
    DZ_perimeter_r = [IWA.*ones(1,(length(DZ_perimeter_th_bottom))) OWA.*ones(1, (length(DZ_perimeter_th_top))) IWA];
    
    %DZ_perimeter_th = [(th_center - (th_range_bottom/2)):0.01:(th_center + (th_range_bottom/2)) (th_center + (th_range_top/2)):-0.01:(th_center - (th_range_top/2))  th_center - (th_range_bottom/2)];
    %DZ_[erimeter_th = [(th_center - th_range_bottom/2)
    %DZ_perimeter_r = [IWA.*ones(1,floor(length(DZ_perimeter_th)/2)) OWA.*ones(1, floor(length(DZ_perimeter_th)/2)) IWA];
    %DZ_perimeter_th = [DZ_perimeter_th th_center - (th_range_bottom/2)];

    Cor_perimeter_x = DZ_perimeter_r.*cos(DZ_perimeter_th);
    Cor_perimeter_y = DZ_perimeter_r.*sin(DZ_perimeter_th);
       
    IWA = IWZScoreReg(1);
    OWA = IWZScoreReg(2);
    
    th_center = 0;
    th_range = deg2rad(angularSize);
    
    th_range_top = th_range - 2*asin(shiftReg/OWA);
    th_range_bottom = th_range - 2*asin(shiftReg/IWA);
    
    DZ_perimeter_th_bottom =  linspace((th_center - (th_range_bottom/2)), (th_center + (th_range_bottom/2)), pts_per_side);
    DZ_perimeter_th_top = linspace((th_center + (th_range_top/2)), (th_center - (th_range_top/2)), pts_per_side);
    %DZ_perimeter_th_top = [(th_center + (th_range_top/2)):-0.001:(th_center - (th_range_top/2))];
    
    DZ_perimeter_th = [DZ_perimeter_th_bottom DZ_perimeter_th_top th_center - (th_range_bottom/2)];
    DZ_perimeter_r = [IWA.*ones(1,(length(DZ_perimeter_th_bottom))) OWA.*ones(1, (length(DZ_perimeter_th_top))) IWA];
    
    %DZ_perimeter_th = [(th_center - (th_range_bottom/2)):0.01:(th_center + (th_range_bottom/2)) (th_center + (th_range_top/2)):-0.01:(th_center - (th_range_top/2))  th_center - (th_range_bottom/2)];
    %DZ_[erimeter_th = [(th_center - th_range_bottom/2)
    %DZ_perimeter_r = [IWA.*ones(1,floor(length(DZ_perimeter_th)/2)) OWA.*ones(1, floor(length(DZ_perimeter_th)/2)) IWA];
    %DZ_perimeter_th = [DZ_perimeter_th th_center - (th_range_bottom/2)];

    IWZ_perimeter_x = DZ_perimeter_r.*cos(DZ_perimeter_th);
    IWZ_perimeter_y = DZ_perimeter_r.*sin(DZ_perimeter_th);
    
    IWA = OWZScoreReg(1);
    OWA = OWZScoreReg(2);
    
    th_center = 0;
    th_range = deg2rad(angularSize);
    
    th_range_top = th_range - 2*asin(shiftReg/OWA);
    th_range_bottom = th_range - 2*asin(shiftReg/IWA);
    
    DZ_perimeter_th_bottom =  linspace((th_center - (th_range_bottom/2)), (th_center + (th_range_bottom/2)), pts_per_side);
    DZ_perimeter_th_top = linspace((th_center + (th_range_top/2)), (th_center - (th_range_top/2)), pts_per_side);
    %DZ_perimeter_th_top = [(th_center + (th_range_top/2)):-0.001:(th_center - (th_range_top/2))];
    
    DZ_perimeter_th = [DZ_perimeter_th_bottom DZ_perimeter_th_top th_center - (th_range_bottom/2)];
    DZ_perimeter_r = [IWA.*ones(1,(length(DZ_perimeter_th_bottom))) OWA.*ones(1, (length(DZ_perimeter_th_top))) IWA];
    
    %DZ_perimeter_th = [(th_center - (th_range_bottom/2)):0.01:(th_center + (th_range_bottom/2)) (th_center + (th_range_top/2)):-0.01:(th_center - (th_range_top/2))  th_center - (th_range_bottom/2)];
    %DZ_[erimeter_th = [(th_center - th_range_bottom/2)
    %DZ_perimeter_r = [IWA.*ones(1,floor(length(DZ_perimeter_th)/2)) OWA.*ones(1, floor(length(DZ_perimeter_th)/2)) IWA];
    %DZ_perimeter_th = [DZ_perimeter_th th_center - (th_range_bottom/2)];

    OWZ_perimeter_x = DZ_perimeter_r.*cos(DZ_perimeter_th);
    OWZ_perimeter_y = DZ_perimeter_r.*sin(DZ_perimeter_th);
    
elseif strcmp(regionType, 'annular')
    [tLamD, rLamD] = cart2pol(imageOut.xxlD, imageOut.yylD);
      
    CorMask = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= CorMaskReg(1) & rLamD < CorMaskReg(2));
    
    ScoreMask = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= ScoreMaskReg(1) & rLamD < ScoreMaskReg(2));
    %ProbeReg = [min(min(USLAMD(corInds))) max(max(USLAMD(corInds))) min(min(VSLAMD(corInds))) max(max(VSLAMD((corInds))))];
    
    IWZscore = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= IWZScoreReg(1) & rLamD < IWZScoreReg(2));

    OWZscore = (tLamD <= deg2rad(angularSize/2) & tLamD >= -deg2rad(angularSize/2) & ...
        rLamD >= OWZScoreReg(1) & rLamD < OWZScoreReg(2));
    
    th_center = 0;
    th_range = deg2rad(angularSize);
    IWA = ScoreMaskReg(1);
    OWA = ScoreMaskReg(2);
    
    DZ_perimeter_th = [(th_center - (th_range/2)):0.001:(th_center + (th_range/2)) (th_center + (th_range/2)):-0.001:(th_center - (th_range/2)) th_center - (th_range/2)];
    DZ_perimeter_r = [IWA.*ones(1,floor(length(DZ_perimeter_th)/2)) OWA.*ones(1, floor(length(DZ_perimeter_th)/2)) IWA];

    Score_perimeter_x = DZ_perimeter_r.*cos(DZ_perimeter_th);
    Score_perimeter_y = DZ_perimeter_r.*sin(DZ_perimeter_th);
    
   	th_center = 0;
    th_range = deg2rad(angularSize);
    IWA = CorMaskReg(1);
    OWA = CorMaskReg(2);

    Cor_perimeter_x = DZ_perimeter_r.*cos(DZ_perimeter_th);
    Cor_perimeter_y = DZ_perimeter_r.*sin(DZ_perimeter_th);
end

%corInds = find(CorMask == 1);

DZcol = imageOut.N;
DZrow = imageOut.N;

area = sum(sum(ScoreMask));
Maskline = reshape(ScoreMask,1,DZcol*DZrow);
cor_ele = find(CorMask~=0);
Maskline_cor_ele = Maskline(cor_ele);
score_ele = find(ScoreMask~=0);
iwz_score_ele = find(IWZscore~=0);
owz_score_ele = find(OWZscore~=0);

Wmat = CorMask;
W = diag(Wmat(cor_ele));

regions.xlD = imageOut.xlD;
regions.ylD = imageOut.ylD;
regions.regionType = regionType;
regions.CorMask = CorMask;
regions.ScoreMask = ScoreMask;
regions.IWZscore = IWZscore;
regions.OWZscore = OWZscore;
%regions.CorInds = corInds;
regions.CorEle = cor_ele;
regions.MasklineCorEle = Maskline_cor_ele;
regions.score_ele = score_ele;
regions.iwz_score_ele = iwz_score_ele;
regions.owz_score_ele = owz_score_ele;
regions.Score_perimeter_x = Score_perimeter_x;
regions.Score_perimeter_y = Score_perimeter_y;
regions.Cor_perimeter_x = Cor_perimeter_x;
regions.Cor_perimeter_y = Cor_perimeter_y;
regions.IWZ_perimeter_x = IWZ_perimeter_x;
regions.IWZ_perimeter_y = IWZ_perimeter_y;
regions.OWZ_perimeter_x = OWZ_perimeter_x;
regions.OWZ_perimeter_y = OWZ_perimeter_y;
regions.Wmat = Wmat; % feathermap
regions.W = diag(Wmat(cor_ele));




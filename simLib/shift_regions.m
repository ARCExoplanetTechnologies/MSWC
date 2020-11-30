function [shiftedRegion] = shift_regions(imageOut, initialRegion, shiftX, shiftY)
    
    shiftXpix = ceil(shiftX * imageOut.samplesPerflD);
    shiftYpix = ceil(shiftY * imageOut.samplesPerflD);
    
    shiftCorMask = circshift(initialRegion.CorMask,shiftXpix,2);
    shiftCorMask = circshift(shiftCorMask,shiftYpix,1);
    
    shiftCorEle = find(shiftCorMask ~= 0);

	shiftScoreMask = circshift(initialRegion.ScoreMask,shiftXpix,2);
    shiftScoreMask = circshift(shiftScoreMask,shiftYpix,1);
    
    shiftScoreEle = find(shiftScoreMask ~= 0);
    
    shiftedRegion = initialRegion;
    shiftedRegion.CorMask = shiftCorMask;
    shiftedRegion.ScoreMask = shiftScoreMask;
    shiftedRegion.CorEle = shiftCorEle;
    shiftedRegion.ScoreEle = shiftScoreEle;
    
    %clear shiftedRegion.IWZ_perimeter_x;% shiftedRegion.IWZ_perimeter_y shiftedRegion.OWZ_perimeter_x shiftedRegion.OWZ_perimeter_y shiftedRegion.Score_perimeter_x shiftedRegion.Score_perimeter_y shiftedRegion.Cor_perimeter_x shiftedRegion.Cor_perimeter_x
    fields = {'IWZ_perimeter_x', 'IWZ_perimeter_y', 'OWZ_perimeter_x', 'OWZ_perimeter_y', 'Cor_perimeter_x', 'Cor_perimeter_y', 'Score_perimeter_x', 'Score_perimeter_y', 'score_ele'};
    shiftedRegion = rmfield(shiftedRegion, fields);

end
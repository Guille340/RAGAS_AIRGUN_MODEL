function isTrue = isBubbleStruct(S)

Bubble = initialiseBubbleStruct();
isTrueCon = false;
isTrueArr = false;
isTrueCal = false;
isTrueInt = false;
isTrueDyn = false;

if isstruct(S)
    % Layer 1
    fieldNamesBub_valid = fieldnames(Bubble);
    fieldNamesBub = fieldnames(S(1));
    
    if all(ismember(fieldNamesBub,fieldNamesBub_valid)) ...
            && all(ismember(fieldNamesBub_valid,fieldNamesBub))
        
        % Layer 2
        fieldNamesCon = fieldnames(S(1).Config);
        fieldNamesCon_valid = fieldnames(Bubble.Config);
        if all(ismember(fieldNamesCon,fieldNamesCon_valid)) ...
            && all(ismember(fieldNamesCon_valid,fieldNamesCon))
            isTrueCon = true;
        end
        
        fieldNamesGun = fieldnames(S(1).Array);
        fieldNamesGun_valid = fieldnames(Bubble.Array);
        if all(ismember(fieldNamesGun,fieldNamesGun_valid)) ...
            && all(ismember(fieldNamesGun_valid,fieldNamesGun))
            isTrueArr = true;
        end
        
        fieldNamesCal = fieldnames(S(1).Calibration);
        fieldNamesCal_valid = fieldnames(Bubble.Calibration);
        if all(ismember(fieldNamesCal,fieldNamesCal_valid)) ...
            && all(ismember(fieldNamesCal_valid,fieldNamesCal))
            isTrueCal = true;
        end
        
        fieldNamesInt = fieldnames(S(1).Interactions);
        fieldNamesInt_valid = fieldnames(Bubble.Interactions);
        if all(ismember(fieldNamesInt,fieldNamesInt_valid)) ...
            && all(ismember(fieldNamesInt_valid,fieldNamesInt))
            isTrueInt = true;
        end
        
        fieldNamesDyn = fieldnames(S(1).Dynamics);
        fieldNamesDyn_valid = fieldnames(Bubble.Dynamics);
        if all(ismember(fieldNamesDyn,fieldNamesDyn_valid)) ...
            && all(ismember(fieldNamesDyn_valid,fieldNamesDyn))
            isTrueDyn = true;
        end
    end
end

isTrue = isTrueCon && isTrueArr && isTrueCal && isTrueInt && isTrueDyn;

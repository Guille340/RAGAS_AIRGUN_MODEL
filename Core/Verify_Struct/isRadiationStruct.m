function isTrue = isRadiationStruct(S)

Radiation = initialiseRadiationStruct();
isTrueCon = false;
isTrueArr = false;
isTrueSig = false;

if isstruct(S)
    % Layer 1
    fieldNamesRad_valid = fieldnames(Radiation);
    fieldNamesRad = fieldnames(S(1));
    
    if all(ismember(fieldNamesRad,fieldNamesRad_valid)) ...
            && all(ismember(fieldNamesRad_valid,fieldNamesRad))
        
        % Layer 2
        fieldNamesCon = fieldnames(S(1).Config);
        fieldNamesCon_valid = fieldnames(Radiation.Config);
        if all(ismember(fieldNamesCon,fieldNamesCon_valid)) ...
            && all(ismember(fieldNamesCon_valid,fieldNamesCon))
            isTrueCon = true;
        end
        
        fieldNamesArr = fieldnames(S(1).Array);
        fieldNamesArr_valid = fieldnames(Radiation.Array);
        if all(ismember(fieldNamesArr,fieldNamesArr_valid)) ...
            && all(ismember(fieldNamesArr_valid,fieldNamesArr))
            isTrueArr = true;
        end
        
        fieldNamesSig = fieldnames(S(1).Signature);
        fieldNamesSig_valid = fieldnames(Radiation.Signature);
        if all(ismember(fieldNamesSig,fieldNamesSig_valid)) ...
            && all(ismember(fieldNamesSig_valid,fieldNamesSig))
            isTrueSig = true;
        end
        
        fieldNamesMet = fieldnames(S(1).Metrics);
        fieldNamesMet_valid = fieldnames(Radiation.Metrics);
        if all(ismember(fieldNamesMet,fieldNamesMet_valid)) ...
            && all(ismember(fieldNamesMet_valid,fieldNamesMet))
            isTrueMet = true;
        end
    end
end

isTrue = isTrueCon && isTrueArr && isTrueSig && isTrueMet;

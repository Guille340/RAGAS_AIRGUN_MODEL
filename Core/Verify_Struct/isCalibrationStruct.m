function isTrue = isCalibrationStruct(S)

isTrue = false;
if isstruct(S)
    Calibration = initialiseCalibrationStruct();
    fieldNames_valid = fieldnames(Calibration);
    fieldNames = fieldnames(S);

    if all(ismember(fieldNames,fieldNames_valid)) ...
            && all(ismember(fieldNames_valid,fieldNames))
        isTrue = true;
    end
end

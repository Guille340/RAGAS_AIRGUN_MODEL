function isTrue = isConfigStruct(S)

Config = initialiseConfigStruct();
fieldNames_valid = fieldnames(Config);
fieldNames = fieldnames(S);

isTrue = false;
if all(ismember(fieldNames,fieldNames_valid)) ...
        && all(ismember(fieldNames_valid,fieldNames))
    isTrue = true;
end


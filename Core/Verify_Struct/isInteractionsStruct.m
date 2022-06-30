function isTrue = isInteractionsStruct(S)

isTrue = false;
if isstruct(S)
    Interactions = initialiseInteractionsStruct();
    fieldNames_valid = fieldnames(Interactions);
    fieldNames = fieldnames(S);

    if all(ismember(fieldNames,fieldNames_valid)) ...
            && all(ismember(fieldNames_valid,fieldNames))
        isTrue = true;
    end
end


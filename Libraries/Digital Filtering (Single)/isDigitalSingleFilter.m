function isTrue = isDigitalSingleFilter(DigitalFilter)

isTrue = false; % initialise output to FALSE

% Verify Field Names in Layer 1
DigitalFilterFields_valid = {'sampleRate','freqCutoff','decimator','target'};
DigitalFilterFields = fieldnames(DigitalFilter(1));
flag_layer1 = all(ismember(DigitalFilterFields_valid,DigitalFilterFields)) ...
    && all(ismember(DigitalFilterFields,DigitalFilterFields_valid));

% Verify Field Names in Layer 2
if flag_layer1
    % Verify Field Names in Field 1 of Layer 2
    DecimatorFields_valid = {'filterType','filterResponse','filterOrder',...
        'filterObject','halfPowerFreqn','decimationFactor','sampleRate',...
        'halfPowerFreq'};
    
    DecimatorFields = fieldnames(DigitalFilter(1).decimator);
    flag1_layer2 = all(ismember(DecimatorFields_valid,...
        DecimatorFields)) && all(ismember(DecimatorFields,...
        DecimatorFields_valid));
    
    TargetFields_valid = {'filterType','filterResponse','filterOrder',...
        'filterObject','halfPowerFreqn1','halfPowerFreqn2','sampleRate',...
        'halfPowerFreq1','halfPowerFreq2'};
    TargetFields = fieldnames(DigitalFilter(1).target);
    flag2_layer2 = all(ismember(TargetFields_valid,...
        TargetFields)) && all(ismember(TargetFields,...
        TargetFields_valid));
    
    % Output is TRUE if All Flags are True
    if flag1_layer2 && flag2_layer2
        isTrue = true;
    end        
end

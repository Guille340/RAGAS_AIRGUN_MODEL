function isTrue = isCoalescent(Bubble)

isTrue = false;
gunIds = [Bubble.Array.gunId];
nGuns = length(gunIds);

for q = 1:nGuns
    % Load Parameters
    gunId_rec = gunIds(q); % id of airgun affected by surrounding airguns
    gunIds_sou = gunIds(gunIds ~= gunId_rec); % id of airguns interacting with current airgun
    
    % Receiver Position
    iGun_rec = q;
    x_rec = Bubble.Array(iGun_rec).gunX; % x position of emitting air gun [m]
    y_rec = Bubble.Array(iGun_rec).gunY; % y position of emitting air gun [m]
    z_rec = Bubble.Array(iGun_rec).gunZ; % z position of emitting air gun [m]

    % Check Bubble Coalescence 
    nGuns_sou = length(gunIds_sou);
    for p = 1:nGuns_sou        
        % Distance between the Two Bubbles
        iGun_sou = find(gunIds == gunIds_sou(p));
        x_sou = Bubble.Array(iGun_sou).gunX; % x position of emitting air gun [m]
        y_sou = Bubble.Array(iGun_sou).gunY; % y position of emitting air gun [m]
        z_sou = Bubble.Array(iGun_sou).gunZ; % z position of emitting air gun [m]
        Rd = sqrt((x_rec - x_sou)^2 + (y_rec - y_sou)^2 + (z_rec - z_sou)^2); % bubble distance [m]

        % Minimum Distance between the Two Bubbles to Avoid Coalescense
        rb_rec = Bubble.Dynamics(iGun_rec).bubbleRadius;
        rb_sou = Bubble.Dynamics(iGun_sou).bubbleRadius;
        Rd_min = max(rb_rec + rb_sou);
        
        % Coalescence Check
        if Rd < Rd_min
            isTrue = true;
        end
    end
end
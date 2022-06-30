
function km = tconductivitymix(ka,kv,mua,muv,Ma,Mv,ma,mv)

%**************************************************************************
%  km = tconductivitymix(ka,kv,mua,muv,Ma,Mv,ma,mv)
%
%  DESCRIPTION: Calculates the thermal conductivity of air or water vapour 
%
%  INPUT VARIABLES
%  - ka: thermal conductivity of dry air [W m-1 K-1]
%  - kv: thermal conductivity of water vapour [W m-1 K-1]
%  - mua: dynamic viscosity of dry air [Pa s]
%  - muv: dynamic viscosity of water vapour [Pa s]
%  - Ma: molecular weight of dry air [kg mol-1]
%  - Mv: molecular weight of water [kg mol-1]
%  - ma: dynamic viscosity of dry air [Pa s]
%  - mv: dynamic viscosity of water vapour [Pa s]
%
%  OUTPUT VARIABLES
%  - km: thermal conductivity of the gas mixture (dry air + water vapour)
%    [W m-1 K-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  FUNCTION CALLS
%  1) km = tconductivitymix(ka,kv,mua,muv,Ma,Mv,ma,mv)
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  26 Sep 2016
%
%**************************************************************************

Aav = (1 + (mua/muv)^0.5*(Mv/Ma)^0.25)^2/(8*(1+Ma/Mv))^0.5;
Ava = (1 + (muv/mua)^0.5*(Ma/Mv)^0.25)^2/(8*(1+Mv/Ma))^0.5;
xa = ma / (ma + mv); % molar fraction of dry air
xv = mv / (ma + mv); % molar fraction of water vapour
kpa = ka/(1 + Aav*xv/xa); % partial thermal conductivity of dry air [W m-1 K-1]
kpv = kv/(1 + Ava*xa/xv); % partial thermal conductivity of water vapor [W m-1 K-1]
km = kpa + kpv; % thermal conductivity of the gas mixture [W m-1 K-1]
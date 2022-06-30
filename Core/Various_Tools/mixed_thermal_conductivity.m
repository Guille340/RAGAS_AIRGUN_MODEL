

T = 400; % temperature [K]
rhov = 1.108; % density [kg m-3]
ma = 0.2; % mass of dry air [kg]
mv = 0.05; % mass of water vapour [kg]
Ma = 0.028970285; % molecular weight of dry air [kg mol-1]
Mv = 0.01801528; % molecular weight of water [kg mol-1]
ka = tconductivity('air',T); % thermal conductivity of air [W m-1 K-1]
kv = tconductivity('h2o',T,rhov); % thermal conductivity of water vapour [W m-1 K-1]
mua = viscosity('air',T);
muv = viscosity('h2o',T,rhov);
km = tconductivitymix(ka,kv,mua,muv,Ma,Mv,ma,mv);
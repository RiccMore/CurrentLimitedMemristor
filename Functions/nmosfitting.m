% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: NMOS fitting function
%
% Args:
% - vGS: gate source voltage [V]
% - vSB: source body voltage [V]
% - vDS: drain source voltage [V]
% - Vtn0: zero bias forward threshold voltage [V]
% - gamma: body effect parameter [sqrt(V)]
% - phi: surface potential [V]
% - Kn: transconductance gain [A/V^2]
% - lambda: channel modulation factor [1/V]
%
% Returns:
% - iD: drain current [A]

function iD = nmosfitting(vGS,vSB,vDS,Vtn0,gamma,phi,Kn,lambda)
    iD = zeros(size(vDS));
    Vtn = Vtn0+gamma*(sqrt(vSB+phi)-sqrt(phi));
    triode = Vtn < vGS & vDS < vGS-Vtn;
    saturation = Vtn < vGS & vDS >= vGS-Vtn;
    iD(triode) = Kn*(1+lambda*vDS(triode)).* ...
                        ((vGS-Vtn(triode)).*vDS(triode)- ...
                        (vDS(triode).^2)/2);
    iD(saturation) = (Kn/2)*(1+lambda*vDS(saturation)).* ...
                     (vGS-Vtn(saturation)).^2;
end
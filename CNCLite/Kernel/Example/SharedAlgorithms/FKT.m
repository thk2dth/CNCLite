function [ wc ] = FKT( mc, mp )
% forward kinematic transformation.
% Input:
%   mc, machine coordinate, [X; Y; Z; A; C]. A, C in rad.
%   mp, machine property, [mx; my; mz]. in mm.
% Output:
%   wc, workpiece coordinate, [x; y; z; i; j; k]
sa = sin( mc(4) );
ca = cos( mc(4) );
sc = sin( mc(5) );
cc = cos( mc(5) );
T = [cc, ca*sc, sa*sc;...
    sc, -ca*cc, -sa*cc;...
    0, -sa, ca];
d = [mc(1)-mp(1); mc(2)+mp(2); mc(3)-mp(3)];
pos = mp + T * d;
ori = [sa*sc; -sa*cc; ca];
wc = [pos; ori];

end


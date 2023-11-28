function stresses = alternative_hertz_stress(r, z, a, nu, p0)

R = abs(r/a);
Z = abs(z/a);

%% Maximum contact pressure
s_0 = 3/2*p0;

%% L>=0
%From implicit e.q. 6 we define L

L = sqrt(0.5*(R^2+Z^2-1+sqrt((R^2+Z^2-1)^2+4*Z^2)));
%(0.5000*((R^2 + Z^2 - 1)^2 + 4*Z^2)^(1/2)

%% The full stress fields normalised by maximum pressure are then given by:

if Z > 0
s_rr_s0 = - L^3*R^2*Z/(L^4+Z^2)/(1+L^2)^2 ...
          - (1-2*nu)*(Z/L/(1+L^2)-1/3/R^2*(1-Z^3/L^3))...
          + Z/L*(L*(1+nu)*atan(1/L)-(1-nu)*L^2/(1+L^2)-2*nu)    ;
s_tt_s0 = -(1-2*nu)/3/R^2*(1-(Z^3/L^3)) ...
          + Z/L*(L*(1+nu)*atan(1/L)-(1-nu)*L^2/(1+L^2)-2*nu)    ;
s_zz_s0 = -Z^3/(L*(L^4+Z^2))                                    ;
s_zr_s0 = L*R*Z^2/(L^4+Z^2)/(1+L^2)                             ;
s_rt_s0 = 0                                                     ;
s_zt_s0 = 0                                                     ;

elseif R<=1 && Z==0
    s_rr_s0 = (1-2*nu)*1/3/R^2*(1-(1-R^2)^(3/2))...
              - (1-R^2)^(0.5)   ;
    
    s_tt_s0 = -(1-2*nu)*1/3/R^2*(1-(1-R^2)^(3/2))...
              - 2*nu*(1-R^2)^0.5;
    s_zz_s0 = -(1-R^2)^0.5      ;
    s_zr_s0 = 0                 ;
    s_rt_s0 = 0                 ;
    s_zt_s0 = 0                 ;
elseif R>1 && Z==0
    s_rr_s0 = (1-2*nu)*1/3/R^2  ;
    s_tt_s0 = -(1-2*nu)/3/R^2   ;
    s_zz_s0 = 0                 ;
    s_zr_s0 = 0                 ;
    s_rt_s0 = 0                 ;
    s_zt_s0 = 0                 ;
else
    error('Condition not found')
end


if R==0
    s_rr_s0 = (1+nu)*Z*atan(1/Z)+(1-2*(1+nu)*(1+Z^2))/2/(1+Z^2) ;

    s_tt_s0 = -(1-2*nu)*0.5*1/(1+Z^2) ...
          +   (Z*(1+nu)*atan(1/Z)-(1-nu)*Z^2/(1+Z^2)-2*nu)      ;   

    s_zz_s0 = -1/(1+Z^2)                                        ;
    s_zr_s0 = 0                                                 ;
    s_rt_s0 = 0                                                 ;
    s_zt_s0 = 0                                                 ;
    
end

stresses = [s_rr_s0   s_rt_s0   s_zr_s0 ;
            s_rt_s0   s_tt_s0   s_zt_s0 ;
            s_zr_s0   s_zt_s0   s_zz_s0];
stresses = stresses*s_0;

end
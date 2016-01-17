function [Bx, By, Bz, wbvec_x, wbvec_y, wbvec_z] = single_mode_bending_energy(xclks,yclks,zclks)
%% Calculates the sum of the bending energies corresponding to individual modes for the three coordinates
global Y_LK_g Y_LK_phi Y_LK_theta Y_PP Y_TT Y_TP wp wt st ct stsq sp spsq cp cpsq ctsq ctcu ctqad r_o
global R
r_o = 1;
%%% calculate bending energy of the sphere
R(:)   =  r_o;Rp  = 0;Rt  = 0;Rpp = 0;Rtt = 0;Rtp = 0;
deltaR = Rtt + ct./st.*Rt + Rpp./stsq;divRsq = Rt.^2 + Rp.^2./stsq;
RDR = R.^2 + divRsq;T2 = R.*divRsq;T3 = Rt.^2.*Rtt + 2.*Rt.*Rp.*Rtp./stsq -Rt.*Rp.^2 .* ct./st.^3 + Rp.^2 .* Rpp./stsq.^2;
I1 = 2.*R -deltaR + (T2 + T3)./(RDR);Q = R.*(RDR).^(1/2);
Swb = 1/16/pi.*abs(sum(sum(wp.*wt.*(I1).^2./Q.*st)));

wb_vec = [];
for ix = 1:length(xclks)
    x = zeros(size(xclks));x(ix) = r_o/5;
    %%%%%%%% Calculate the surface/shape properties
    R(:)   =  r_o + Y_LK_g*x;Rp  = Y_LK_phi*x;Rt  = Y_LK_theta*x;Rpp = Y_PP*x;Rtt = Y_TT*x;Rtp = Y_TP*x;
    deltaR = Rtt + ct./st.*Rt + Rpp./stsq;
    divRsq = Rt.^2 + Rp.^2./stsq;
    RDR = R.^2 + divRsq;
    T2 = R.*divRsq;
    T3 = Rt.^2.*Rtt + 2.*Rt.*Rp.*Rtp./stsq -Rt.*Rp.^2 .* ct./st.^3 + Rp.^2 .* Rpp./stsq.^2;
    I1 = 2.*R -deltaR + (T2 + T3)./(RDR);
    Q = R.*(RDR).^(1/2);
    wb = 1/16/pi.*abs(sum(sum(wp.*wt.*(I1).^2./Q.*st)));
    wb_vec(ix) = wb;
end
wbvec_x = wb_vec-Swb;
Bx = sum(wbvec_x);
wb_vec = [];
for ix = 1:length(yclks)
    x = zeros(size(yclks));x(ix) = r_o/5;
    %%%%%%%% Calculate the surface/shape properties
    R(:)   =  r_o + Y_LK_g*x;Rp  = Y_LK_phi*x;Rt  = Y_LK_theta*x;Rpp = Y_PP*x;Rtt = Y_TT*x;Rtp = Y_TP*x;
    deltaR = Rtt + ct./st.*Rt + Rpp./stsq;
    divRsq = Rt.^2 + Rp.^2./stsq;
    RDR = R.^2 + divRsq;
    T2 = R.*divRsq;
    T3 = Rt.^2.*Rtt + 2.*Rt.*Rp.*Rtp./stsq -Rt.*Rp.^2 .* ct./st.^3 + Rp.^2 .* Rpp./stsq.^2;
    I1 = 2.*R -deltaR + (T2 + T3)./(RDR);
    Q = R.*(RDR).^(1/2);
    wb = 1/16/pi.*abs(sum(sum(wp.*wt.*(I1).^2./Q.*st)));
    wb_vec(ix) = wb;
end
wbvec_y = (wb_vec-Swb);
By = sum(wbvec_y);

wb_vec = [];
for ix = 1:length(zclks)
    x = zeros(size(zclks));x(ix) = r_o/5;
    %%%%%%%% Calculate the surface/shape properties
    R(:)   =  r_o + Y_LK_g*x;Rp  = Y_LK_phi*x;Rt  = Y_LK_theta*x;Rpp = Y_PP*x;Rtt = Y_TT*x;Rtp = Y_TP*x;
    %%% Calculate some expressions only once;
    deltaR = Rtt + ct./st.*Rt + Rpp./stsq;
    divRsq = Rt.^2 + Rp.^2./stsq;
    RDR = R.^2 + divRsq;
    T2 = R.*divRsq;
    T3 = Rt.^2.*Rtt + 2.*Rt.*Rp.*Rtp./stsq -Rt.*Rp.^2 .* ct./st.^3 + Rp.^2 .* Rpp./stsq.^2;
    I1 = 2.*R -deltaR + (T2 + T3)./(RDR);
    Q = R.*(RDR).^(1/2);
    wb = 1/16/pi.*abs(sum(sum(wp.*wt.*(I1).^2./Q.*st)));
    wb_vec(ix) = wb;
end
wbvec_z = (wb_vec-Swb);
Bz = sum(wbvec_z);


% for ix = 1:length(xclks)
%     x = zeros(size(xclks));x(ix) = xclks(ix);
%     %%%%%%%% Calculate the surface/shape properties
%     R(:)   =  r_o + Y_LK_g*x;Rp  = Y_LK_phi*x;Rt  = Y_LK_theta*x;Rpp = Y_PP*x;Rtt = Y_TT*x;Rtp = Y_TP*x;
%     %%% Calculate some expressions only once;
%     Rsq(:) = R.^2; Rcu(:) = R.^3;Rtsq(:) = Rt.^2; Rtcu(:) = Rt.^3;Rpsq(:) = Rp.^2; Rpcu(:) = Rp.^3;
%     blk_I(:) = (Rpsq+Rsq+Rtsq-Rtsq.*ctsq-Rsq.*ctsq);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     S(:) = R.*(blk_I).^(1./2);
%     E(:) = (Rtsq.*stsq+2.*Rt.*ct.*R.*st+Rsq.*ctsq).*spsq+Rtsq.*stsq.*cpsq+2.*Rt.*st.*cpsq.*R.*ct+Rsq.*ctsq.*cpsq-2.*Rt.*ct.*R.*st+Rsq.*stsq+Rtsq.*ctsq;
%     F(:) = Rp.*(stsq.*cpsq.*Rt+st.*cpsq.*R.*ct+stsq.*spsq.*Rt+st.*spsq.*R.*ct+ctsq.*Rt-ct.*R.*st);
%     G(:) = stsq.*Rpsq.*cpsq+stsq.*Rsq.*spsq+stsq.*Rpsq.*spsq+stsq.*Rsq.*cpsq+ctsq.*Rpsq;
%     LS(:) = (Rt.*st.*cp+R.*cp.*ct).*(Rtp.*sp.*Rsq-Rtt.*ct.*st.*cp.*Rpsq+Rtt.*ctcu.*st.*cp.*Rsq+Rsq.*cp.*Rt+2.*Rtcu.*ctqad.*cp+2.*Rtcu.*cp-Rtt.*ct.*st.*cp.*Rsq+Rtp.*sp.*Rtsq+Rcu.*cp.*ct.*st-Rcu.*cp.*ctcu.*st+2.*R.*cp.*Rt.*ctsq.*Rtt+Rsq.*cp.*ctqad.*Rt-Rp.*sp.*R.*Rt-2.*Rsq.*cp.*ctsq.*Rt+Rt.*ct.*st.*cp.*Rp.*Rtp+R.*cp.*ctsq.*Rp.*Rtp-Rp.*sp.*Rsq.*ct.*st-3.*Rt.*ctsq.*cp.*Rpsq-R.*cp.*Rt.*Rtt-2.*Rtsq.*ctcu.*st.*cp.*R-R.*cp.*ctqad.*Rt.*Rtt-R.*cp.*Rp.*Rtp-Rtp.*sp.*Rsq.*ctsq-Rtp.*sp.*Rtsq.*ctsq-Rp.*sp.*Rt.*Rtt+2.*R.*cp.*ct.*st.*Rpsq-Rp.*sp.*Rtsq.*ct.*st-4.*Rtcu.*cp.*ctsq+Rp.*sp.*Rt.*ctsq.*Rtt+2.*Rtsq.*ct.*st.*cp.*R+2.*Rt.*cp.*Rpsq+Rp.*sp.*R.*ctsq.*Rt)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(Rpsq+Rsq+Rtsq-Rtsq.*ctsq-Rsq.*ctsq).^(1./2)-(Rt.*st.*sp+R.*ct.*sp).*(-4.*Rtcu.*ctsq.*sp+2.*sp.*Rtcu.*ctqad+2.*sp.*Rt.*Rpsq+2.*R.*sp.*Rt.*ctsq.*Rtt-cp.*Rtp.*Rsq-cp.*Rtp.*Rtsq-Rtt.*ct.*st.*sp.*Rpsq+Rsq.*sp.*Rt+2.*Rtsq.*ct.*st.*sp.*R+cp.*Rtp.*Rtsq.*ctsq+Rtt.*ctcu.*st.*sp.*Rsq+Rt.*ct.*st.*sp.*Rp.*Rtp-Rcu.*sp.*ctcu.*st+Rp.*sp.*R.*ctsq.*Rtp-Rtt.*ct.*st.*sp.*Rsq+cp.*Rtp.*Rsq.*ctsq+Rp.*cp.*Rtsq.*ct.*st+Rp.*cp.*Rt.*Rtt-3.*sp.*Rt.*ctsq.*Rpsq-Rp.*sp.*R.*Rtp+2.*R.*ct.*sp.*st.*Rpsq-R.*sp.*Rt.*Rtt-2.*Rsq.*sp.*ctsq.*Rt+Rsq.*sp.*ctqad.*Rt+2.*sp.*Rtcu-Rp.*cp.*Rt.*ctsq.*Rtt+Rp.*cp.*R.*Rt+Rcu.*sp.*ct.*st-R.*cp.*ctsq.*Rp.*Rt-2.*Rtsq.*ctcu.*st.*sp.*R+Rp.*cp.*Rsq.*ct.*st-R.*sp.*ctqad.*Rt.*Rtt)./(blk_I).^(3./2)-(Rt.*ct-R.*st).*(Rcu-Rt.*Rsq.*ct.*st-Rt.*ctsq.*Rp.*Rtp+Rt.*ctcu.*Rsq.*st-3.*Rt.*ct.*st.*Rpsq+2.*Rtt.*Rsq.*ctsq+R.*st.*ct.*Rt.*Rtt-R.*st.*ctcu.*Rt.*Rtt+2.*Rtcu.*ctcu.*st+2.*Rtsq.*ctqad.*R+Rtt.*ctsq.*Rpsq-Rtt.*ctqad.*Rsq-2.*R.*ctsq.*Rpsq+2.*R.*Rtsq+R.*st.*ct.*Rp.*Rtp-Rtt.*Rpsq-Rtt.*Rsq-2.*Rcu.*ctsq+Rcu.*ctqad+R.*Rpsq+Rt.*Rp.*Rtp-2.*Rtcu.*ct.*st-4.*R.*ctsq.*Rtsq)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(blk_I).^(1./2);
%     M(:) = 1./2.*(-Rcu.*st.*Rtp-ctcu.*Rcu.*cpsq.*Rp-ctqad.*R.*st.*spsq.*Rp.*Rt.*Rtt+ctqad.*Rcu.*st.*cpsq.*Rtp+ct.^5.*R.*cpsq.*Rp.*Rtsq-ct.^5.*Rsq.*cpsq.*Rt.*Rtp+ctqad.*R.*Rtsq.*st.*Rtp-ctqad.*R.*Rtsq.*st.*cpsq.*Rtp-ctcu.*R.*cpsq.*Rpcu+ct.*Rsq.*Rt.*Rtp+ctcu.*Rsq.*Rt.*stsq.*Rtp.*cpsq-2.*R.*Rtsq.*st.*cpsq.*Rtp+ctcu.*Rsq.*Rp.*spsq.*Rtt.*stsq+2.*ctsq.*Rcu.*st.*Rtp+ctcu.*R.*stsq.*Rp.*Rtsq-ctcu.*Rcu.*spsq.*Rp.*stsq+2.*ctcu.*Rsq.*cpsq.*Rt.*Rtp-3.*ctcu.*R.*cpsq.*Rp.*Rtsq-ctcu.*Rsq.*stsq.*Rt.*Rtp+ctcu.*Rsq.*cpsq.*Rp.*Rpp-ctsq.*R.*Rtsq.*st.*Rtp+ct.*Rcu.*cpsq.*Rp-ct.*Rtsq.*Rp.*Rpp+3.*ctsq.*R.*Rtsq.*st.*cpsq.*Rtp+2.*ct.*R.*cpsq.*Rpcu+2.*ctsq.*R.*Rt.*st.*cpsq.*Rp.*Rpp-2.*ctsq.*R.*st.*Rt.*Rp.*Rpp-R.*Rt.*st.*cpsq.*Rp.*Rpp-ct.*Rsq.*Rt.*stsq.*Rtp.*cpsq-2.*ctcu.*Rsq.*Rt.*Rtp-ct.*R.*stsq.*Rp.*Rtsq-ct.*Rsq.*Rp.*spsq.*Rtt.*stsq-R.*st.*Rpsq.*spsq.*Rtp-ct.*Rpcu.*spsq.*Rtt.*stsq-ct.*Rsq.*cpsq.*Rt.*Rtp+2.*ct.*R.*cpsq.*Rp.*Rtsq+ct.*Rtsq.*stsq.*cpsq.*Rp.*Rpp-ct.*R.*stsq.*Rpcu-R.*st.*Rtp.*Rpsq+2.*ct.*R.*Rpcu.*spsq.*stsq+ct.*Rsq.*stsq.*Rp.*Rpp+ct.*Rsq.*stsq.*Rt.*Rtp-ct.*Rsq.*cpsq.*Rp.*Rpp-spsq.*Rsq.*ct.*Rt.*Rtp-spsq.*Rcu.*ctcu.*Rp+spsq.*Rcu.*ct.*Rp+2.*spsq.*R.*ct.*Rpcu-spsq.*R.*ctcu.*Rpcu-stsq.*Rp.*cpsq.*Rtt.*ct.*Rsq+2.*stsq.*Rp.*cpsq.*Rcu.*ct+ctcu.*Rtsq.*Rp.*Rpp-3.*spsq.*R.*ctcu.*Rp.*Rtsq-spsq.*Rsq.*ct.^5.*Rt.*Rtp-spsq.*Rsq.*ct.*Rp.*Rpp+spsq.*Rsq.*ctcu.*Rp.*Rpp-ctqad.*Rcu.*st.*Rtp-Rcu.*st.*cpsq.*Rtp-3.*stsq.*Rp.*cpsq.*Rtsq.*ctcu.*R-st.*Rp.*cpsq.*R.*ctqad.*Rt.*Rtt+2.*stsq.*Rpcu.*cpsq.*R.*ct+4.*Rt.*ctsq.*st.*Rpcu-3.*Rtsq.*ct.*Rp.*R-3.*Rtsq.*ct.^5.*R.*Rp+6.*Rtsq.*ctcu.*Rp.*R-3.*Rtcu.*ctqad.*st.*Rp+3.*Rtcu.*ctsq.*Rp.*st-4.*spsq.*Rt.*st.*Rpcu.*ctsq-7.*spsq.*Rtcu.*st.*Rp.*ctsq-spsq.*Rtsq.*st.*R.*ctqad.*Rtp+3.*spsq.*Rtsq.*st.*R.*ctsq.*Rtp+ct.^5.*Rsq.*Rt.*Rtp-spsq.*Rt.*stsq.*Rtp.*ct.*Rsq+spsq.*Rt.*stsq.*Rtp.*ctcu.*Rsq+spsq.*Rtsq.*stsq.*ct.*Rp.*Rpp+4.*st.*spsq.*Rp.*Rtcu+4.*Rtcu.*st.*cpsq.*Rp-stsq.*Rp.*cpsq.*Rcu.*ctcu+st.*Rp.*cpsq.*R.*Rt.*ctsq.*Rtt-4.*st.*Rpcu.*cpsq.*Rt.*ctsq+4.*stsq.*Rp.*cpsq.*Rtsq.*ct.*R-spsq.*Rt.*st.*R.*Rp.*Rpp-2.*spsq.*Rtsq.*st.*R.*Rtp+4.*spsq.*Rtsq.*stsq.*ct.*R.*Rp+2.*spsq.*Rt.*st.*R.*ctsq.*Rp.*Rpp-3.*spsq.*Rtsq.*stsq.*ctcu.*R.*Rp+spsq.*Rcu.*ctqad.*Rtp.*st+2.*spsq.*Rsq.*ctcu.*Rt.*Rtp+spsq.*R.*ct.^5.*Rp.*Rtsq+2.*spsq.*R.*ct.*Rp.*Rtsq+3.*spsq.*Rtcu.*st.*Rp.*ctqad+4.*spsq.*Rt.*st.*Rpcu-ct.*Rp.*Rcu-ctcu.*Rpcu.*Rtt+2.*ctcu.*Rpcu.*R+ct.*Rpcu.*Rtt+2.*ctcu.*Rp.*Rcu-ct.^5.*Rp.*Rcu-ct.*Rpcu.*R-3.*st.*Rsq.*spsq.*Rp.*ctsq.*Rt-3.*st.*Rp.*cpsq.*Rsq.*ctsq.*Rt+3.*st.*Rp.*cpsq.*Rtcu.*ctqad-st.*Rpsq.*cpsq.*R.*Rtp-7.*st.*Rp.*cpsq.*Rtcu.*ctsq-Rsq.*st.*Rp.*ctsq.*Rt-2.*ctcu.*Rp.*Rtt.*Rsq-ctsq.*Rp.*R.*st.*Rt.*Rtt+ctqad.*Rp.*R.*st.*Rt.*Rtt+ct.^5.*Rp.*Rtt.*Rsq+ct.*Rp.*Rtt.*Rsq+4.*st.*Rpcu.*cpsq.*Rt-st.*Rcu.*spsq.*Rtp+Rsq.*st.*Rp.*Rt+R.*st.*Rt.*Rp.*Rpp+3.*Rsq.*Rt.*st.*cpsq.*Rp+3.*Rsq.*st.*spsq.*Rp.*Rt-stsq.*Rpcu.*cpsq.*Rtt.*ct+stsq.*Rp.*cpsq.*Rtt.*ctcu.*Rsq+2.*stsq.*Rcu.*spsq.*Rp.*ct+st.*R.*spsq.*Rp.*Rt.*ctsq.*Rtt)./(blk_I).^(3./2);
%     N(:) =(Rp.*st.*cp-R.*st.*sp).*(3.*Rp.*cp.*Rtsq.*ctsq-Rtp.*ctcu.*st.*cp.*Rsq+2.*Rpsq.*sp.*R+R.*sp.*Rtsq+Rpcu.*ctsq.*cp-Rsq.*cp.*Rp+Rcu.*sp.*ctqad+Rtp.*ct.*st.*cp.*Rpsq-2.*Rp.*cp.*Rtsq+Rtcu.*ctcu.*st.*sp-Rpp.*sp.*Rtsq+R.*sp.*ctqad.*Rtsq+Rcu.*sp+Rp.*sp.*Rt.*Rtp-Rpp.*sp.*Rsq+R.*cp.*ctqad.*Rt.*Rtp-Rtcu.*ct.*st.*sp-Rp.*ctqad.*cp.*Rtsq-Rp.*sp.*Rt.*ctsq.*Rtp-2.*R.*cp.*Rt.*ctsq.*Rtp-2.*Rcu.*sp.*ctsq+Rt.*ctcu.*st.*sp.*Rsq+R.*cp.*Rp.*Rpp-Rt.*ct.*st.*cp.*Rp.*Rpp-Rt.*ct.*st.*cp.*R.*Rp+R.*cp.*Rt.*Rtp+Rpp.*sp.*Rsq.*ctsq-2.*Rpsq.*sp.*R.*ctsq-2.*R.*sp.*Rtsq.*ctsq+Rtp.*ct.*st.*cp.*Rsq-Rt.*ct.*st.*sp.*Rpsq+Rpp.*sp.*Rtsq.*ctsq+Rt.*ctcu.*st.*cp.*R.*Rp+Rsq.*cp.*ctsq.*Rp-2.*Rpcu.*cp-R.*cp.*ctsq.*Rp.*Rpp-Rt.*ct.*st.*sp.*Rsq)./(blk_I).^(3./2)+(Rp.*st.*sp+R.*st.*cp).*(-Rtp.*ct.*st.*sp.*Rsq+Rtp.*ctcu.*st.*sp.*Rsq-2.*Rcu.*cp.*ctsq-Rt.*ct.*st.*cp.*Rpsq+Rp.*cp.*Rt.*Rtp-Rtcu.*ct.*st.*cp+2.*R.*sp.*Rt.*ctsq.*Rtp+2.*Rpsq.*cp.*R+Rp.*sp.*ctqad.*Rtsq+2.*Rp.*sp.*Rtsq-Rtp.*ct.*st.*sp.*Rpsq-Rsq.*sp.*ctsq.*Rp-2.*Rpsq.*cp.*R.*ctsq+R.*cp.*Rtsq-Rpp.*cp.*Rtsq+Rsq.*sp.*Rp+Rt.*ct.*st.*sp.*Rp.*Rpp+2.*Rpcu.*sp-Rpcu.*sp.*ctsq-3.*Rp.*sp.*Rtsq.*ctsq-R.*sp.*ctqad.*Rt.*Rtp+Rtcu.*ctcu.*st.*cp+Rpp.*cp.*Rsq.*ctsq+Rcu.*cp-Rt.*ct.*st.*cp.*Rsq-Rpp.*cp.*Rsq+Rcu.*cp.*ctqad-R.*sp.*Rp.*Rpp-2.*R.*cp.*Rtsq.*ctsq-R.*sp.*Rt.*Rtp+R.*cp.*ctqad.*Rtsq-Rp.*cp.*Rt.*ctsq.*Rtp+Rt.*ct.*st.*sp.*R.*Rp+R.*sp.*ctsq.*Rp.*Rpp-Rt.*ctcu.*st.*sp.*R.*Rp+Rt.*ctcu.*st.*cp.*Rsq+Rpp.*cp.*Rtsq.*ctsq)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(blk_I).^(1./2)-ct.*Rp.*(st.*Rpcu.*ct+Rtp.*ctqad.*Rsq-Rtp.*ctsq.*Rpsq-2.*Rtp.*Rsq.*ctsq-Rp.*R.*Rt-Rt.*Rp.*Rpp-Rt.*ctqad.*R.*Rp+Rt.*ctsq.*Rp.*Rpp+2.*Rp.*R.*ctsq.*Rt+Rtp.*Rpsq+Rtp.*Rsq-st.*Rp.*ctcu.*Rtsq+Rp.*Rtsq.*ct.*st+R.*st.*ctcu.*Rt.*Rtp-R.*st.*ct.*Rt.*Rtp-R.*st.*ct.*Rp.*Rpp)./(blk_I).^(3./2);
%     H(:) = (E.*N + G.*LS - 2.*F.*M)./(2*(E.*G - F.^2));  % local mean curvature
%     %%%%%%%% Calculate the energy
%     wb = abs((sum(sum(wp.*wt.*(2.*H).^2.*S))))/16/pi;
%     wb_vec(ix) = wb;
% end
% Bx = sum(wb_vec);
% wb_vec = [];
% for ix = 1:length(yclks)
%     x = zeros(size(yclks));x(ix) = yclks(ix);
%     %%%%%%%% Calculate the surface/shape properties
%     R(:)   =  r_o + Y_LK_g*x;Rp  = Y_LK_phi*x;Rt  = Y_LK_theta*x;Rpp = Y_PP*x;Rtt = Y_TT*x;Rtp = Y_TP*x;
%     %%% Calculate some expressions only once;
%     Rsq(:) = R.^2; Rcu(:) = R.^3;Rtsq(:) = Rt.^2; Rtcu(:) = Rt.^3;Rpsq(:) = Rp.^2; Rpcu(:) = Rp.^3;
%     blk_I(:) = (Rpsq+Rsq+Rtsq-Rtsq.*ctsq-Rsq.*ctsq);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     S(:) = R.*(blk_I).^(1./2);
%     E(:) = (Rtsq.*stsq+2.*Rt.*ct.*R.*st+Rsq.*ctsq).*spsq+Rtsq.*stsq.*cpsq+2.*Rt.*st.*cpsq.*R.*ct+Rsq.*ctsq.*cpsq-2.*Rt.*ct.*R.*st+Rsq.*stsq+Rtsq.*ctsq;
%     F(:) = Rp.*(stsq.*cpsq.*Rt+st.*cpsq.*R.*ct+stsq.*spsq.*Rt+st.*spsq.*R.*ct+ctsq.*Rt-ct.*R.*st);
%     G(:) = stsq.*Rpsq.*cpsq+stsq.*Rsq.*spsq+stsq.*Rpsq.*spsq+stsq.*Rsq.*cpsq+ctsq.*Rpsq;
%     LS(:) = (Rt.*st.*cp+R.*cp.*ct).*(Rtp.*sp.*Rsq-Rtt.*ct.*st.*cp.*Rpsq+Rtt.*ctcu.*st.*cp.*Rsq+Rsq.*cp.*Rt+2.*Rtcu.*ctqad.*cp+2.*Rtcu.*cp-Rtt.*ct.*st.*cp.*Rsq+Rtp.*sp.*Rtsq+Rcu.*cp.*ct.*st-Rcu.*cp.*ctcu.*st+2.*R.*cp.*Rt.*ctsq.*Rtt+Rsq.*cp.*ctqad.*Rt-Rp.*sp.*R.*Rt-2.*Rsq.*cp.*ctsq.*Rt+Rt.*ct.*st.*cp.*Rp.*Rtp+R.*cp.*ctsq.*Rp.*Rtp-Rp.*sp.*Rsq.*ct.*st-3.*Rt.*ctsq.*cp.*Rpsq-R.*cp.*Rt.*Rtt-2.*Rtsq.*ctcu.*st.*cp.*R-R.*cp.*ctqad.*Rt.*Rtt-R.*cp.*Rp.*Rtp-Rtp.*sp.*Rsq.*ctsq-Rtp.*sp.*Rtsq.*ctsq-Rp.*sp.*Rt.*Rtt+2.*R.*cp.*ct.*st.*Rpsq-Rp.*sp.*Rtsq.*ct.*st-4.*Rtcu.*cp.*ctsq+Rp.*sp.*Rt.*ctsq.*Rtt+2.*Rtsq.*ct.*st.*cp.*R+2.*Rt.*cp.*Rpsq+Rp.*sp.*R.*ctsq.*Rt)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(Rpsq+Rsq+Rtsq-Rtsq.*ctsq-Rsq.*ctsq).^(1./2)-(Rt.*st.*sp+R.*ct.*sp).*(-4.*Rtcu.*ctsq.*sp+2.*sp.*Rtcu.*ctqad+2.*sp.*Rt.*Rpsq+2.*R.*sp.*Rt.*ctsq.*Rtt-cp.*Rtp.*Rsq-cp.*Rtp.*Rtsq-Rtt.*ct.*st.*sp.*Rpsq+Rsq.*sp.*Rt+2.*Rtsq.*ct.*st.*sp.*R+cp.*Rtp.*Rtsq.*ctsq+Rtt.*ctcu.*st.*sp.*Rsq+Rt.*ct.*st.*sp.*Rp.*Rtp-Rcu.*sp.*ctcu.*st+Rp.*sp.*R.*ctsq.*Rtp-Rtt.*ct.*st.*sp.*Rsq+cp.*Rtp.*Rsq.*ctsq+Rp.*cp.*Rtsq.*ct.*st+Rp.*cp.*Rt.*Rtt-3.*sp.*Rt.*ctsq.*Rpsq-Rp.*sp.*R.*Rtp+2.*R.*ct.*sp.*st.*Rpsq-R.*sp.*Rt.*Rtt-2.*Rsq.*sp.*ctsq.*Rt+Rsq.*sp.*ctqad.*Rt+2.*sp.*Rtcu-Rp.*cp.*Rt.*ctsq.*Rtt+Rp.*cp.*R.*Rt+Rcu.*sp.*ct.*st-R.*cp.*ctsq.*Rp.*Rt-2.*Rtsq.*ctcu.*st.*sp.*R+Rp.*cp.*Rsq.*ct.*st-R.*sp.*ctqad.*Rt.*Rtt)./(blk_I).^(3./2)-(Rt.*ct-R.*st).*(Rcu-Rt.*Rsq.*ct.*st-Rt.*ctsq.*Rp.*Rtp+Rt.*ctcu.*Rsq.*st-3.*Rt.*ct.*st.*Rpsq+2.*Rtt.*Rsq.*ctsq+R.*st.*ct.*Rt.*Rtt-R.*st.*ctcu.*Rt.*Rtt+2.*Rtcu.*ctcu.*st+2.*Rtsq.*ctqad.*R+Rtt.*ctsq.*Rpsq-Rtt.*ctqad.*Rsq-2.*R.*ctsq.*Rpsq+2.*R.*Rtsq+R.*st.*ct.*Rp.*Rtp-Rtt.*Rpsq-Rtt.*Rsq-2.*Rcu.*ctsq+Rcu.*ctqad+R.*Rpsq+Rt.*Rp.*Rtp-2.*Rtcu.*ct.*st-4.*R.*ctsq.*Rtsq)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(blk_I).^(1./2);
%     M(:) = 1./2.*(-Rcu.*st.*Rtp-ctcu.*Rcu.*cpsq.*Rp-ctqad.*R.*st.*spsq.*Rp.*Rt.*Rtt+ctqad.*Rcu.*st.*cpsq.*Rtp+ct.^5.*R.*cpsq.*Rp.*Rtsq-ct.^5.*Rsq.*cpsq.*Rt.*Rtp+ctqad.*R.*Rtsq.*st.*Rtp-ctqad.*R.*Rtsq.*st.*cpsq.*Rtp-ctcu.*R.*cpsq.*Rpcu+ct.*Rsq.*Rt.*Rtp+ctcu.*Rsq.*Rt.*stsq.*Rtp.*cpsq-2.*R.*Rtsq.*st.*cpsq.*Rtp+ctcu.*Rsq.*Rp.*spsq.*Rtt.*stsq+2.*ctsq.*Rcu.*st.*Rtp+ctcu.*R.*stsq.*Rp.*Rtsq-ctcu.*Rcu.*spsq.*Rp.*stsq+2.*ctcu.*Rsq.*cpsq.*Rt.*Rtp-3.*ctcu.*R.*cpsq.*Rp.*Rtsq-ctcu.*Rsq.*stsq.*Rt.*Rtp+ctcu.*Rsq.*cpsq.*Rp.*Rpp-ctsq.*R.*Rtsq.*st.*Rtp+ct.*Rcu.*cpsq.*Rp-ct.*Rtsq.*Rp.*Rpp+3.*ctsq.*R.*Rtsq.*st.*cpsq.*Rtp+2.*ct.*R.*cpsq.*Rpcu+2.*ctsq.*R.*Rt.*st.*cpsq.*Rp.*Rpp-2.*ctsq.*R.*st.*Rt.*Rp.*Rpp-R.*Rt.*st.*cpsq.*Rp.*Rpp-ct.*Rsq.*Rt.*stsq.*Rtp.*cpsq-2.*ctcu.*Rsq.*Rt.*Rtp-ct.*R.*stsq.*Rp.*Rtsq-ct.*Rsq.*Rp.*spsq.*Rtt.*stsq-R.*st.*Rpsq.*spsq.*Rtp-ct.*Rpcu.*spsq.*Rtt.*stsq-ct.*Rsq.*cpsq.*Rt.*Rtp+2.*ct.*R.*cpsq.*Rp.*Rtsq+ct.*Rtsq.*stsq.*cpsq.*Rp.*Rpp-ct.*R.*stsq.*Rpcu-R.*st.*Rtp.*Rpsq+2.*ct.*R.*Rpcu.*spsq.*stsq+ct.*Rsq.*stsq.*Rp.*Rpp+ct.*Rsq.*stsq.*Rt.*Rtp-ct.*Rsq.*cpsq.*Rp.*Rpp-spsq.*Rsq.*ct.*Rt.*Rtp-spsq.*Rcu.*ctcu.*Rp+spsq.*Rcu.*ct.*Rp+2.*spsq.*R.*ct.*Rpcu-spsq.*R.*ctcu.*Rpcu-stsq.*Rp.*cpsq.*Rtt.*ct.*Rsq+2.*stsq.*Rp.*cpsq.*Rcu.*ct+ctcu.*Rtsq.*Rp.*Rpp-3.*spsq.*R.*ctcu.*Rp.*Rtsq-spsq.*Rsq.*ct.^5.*Rt.*Rtp-spsq.*Rsq.*ct.*Rp.*Rpp+spsq.*Rsq.*ctcu.*Rp.*Rpp-ctqad.*Rcu.*st.*Rtp-Rcu.*st.*cpsq.*Rtp-3.*stsq.*Rp.*cpsq.*Rtsq.*ctcu.*R-st.*Rp.*cpsq.*R.*ctqad.*Rt.*Rtt+2.*stsq.*Rpcu.*cpsq.*R.*ct+4.*Rt.*ctsq.*st.*Rpcu-3.*Rtsq.*ct.*Rp.*R-3.*Rtsq.*ct.^5.*R.*Rp+6.*Rtsq.*ctcu.*Rp.*R-3.*Rtcu.*ctqad.*st.*Rp+3.*Rtcu.*ctsq.*Rp.*st-4.*spsq.*Rt.*st.*Rpcu.*ctsq-7.*spsq.*Rtcu.*st.*Rp.*ctsq-spsq.*Rtsq.*st.*R.*ctqad.*Rtp+3.*spsq.*Rtsq.*st.*R.*ctsq.*Rtp+ct.^5.*Rsq.*Rt.*Rtp-spsq.*Rt.*stsq.*Rtp.*ct.*Rsq+spsq.*Rt.*stsq.*Rtp.*ctcu.*Rsq+spsq.*Rtsq.*stsq.*ct.*Rp.*Rpp+4.*st.*spsq.*Rp.*Rtcu+4.*Rtcu.*st.*cpsq.*Rp-stsq.*Rp.*cpsq.*Rcu.*ctcu+st.*Rp.*cpsq.*R.*Rt.*ctsq.*Rtt-4.*st.*Rpcu.*cpsq.*Rt.*ctsq+4.*stsq.*Rp.*cpsq.*Rtsq.*ct.*R-spsq.*Rt.*st.*R.*Rp.*Rpp-2.*spsq.*Rtsq.*st.*R.*Rtp+4.*spsq.*Rtsq.*stsq.*ct.*R.*Rp+2.*spsq.*Rt.*st.*R.*ctsq.*Rp.*Rpp-3.*spsq.*Rtsq.*stsq.*ctcu.*R.*Rp+spsq.*Rcu.*ctqad.*Rtp.*st+2.*spsq.*Rsq.*ctcu.*Rt.*Rtp+spsq.*R.*ct.^5.*Rp.*Rtsq+2.*spsq.*R.*ct.*Rp.*Rtsq+3.*spsq.*Rtcu.*st.*Rp.*ctqad+4.*spsq.*Rt.*st.*Rpcu-ct.*Rp.*Rcu-ctcu.*Rpcu.*Rtt+2.*ctcu.*Rpcu.*R+ct.*Rpcu.*Rtt+2.*ctcu.*Rp.*Rcu-ct.^5.*Rp.*Rcu-ct.*Rpcu.*R-3.*st.*Rsq.*spsq.*Rp.*ctsq.*Rt-3.*st.*Rp.*cpsq.*Rsq.*ctsq.*Rt+3.*st.*Rp.*cpsq.*Rtcu.*ctqad-st.*Rpsq.*cpsq.*R.*Rtp-7.*st.*Rp.*cpsq.*Rtcu.*ctsq-Rsq.*st.*Rp.*ctsq.*Rt-2.*ctcu.*Rp.*Rtt.*Rsq-ctsq.*Rp.*R.*st.*Rt.*Rtt+ctqad.*Rp.*R.*st.*Rt.*Rtt+ct.^5.*Rp.*Rtt.*Rsq+ct.*Rp.*Rtt.*Rsq+4.*st.*Rpcu.*cpsq.*Rt-st.*Rcu.*spsq.*Rtp+Rsq.*st.*Rp.*Rt+R.*st.*Rt.*Rp.*Rpp+3.*Rsq.*Rt.*st.*cpsq.*Rp+3.*Rsq.*st.*spsq.*Rp.*Rt-stsq.*Rpcu.*cpsq.*Rtt.*ct+stsq.*Rp.*cpsq.*Rtt.*ctcu.*Rsq+2.*stsq.*Rcu.*spsq.*Rp.*ct+st.*R.*spsq.*Rp.*Rt.*ctsq.*Rtt)./(blk_I).^(3./2);
%     N(:) =(Rp.*st.*cp-R.*st.*sp).*(3.*Rp.*cp.*Rtsq.*ctsq-Rtp.*ctcu.*st.*cp.*Rsq+2.*Rpsq.*sp.*R+R.*sp.*Rtsq+Rpcu.*ctsq.*cp-Rsq.*cp.*Rp+Rcu.*sp.*ctqad+Rtp.*ct.*st.*cp.*Rpsq-2.*Rp.*cp.*Rtsq+Rtcu.*ctcu.*st.*sp-Rpp.*sp.*Rtsq+R.*sp.*ctqad.*Rtsq+Rcu.*sp+Rp.*sp.*Rt.*Rtp-Rpp.*sp.*Rsq+R.*cp.*ctqad.*Rt.*Rtp-Rtcu.*ct.*st.*sp-Rp.*ctqad.*cp.*Rtsq-Rp.*sp.*Rt.*ctsq.*Rtp-2.*R.*cp.*Rt.*ctsq.*Rtp-2.*Rcu.*sp.*ctsq+Rt.*ctcu.*st.*sp.*Rsq+R.*cp.*Rp.*Rpp-Rt.*ct.*st.*cp.*Rp.*Rpp-Rt.*ct.*st.*cp.*R.*Rp+R.*cp.*Rt.*Rtp+Rpp.*sp.*Rsq.*ctsq-2.*Rpsq.*sp.*R.*ctsq-2.*R.*sp.*Rtsq.*ctsq+Rtp.*ct.*st.*cp.*Rsq-Rt.*ct.*st.*sp.*Rpsq+Rpp.*sp.*Rtsq.*ctsq+Rt.*ctcu.*st.*cp.*R.*Rp+Rsq.*cp.*ctsq.*Rp-2.*Rpcu.*cp-R.*cp.*ctsq.*Rp.*Rpp-Rt.*ct.*st.*sp.*Rsq)./(blk_I).^(3./2)+(Rp.*st.*sp+R.*st.*cp).*(-Rtp.*ct.*st.*sp.*Rsq+Rtp.*ctcu.*st.*sp.*Rsq-2.*Rcu.*cp.*ctsq-Rt.*ct.*st.*cp.*Rpsq+Rp.*cp.*Rt.*Rtp-Rtcu.*ct.*st.*cp+2.*R.*sp.*Rt.*ctsq.*Rtp+2.*Rpsq.*cp.*R+Rp.*sp.*ctqad.*Rtsq+2.*Rp.*sp.*Rtsq-Rtp.*ct.*st.*sp.*Rpsq-Rsq.*sp.*ctsq.*Rp-2.*Rpsq.*cp.*R.*ctsq+R.*cp.*Rtsq-Rpp.*cp.*Rtsq+Rsq.*sp.*Rp+Rt.*ct.*st.*sp.*Rp.*Rpp+2.*Rpcu.*sp-Rpcu.*sp.*ctsq-3.*Rp.*sp.*Rtsq.*ctsq-R.*sp.*ctqad.*Rt.*Rtp+Rtcu.*ctcu.*st.*cp+Rpp.*cp.*Rsq.*ctsq+Rcu.*cp-Rt.*ct.*st.*cp.*Rsq-Rpp.*cp.*Rsq+Rcu.*cp.*ctqad-R.*sp.*Rp.*Rpp-2.*R.*cp.*Rtsq.*ctsq-R.*sp.*Rt.*Rtp+R.*cp.*ctqad.*Rtsq-Rp.*cp.*Rt.*ctsq.*Rtp+Rt.*ct.*st.*sp.*R.*Rp+R.*sp.*ctsq.*Rp.*Rpp-Rt.*ctcu.*st.*sp.*R.*Rp+Rt.*ctcu.*st.*cp.*Rsq+Rpp.*cp.*Rtsq.*ctsq)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(blk_I).^(1./2)-ct.*Rp.*(st.*Rpcu.*ct+Rtp.*ctqad.*Rsq-Rtp.*ctsq.*Rpsq-2.*Rtp.*Rsq.*ctsq-Rp.*R.*Rt-Rt.*Rp.*Rpp-Rt.*ctqad.*R.*Rp+Rt.*ctsq.*Rp.*Rpp+2.*Rp.*R.*ctsq.*Rt+Rtp.*Rpsq+Rtp.*Rsq-st.*Rp.*ctcu.*Rtsq+Rp.*Rtsq.*ct.*st+R.*st.*ctcu.*Rt.*Rtp-R.*st.*ct.*Rt.*Rtp-R.*st.*ct.*Rp.*Rpp)./(blk_I).^(3./2);
%     H(:) = (E.*N + G.*LS - 2.*F.*M)./(2*(E.*G - F.^2));  % local mean curvature
%     %%%%%%%% Calculate the energy
%     wb = abs((sum(sum(wp.*wt.*(2.*H).^2.*S))))/16/pi;
%     wb_vec(ix) = wb;
% end
% By = sum(wb_vec);
% wb_vec = [];
% for ix = 1:length(zclks)
%     x = zeros(size(zclks));x(ix) = zclks(ix);
%     %%%%%%%% Calculate the surface/shape properties
%     R(:)   =  r_o + Y_LK_g*x;Rp  = Y_LK_phi*x;Rt  = Y_LK_theta*x;Rpp = Y_PP*x;Rtt = Y_TT*x;Rtp = Y_TP*x;
%     %%% Calculate some expressions only once;
%     Rsq(:) = R.^2; Rcu(:) = R.^3;Rtsq(:) = Rt.^2; Rtcu(:) = Rt.^3;Rpsq(:) = Rp.^2; Rpcu(:) = Rp.^3;
%     blk_I(:) = (Rpsq+Rsq+Rtsq-Rtsq.*ctsq-Rsq.*ctsq);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     S(:) = R.*(blk_I).^(1./2);
%     E(:) = (Rtsq.*stsq+2.*Rt.*ct.*R.*st+Rsq.*ctsq).*spsq+Rtsq.*stsq.*cpsq+2.*Rt.*st.*cpsq.*R.*ct+Rsq.*ctsq.*cpsq-2.*Rt.*ct.*R.*st+Rsq.*stsq+Rtsq.*ctsq;
%     F(:) = Rp.*(stsq.*cpsq.*Rt+st.*cpsq.*R.*ct+stsq.*spsq.*Rt+st.*spsq.*R.*ct+ctsq.*Rt-ct.*R.*st);
%     G(:) = stsq.*Rpsq.*cpsq+stsq.*Rsq.*spsq+stsq.*Rpsq.*spsq+stsq.*Rsq.*cpsq+ctsq.*Rpsq;
%     LS(:) = (Rt.*st.*cp+R.*cp.*ct).*(Rtp.*sp.*Rsq-Rtt.*ct.*st.*cp.*Rpsq+Rtt.*ctcu.*st.*cp.*Rsq+Rsq.*cp.*Rt+2.*Rtcu.*ctqad.*cp+2.*Rtcu.*cp-Rtt.*ct.*st.*cp.*Rsq+Rtp.*sp.*Rtsq+Rcu.*cp.*ct.*st-Rcu.*cp.*ctcu.*st+2.*R.*cp.*Rt.*ctsq.*Rtt+Rsq.*cp.*ctqad.*Rt-Rp.*sp.*R.*Rt-2.*Rsq.*cp.*ctsq.*Rt+Rt.*ct.*st.*cp.*Rp.*Rtp+R.*cp.*ctsq.*Rp.*Rtp-Rp.*sp.*Rsq.*ct.*st-3.*Rt.*ctsq.*cp.*Rpsq-R.*cp.*Rt.*Rtt-2.*Rtsq.*ctcu.*st.*cp.*R-R.*cp.*ctqad.*Rt.*Rtt-R.*cp.*Rp.*Rtp-Rtp.*sp.*Rsq.*ctsq-Rtp.*sp.*Rtsq.*ctsq-Rp.*sp.*Rt.*Rtt+2.*R.*cp.*ct.*st.*Rpsq-Rp.*sp.*Rtsq.*ct.*st-4.*Rtcu.*cp.*ctsq+Rp.*sp.*Rt.*ctsq.*Rtt+2.*Rtsq.*ct.*st.*cp.*R+2.*Rt.*cp.*Rpsq+Rp.*sp.*R.*ctsq.*Rt)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(Rpsq+Rsq+Rtsq-Rtsq.*ctsq-Rsq.*ctsq).^(1./2)-(Rt.*st.*sp+R.*ct.*sp).*(-4.*Rtcu.*ctsq.*sp+2.*sp.*Rtcu.*ctqad+2.*sp.*Rt.*Rpsq+2.*R.*sp.*Rt.*ctsq.*Rtt-cp.*Rtp.*Rsq-cp.*Rtp.*Rtsq-Rtt.*ct.*st.*sp.*Rpsq+Rsq.*sp.*Rt+2.*Rtsq.*ct.*st.*sp.*R+cp.*Rtp.*Rtsq.*ctsq+Rtt.*ctcu.*st.*sp.*Rsq+Rt.*ct.*st.*sp.*Rp.*Rtp-Rcu.*sp.*ctcu.*st+Rp.*sp.*R.*ctsq.*Rtp-Rtt.*ct.*st.*sp.*Rsq+cp.*Rtp.*Rsq.*ctsq+Rp.*cp.*Rtsq.*ct.*st+Rp.*cp.*Rt.*Rtt-3.*sp.*Rt.*ctsq.*Rpsq-Rp.*sp.*R.*Rtp+2.*R.*ct.*sp.*st.*Rpsq-R.*sp.*Rt.*Rtt-2.*Rsq.*sp.*ctsq.*Rt+Rsq.*sp.*ctqad.*Rt+2.*sp.*Rtcu-Rp.*cp.*Rt.*ctsq.*Rtt+Rp.*cp.*R.*Rt+Rcu.*sp.*ct.*st-R.*cp.*ctsq.*Rp.*Rt-2.*Rtsq.*ctcu.*st.*sp.*R+Rp.*cp.*Rsq.*ct.*st-R.*sp.*ctqad.*Rt.*Rtt)./(blk_I).^(3./2)-(Rt.*ct-R.*st).*(Rcu-Rt.*Rsq.*ct.*st-Rt.*ctsq.*Rp.*Rtp+Rt.*ctcu.*Rsq.*st-3.*Rt.*ct.*st.*Rpsq+2.*Rtt.*Rsq.*ctsq+R.*st.*ct.*Rt.*Rtt-R.*st.*ctcu.*Rt.*Rtt+2.*Rtcu.*ctcu.*st+2.*Rtsq.*ctqad.*R+Rtt.*ctsq.*Rpsq-Rtt.*ctqad.*Rsq-2.*R.*ctsq.*Rpsq+2.*R.*Rtsq+R.*st.*ct.*Rp.*Rtp-Rtt.*Rpsq-Rtt.*Rsq-2.*Rcu.*ctsq+Rcu.*ctqad+R.*Rpsq+Rt.*Rp.*Rtp-2.*Rtcu.*ct.*st-4.*R.*ctsq.*Rtsq)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(blk_I).^(1./2);
%     M(:) = 1./2.*(-Rcu.*st.*Rtp-ctcu.*Rcu.*cpsq.*Rp-ctqad.*R.*st.*spsq.*Rp.*Rt.*Rtt+ctqad.*Rcu.*st.*cpsq.*Rtp+ct.^5.*R.*cpsq.*Rp.*Rtsq-ct.^5.*Rsq.*cpsq.*Rt.*Rtp+ctqad.*R.*Rtsq.*st.*Rtp-ctqad.*R.*Rtsq.*st.*cpsq.*Rtp-ctcu.*R.*cpsq.*Rpcu+ct.*Rsq.*Rt.*Rtp+ctcu.*Rsq.*Rt.*stsq.*Rtp.*cpsq-2.*R.*Rtsq.*st.*cpsq.*Rtp+ctcu.*Rsq.*Rp.*spsq.*Rtt.*stsq+2.*ctsq.*Rcu.*st.*Rtp+ctcu.*R.*stsq.*Rp.*Rtsq-ctcu.*Rcu.*spsq.*Rp.*stsq+2.*ctcu.*Rsq.*cpsq.*Rt.*Rtp-3.*ctcu.*R.*cpsq.*Rp.*Rtsq-ctcu.*Rsq.*stsq.*Rt.*Rtp+ctcu.*Rsq.*cpsq.*Rp.*Rpp-ctsq.*R.*Rtsq.*st.*Rtp+ct.*Rcu.*cpsq.*Rp-ct.*Rtsq.*Rp.*Rpp+3.*ctsq.*R.*Rtsq.*st.*cpsq.*Rtp+2.*ct.*R.*cpsq.*Rpcu+2.*ctsq.*R.*Rt.*st.*cpsq.*Rp.*Rpp-2.*ctsq.*R.*st.*Rt.*Rp.*Rpp-R.*Rt.*st.*cpsq.*Rp.*Rpp-ct.*Rsq.*Rt.*stsq.*Rtp.*cpsq-2.*ctcu.*Rsq.*Rt.*Rtp-ct.*R.*stsq.*Rp.*Rtsq-ct.*Rsq.*Rp.*spsq.*Rtt.*stsq-R.*st.*Rpsq.*spsq.*Rtp-ct.*Rpcu.*spsq.*Rtt.*stsq-ct.*Rsq.*cpsq.*Rt.*Rtp+2.*ct.*R.*cpsq.*Rp.*Rtsq+ct.*Rtsq.*stsq.*cpsq.*Rp.*Rpp-ct.*R.*stsq.*Rpcu-R.*st.*Rtp.*Rpsq+2.*ct.*R.*Rpcu.*spsq.*stsq+ct.*Rsq.*stsq.*Rp.*Rpp+ct.*Rsq.*stsq.*Rt.*Rtp-ct.*Rsq.*cpsq.*Rp.*Rpp-spsq.*Rsq.*ct.*Rt.*Rtp-spsq.*Rcu.*ctcu.*Rp+spsq.*Rcu.*ct.*Rp+2.*spsq.*R.*ct.*Rpcu-spsq.*R.*ctcu.*Rpcu-stsq.*Rp.*cpsq.*Rtt.*ct.*Rsq+2.*stsq.*Rp.*cpsq.*Rcu.*ct+ctcu.*Rtsq.*Rp.*Rpp-3.*spsq.*R.*ctcu.*Rp.*Rtsq-spsq.*Rsq.*ct.^5.*Rt.*Rtp-spsq.*Rsq.*ct.*Rp.*Rpp+spsq.*Rsq.*ctcu.*Rp.*Rpp-ctqad.*Rcu.*st.*Rtp-Rcu.*st.*cpsq.*Rtp-3.*stsq.*Rp.*cpsq.*Rtsq.*ctcu.*R-st.*Rp.*cpsq.*R.*ctqad.*Rt.*Rtt+2.*stsq.*Rpcu.*cpsq.*R.*ct+4.*Rt.*ctsq.*st.*Rpcu-3.*Rtsq.*ct.*Rp.*R-3.*Rtsq.*ct.^5.*R.*Rp+6.*Rtsq.*ctcu.*Rp.*R-3.*Rtcu.*ctqad.*st.*Rp+3.*Rtcu.*ctsq.*Rp.*st-4.*spsq.*Rt.*st.*Rpcu.*ctsq-7.*spsq.*Rtcu.*st.*Rp.*ctsq-spsq.*Rtsq.*st.*R.*ctqad.*Rtp+3.*spsq.*Rtsq.*st.*R.*ctsq.*Rtp+ct.^5.*Rsq.*Rt.*Rtp-spsq.*Rt.*stsq.*Rtp.*ct.*Rsq+spsq.*Rt.*stsq.*Rtp.*ctcu.*Rsq+spsq.*Rtsq.*stsq.*ct.*Rp.*Rpp+4.*st.*spsq.*Rp.*Rtcu+4.*Rtcu.*st.*cpsq.*Rp-stsq.*Rp.*cpsq.*Rcu.*ctcu+st.*Rp.*cpsq.*R.*Rt.*ctsq.*Rtt-4.*st.*Rpcu.*cpsq.*Rt.*ctsq+4.*stsq.*Rp.*cpsq.*Rtsq.*ct.*R-spsq.*Rt.*st.*R.*Rp.*Rpp-2.*spsq.*Rtsq.*st.*R.*Rtp+4.*spsq.*Rtsq.*stsq.*ct.*R.*Rp+2.*spsq.*Rt.*st.*R.*ctsq.*Rp.*Rpp-3.*spsq.*Rtsq.*stsq.*ctcu.*R.*Rp+spsq.*Rcu.*ctqad.*Rtp.*st+2.*spsq.*Rsq.*ctcu.*Rt.*Rtp+spsq.*R.*ct.^5.*Rp.*Rtsq+2.*spsq.*R.*ct.*Rp.*Rtsq+3.*spsq.*Rtcu.*st.*Rp.*ctqad+4.*spsq.*Rt.*st.*Rpcu-ct.*Rp.*Rcu-ctcu.*Rpcu.*Rtt+2.*ctcu.*Rpcu.*R+ct.*Rpcu.*Rtt+2.*ctcu.*Rp.*Rcu-ct.^5.*Rp.*Rcu-ct.*Rpcu.*R-3.*st.*Rsq.*spsq.*Rp.*ctsq.*Rt-3.*st.*Rp.*cpsq.*Rsq.*ctsq.*Rt+3.*st.*Rp.*cpsq.*Rtcu.*ctqad-st.*Rpsq.*cpsq.*R.*Rtp-7.*st.*Rp.*cpsq.*Rtcu.*ctsq-Rsq.*st.*Rp.*ctsq.*Rt-2.*ctcu.*Rp.*Rtt.*Rsq-ctsq.*Rp.*R.*st.*Rt.*Rtt+ctqad.*Rp.*R.*st.*Rt.*Rtt+ct.^5.*Rp.*Rtt.*Rsq+ct.*Rp.*Rtt.*Rsq+4.*st.*Rpcu.*cpsq.*Rt-st.*Rcu.*spsq.*Rtp+Rsq.*st.*Rp.*Rt+R.*st.*Rt.*Rp.*Rpp+3.*Rsq.*Rt.*st.*cpsq.*Rp+3.*Rsq.*st.*spsq.*Rp.*Rt-stsq.*Rpcu.*cpsq.*Rtt.*ct+stsq.*Rp.*cpsq.*Rtt.*ctcu.*Rsq+2.*stsq.*Rcu.*spsq.*Rp.*ct+st.*R.*spsq.*Rp.*Rt.*ctsq.*Rtt)./(blk_I).^(3./2);
%     N(:) =(Rp.*st.*cp-R.*st.*sp).*(3.*Rp.*cp.*Rtsq.*ctsq-Rtp.*ctcu.*st.*cp.*Rsq+2.*Rpsq.*sp.*R+R.*sp.*Rtsq+Rpcu.*ctsq.*cp-Rsq.*cp.*Rp+Rcu.*sp.*ctqad+Rtp.*ct.*st.*cp.*Rpsq-2.*Rp.*cp.*Rtsq+Rtcu.*ctcu.*st.*sp-Rpp.*sp.*Rtsq+R.*sp.*ctqad.*Rtsq+Rcu.*sp+Rp.*sp.*Rt.*Rtp-Rpp.*sp.*Rsq+R.*cp.*ctqad.*Rt.*Rtp-Rtcu.*ct.*st.*sp-Rp.*ctqad.*cp.*Rtsq-Rp.*sp.*Rt.*ctsq.*Rtp-2.*R.*cp.*Rt.*ctsq.*Rtp-2.*Rcu.*sp.*ctsq+Rt.*ctcu.*st.*sp.*Rsq+R.*cp.*Rp.*Rpp-Rt.*ct.*st.*cp.*Rp.*Rpp-Rt.*ct.*st.*cp.*R.*Rp+R.*cp.*Rt.*Rtp+Rpp.*sp.*Rsq.*ctsq-2.*Rpsq.*sp.*R.*ctsq-2.*R.*sp.*Rtsq.*ctsq+Rtp.*ct.*st.*cp.*Rsq-Rt.*ct.*st.*sp.*Rpsq+Rpp.*sp.*Rtsq.*ctsq+Rt.*ctcu.*st.*cp.*R.*Rp+Rsq.*cp.*ctsq.*Rp-2.*Rpcu.*cp-R.*cp.*ctsq.*Rp.*Rpp-Rt.*ct.*st.*sp.*Rsq)./(blk_I).^(3./2)+(Rp.*st.*sp+R.*st.*cp).*(-Rtp.*ct.*st.*sp.*Rsq+Rtp.*ctcu.*st.*sp.*Rsq-2.*Rcu.*cp.*ctsq-Rt.*ct.*st.*cp.*Rpsq+Rp.*cp.*Rt.*Rtp-Rtcu.*ct.*st.*cp+2.*R.*sp.*Rt.*ctsq.*Rtp+2.*Rpsq.*cp.*R+Rp.*sp.*ctqad.*Rtsq+2.*Rp.*sp.*Rtsq-Rtp.*ct.*st.*sp.*Rpsq-Rsq.*sp.*ctsq.*Rp-2.*Rpsq.*cp.*R.*ctsq+R.*cp.*Rtsq-Rpp.*cp.*Rtsq+Rsq.*sp.*Rp+Rt.*ct.*st.*sp.*Rp.*Rpp+2.*Rpcu.*sp-Rpcu.*sp.*ctsq-3.*Rp.*sp.*Rtsq.*ctsq-R.*sp.*ctqad.*Rt.*Rtp+Rtcu.*ctcu.*st.*cp+Rpp.*cp.*Rsq.*ctsq+Rcu.*cp-Rt.*ct.*st.*cp.*Rsq-Rpp.*cp.*Rsq+Rcu.*cp.*ctqad-R.*sp.*Rp.*Rpp-2.*R.*cp.*Rtsq.*ctsq-R.*sp.*Rt.*Rtp+R.*cp.*ctqad.*Rtsq-Rp.*cp.*Rt.*ctsq.*Rtp+Rt.*ct.*st.*sp.*R.*Rp+R.*sp.*ctsq.*Rp.*Rpp-Rt.*ctcu.*st.*sp.*R.*Rp+Rt.*ctcu.*st.*cp.*Rsq+Rpp.*cp.*Rtsq.*ctsq)./(Rtsq.*ctsq+Rsq.*ctsq-Rsq-Rpsq-Rtsq)./(blk_I).^(1./2)-ct.*Rp.*(st.*Rpcu.*ct+Rtp.*ctqad.*Rsq-Rtp.*ctsq.*Rpsq-2.*Rtp.*Rsq.*ctsq-Rp.*R.*Rt-Rt.*Rp.*Rpp-Rt.*ctqad.*R.*Rp+Rt.*ctsq.*Rp.*Rpp+2.*Rp.*R.*ctsq.*Rt+Rtp.*Rpsq+Rtp.*Rsq-st.*Rp.*ctcu.*Rtsq+Rp.*Rtsq.*ct.*st+R.*st.*ctcu.*Rt.*Rtp-R.*st.*ct.*Rt.*Rtp-R.*st.*ct.*Rp.*Rpp)./(blk_I).^(3./2);
%     H(:) = (E.*N + G.*LS - 2.*F.*M)./(2*(E.*G - F.^2));  % local mean curvature
%     %%%%%%%% Calculate the energy
%     wb = abs((sum(sum(wp.*wt.*(2.*H).^2.*S))))/16/pi;
%     wb_vec(ix) = wb;
% end
% Bz = sum(wb_vec);




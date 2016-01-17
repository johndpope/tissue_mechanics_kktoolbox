function X_odd = select_odd_L_forward(X_o)
% takes as input the usual X_o vector and outputs a vector that has the
% even L values removed (including L = 0)

X_odd = X_o;
% % % [xc yc zc] = get_xyz_clks(X_o);
% % % L_max = get_L_max(X_o);
% % % last = 0;
% % % xco = [];
% % % yco = [];
% % % zco = [];
% % % for ix = 0:L_max
% % %     start = last+1; finish = (ix + 1).^2;
% % %     vec = start:finish;
% % %     last = finish;
% % %     %vec
% % %     if rem(ix,2)~=0,        % check for odd ix
% % %         xco = [xco xc(vec)'];
% % %         yco = [yco yc(vec)'];
% % %         zco = [zco zc(vec)'];
% % %     end
% % % end
% % % X_odd = [xco(:)' yco(:)' zco(:)'];

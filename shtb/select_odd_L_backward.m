function X_o = select_odd_L_backward(X_odd, L_max)
% takes as input a vector that has the
% even L values removed (including L = 0), and outputs the usual X_o vector
% length, with the even L values set to 0

X_o = X_odd;

% % % nc = round(length(X_odd)/3);
% % % xco = X_odd(1:nc);
% % % yco = X_odd(nc+1:2*nc);
% % % zco = X_odd(2*nc+1:3*nc);
% % % 
% % % X_o = zeros(3*(L_max + 1)^2,1);
% % % [xc yc zc] = get_xyz_clks(X_o);
% % % 
% % % last = 0;
% % % vec = [];
% % % pos = 1;
% % % for ix = 0:L_max
% % %     start = last+1; finish = (ix + 1).^2;
% % %     vec_ix = start:finish;
% % %     last = finish;
% % %     if rem(ix,2)~=0 && ix~=0,     % i.e. if ix is odd
% % %         xc(vec_ix) = xco(pos:pos+length(vec_ix)-1);
% % %         yc(vec_ix) = yco(pos:pos+length(vec_ix)-1);
% % %         zc(vec_ix) = zco(pos:pos+length(vec_ix)-1);
% % %         pos = pos + length(vec_ix);
% % %     end
% % % end
% % % 
% % % X_o = [xc(:)' yc(:)' zc(:)'];

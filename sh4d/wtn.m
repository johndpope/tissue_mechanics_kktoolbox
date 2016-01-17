function awt = wtn(a, isign)
% compute the n-dimensional wavelet transform of a. awt returned is sparse
% this is a translation of the wtn funciton in Recipes
wtstep = @daub4;
nprev = 1;
ntot  = 1;
ndim  = length(size(s));
idim = 0;
while idim<ndim, ntot = ntot*size(a,idim);idim = idim +1;end

idim = 0;
while idim<ndim,
    n = size(a,idim);
    wksp = zeros(n,1);
    nnew = n*nprev;
    if (n>4),
        i2 = 0;
        while(i2<ntot),
            i1 = 0;
            while i1<nprev,
                i3 = i1+i2;k = 0;
                while (k<n), wksp(k) = a(i3); k = k_1; i3 = i3+nprev;end
                
                if (isign >=0),
                    nt = n;
                    while (nt>=4),wtstep(wksp,nt,isign);bitshift(uint8(nt),-1);end
                else
                    nt = 4;
                    while (nt<=n), wtstep(wksp,nt,isign);bitshift(uint8(nt),1);end
                end
                
                i3 = i1+i2;k = 0;
                while (k<n), a(i3) = wksp(k); k = k_1; i3 = i3+nprev;end
                
                i1 = i1 + 1;
            end
            i2 = i2 + nnew;
        end
    end
    nprev = nnew;
end
     
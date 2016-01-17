function [r c all_Xdr cm m] = match_objects_SHP_com(T, Td, dangle)
%% match the objects defined by the pair of coordinate sets rotated with the specified
%% rotation (angle)

gdim = 30;

[r c all_Xdr cm m] = ...
    match_objects_rendered(...
    get_all_X_from_T(T,gdim),...
    get_all_X_from_T(Td,gdim),dangle);


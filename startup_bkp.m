clc;
if ~isdeployed
    default_home = 'Z:\Khaled\mcode\kktoolbox';
    secondary_home = 'C:\Khaled\m_scripts\kktoolbox';
    if exist(default_home)==7
        addpath(genpath(default_home));
        disp(['Using: ' default_home]);
    else
        addpath(genpath(secondary_home));
        disp(['Using: ' secondary_home]);
    end
    cd C:\Khaled\m_scripts
    path
end

function mfn=cellmodelname(mmod)

if (mmod==1)
    mfn='01_Aliev-Panfilov.in';
elseif (mmod==2)
    mfn='02_Rogers-McCulloch.in';
elseif (mmod==3)
    mfn='03_Beeler-Reuter.in';
else
    error('... membrane model not available')
end

mfn=strcat('./template_files/',mfn);

end
ThePath='/opt/matlab/netcdf/';

javaaddpath([ThePath 'netcdfAll-4.2.jar']);
javaaddpath([ThePath 'mexcdf/snctools/classes']);
addpath([ThePath 'mexcdf/mexnc']);
addpath([ThePath 'mexcdf/snctools']);

eFile='field.nc';

TheX=nc_varget(eFile, 'x');
TheY=nc_varget(eFile, 'y');

ListTrig=nc_varget(eFile, 'ele');

LTime=nc_varget(eFile, 'ocean_time');
PrintOpt1='-dpng';
PrintOpt2='-r200';
TheExtension='.png';

nbTime=size(LTime,1);

fem_grid_struct.e=ListTrig;
fem_grid_struct.x=TheX;
fem_grid_struct.y=TheY;

len=size(TheX, 1);

VarName='HS';

for iTime=1:nbTime
    disp(['iTime=' num2str(iTime) ' / ' num2str(nbTime)]);
    eVAR=nc_varget(eFile, VarName, [iTime-1, 0], [1, len]);
    disp(['min=' num2str(min(eVAR(:))) ' max=' num2str(max(eVAR(:)))]);
    maxVAR=max(eVAR(:));
    minVAR=min(eVAR(:));
    disp(['VAR=' VarName ' min=' num2str(minVAR) '  max=' num2str(maxVAR)]);
    %
    hold on;
    title([VarName ' iTime=' num2str(iTime)]);
    colormesh2d(fem_grid_struct,eVAR);
    caxis([minVAR maxVAR]);
    colorbar;
    shading interp;
    FileName=['PLOT_' VarName '_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
end;

TheRed  =[1 0 0];
TheGreen=[0 1 0];
TheBlue =[0 0 1];
TheBlack=[0,0,0];


ThePath='/opt/matlab/netcdf/';

javaaddpath([ThePath 'netcdfAll-4.2.jar']);
javaaddpath([ThePath 'mexcdf/snctools/classes']);
addpath([ThePath 'mexcdf/mexnc']);
addpath([ThePath 'mexcdf/snctools']);

eFileHist ='history_0001.nc';

TheX=nc_varget(eFileHist, 'x');
TheY=nc_varget(eFileHist, 'y');

ListTrig=nc_varget(eFileHist, 'ele');

LTime=nc_varget(eFileHist, 'ocean_time');
PrintOpt1='-dpng';
PrintOpt2='-r200';
TheExtension='.png';

nbTime=size(LTime,1);

fem_grid_struct.e=ListTrig;
fem_grid_struct.x=TheX;
fem_grid_struct.y=TheY;

len=size(TheX, 1);

minX=min(TheX);
maxX=max(TheX);
minY=min(TheY);
maxY=max(TheY);

nbSub=100;
ListX=zeros(nbSub,1);
IX=zeros(nbSub,1);
for i=1:nbSub
    eX=minX + ((i-1)/(nbSub-1)) * (maxX - minX);
    eY=(minY + maxY) /2;
    [index,distance]=nearxy(TheX,TheY,eX,eY);
%    disp(['i=' num2str(i) ' eX=' num2str(eX) ...
%          ' index=' num2str(index)]);
    IX(i,1)=index;
end;


%
% Now the analytical computation
%

TheOpt='volume conservation';
%TheOpt='zero infinity';
DoDynBathy=0;
ChoiceDep=1;
N=400;



if (ChoiceDep == 1)
  aDep=0.45;
  TheLength=9;
  Ntot=TheLength*N;
  deltaX=1/N;
  h=TheLength/2;
  nTot=N*h;
  ListDep=aDep*ones(Ntot, 1);
  for i=1:nTot
    ListDep(i,1)=aDep*(i/nTot);
  end;
elseif (ChoiceDep == 2)
  aDep=0.45;
  TheLength=4.5;
  Ntot=TheLength*N;
  deltaX=1/N;
  h=TheLength;
  nTot=N*h;
  ListDep=aDep*ones(Ntot, 1);
  for i=1:nTot
    ListDep(i,1)=aDep*(i/nTot);
  end;
else
  disp('Please put what you have in mind');
  error('Please correct');
end;
ListX=zeros(Ntot,1);
for i=1:Ntot
  ListX(i,1)=i/N;
end;
ListXmid=(ListX(2:Ntot,1) + ListX(1:Ntot-1,1))/2;
ListZeta=zeros(Ntot,1);
ListK=zeros(Ntot,1);
ListHmono=zeros(Ntot,1);
T=1.5;
omega=2*pi/T;
ListOmega=omega*ones(Ntot,1);
gAccel=9.81;
%alphaBreaking=0.415;
%alphaBreakingHS=0.78;
%alphaBreakingHS=0.82*sqrt(2);
alphaBreakingHS=0.81;
%alphaBreakingHS=0.71;
%alphaBreakingHS=2*alphaBreaking;

alphaBreaking = alphaBreakingHS/sqrt(2);
TheDifferential=zeros(Ntot-1,1);
%aHS=0.256;
%aHS=2*0.09*sqrt(2);
aHS=0.1801;
%aHS=0.09*sqrt(2);
MaxErr=10^(-14);
iter=0;
nbIterMax=400;


while(1)
  if (DoDynBathy == 1)
    DynBathy=ListDep+ListZeta;
  else
    DynBathy=ListDep;
  end;
  iter=iter+1;
  ListKexact=GetKvectorExact(ListOmega, DynBathy);
  kD=ListKexact.*DynBathy;
  Lwave=2*pi./ListKexact;
  nNumber=(1/2) + kD./sinh(2*kD);
  cPhase=omega./ListKexact;
  cGroup=cPhase.*nNumber;
  %
  Cst=aHS*aHS*cGroup(Ntot,1);
  ListHnobreak=sqrt(Cst./cGroup);
  MaxBreaking=alphaBreakingHS*DynBathy;
  TheDiff=abs(MaxBreaking - ListHnobreak);
  minDiff=min(TheDiff);
  K=find(TheDiff == minDiff);
  aV=K(1,1);
  BreakingPoint=ListX(aV);
  ListHsignificant=(min([ListHnobreak'; MaxBreaking']))';
  ListHmono=ListHsignificant/sqrt(2);
  %
  eK=2*nNumber-1/2;
  ListSxx=(1/16)*(ListHsignificant.*ListHsignificant.*eK);
  for i=2:Ntot
    TheDifferential(i-1,1)=(ListSxx(i,1) - ListSxx(i-1,1))/deltaX;
  end;
  ListZetaNew=zeros(Ntot,1);
  ListZetaNew(Ntot,1)=0;
  for i=2:Ntot
    j=Ntot+1-i;
    Dmid=(DynBathy(j,1)+DynBathy(j+1,1))/2;
    TheDiff=-(1/Dmid)*TheDifferential(j,1);
    ListZetaNew(j,1)=ListZetaNew(j+1,1)-TheDiff*deltaX;
  end;
  if (UTIL_IsStringEqual(TheOpt, 'volume conservation') == 1)
    TheSum=sum(ListZetaNew);
    ListZetaNew=ListZetaNew - TheSum/Ntot;
  end;
  DZetaDx=zeros(Ntot-1,1);
  for i=2:Ntot
    dMid=(DynBathy(i-1,1)+DynBathy(i,1))/2;
    ZetaDx=(ListZetaNew(i,1) - ListZetaNew(i-1,1))/deltaX;
    DZetaDx(i-1,1)=dMid*ZetaDx;
  end;
  absDiff=sum(abs(ListZetaNew(:) - ListZeta(:)));
  disp(['iter=' num2str(iter) ' absDiff=' num2str(absDiff)]);
  ListZeta=ListZetaNew;
  minZeta=min(ListZeta);
  aIncrease=max(ListHsignificant) - aHS;
  disp(['minZETA=' num2str(minZeta) ...
	' breakPoint=' num2str(BreakingPoint) ...
	' aIncrease=' num2str(aIncrease)]);
  if (absDiff < MaxErr)
    break;
  end;
  if (iter > nbIterMax)
    break;
  end;
end;


delta_ListX=max(ListX) - min(ListX);
delta_TheX =max(TheX)  - min(TheX);

eScal_X=delta_ListX / delta_TheX;


for iTime=1:nbTime
    disp(['iTime=' num2str(iTime) ' / ' num2str(nbTime)]);
    eVAR_hs =nc_varget(eFileHist, 'HS' , [iTime-1, 0], [1, len]);
    eVAR_stp=nc_varget(eFileHist, 'ZETA_SETUP', [iTime-1, 0], [1, len]);
    eVAR_dpt=nc_varget(eFileHist, 'DEP', [iTime-1, 0], [1, len]);
    eVAR_t01=nc_varget(eFileHist, 'TM01', [iTime-1, 0], [1, len]);
    eVAR_t02=nc_varget(eFileHist, 'TM02', [iTime-1, 0], [1, len]);

    eVAR_sxx=nc_varget(eFileHist, 'RSXX', [iTime-1, 0], [1, len]);
    eVAR_sxy=nc_varget(eFileHist, 'RSXY', [iTime-1, 0], [1, len]);
    eVAR_syy=nc_varget(eFileHist, 'RSYY', [iTime-1, 0], [1, len]);
    %
    hold on;
    TheTitle=['HS longitudinal plot'];
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_hs(IX), 'Color', TheBlue);
    line(ListX, ListHsignificant, 'Color', TheRed);
    FileName=['LINE_HS_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
    %
    hold on;
    TheTitle=['STP longitudinal plot'];
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_stp(IX), 'Color', TheBlue);
    line(ListX, ListZeta, 'Color', TheRed);
    FileName=['LINE_stp_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
    %
    hold on;
    TheTitle='DPT longitudinal plot';
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_dpt(IX), 'Color', TheBlue);
    axis([0 9 0 0.50]);
    line(ListX, ListDep, 'Color', TheRed);
    FileName=['LINE_dpt_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
    %
    hold on;
    TheTitle='T01 longitudinal plot';
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_t01(IX), 'Color', TheBlue);
    FileName=['LINE_t01_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
    %
    hold on;
    TheTitle='T02 longitudinal plot';
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_t02(IX), 'Color', TheBlue);
    FileName=['LINE_t02_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;

    %
    hold on;
    TheTitle='SXX longitudinal plot';
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_sxx(IX), 'Color', TheBlue);
    FileName=['LINE_sxx_' StringNumber(iTime, 4) TheExtension];
    line(ListX, ListSxx, 'Color', TheRed);
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
    %
    hold on;
    TheTitle='SXY longitudinal plot';
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_sxy(IX), 'Color', TheBlue);
    FileName=['LINE_sxy_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
    %
    hold on;
    TheTitle='SYY longitudinal plot';
    title(TheTitle);
    line(eScal_X*TheX(IX), eVAR_syy(IX), 'Color', TheBlue);
    FileName=['LINE_syy_' StringNumber(iTime, 4) TheExtension];
    disp(['FileName=' FileName]);
    print(PrintOpt1, PrintOpt2, FileName);
    hold off;
    clf;
end;
eVAR_hsB=eVAR_hs(IX);
eMaxHS=max(eVAR_hsB);
K=find(eVAR_hsB == eMaxHS);
idx=K(1,1);

eDep_break=eVAR_dpt(IX(idx));
eHS_break =eVAR_hs(IX(idx));
alphaBreaking=eHS_break / eDep_break;
disp(['alphaBreaking=' num2str(alphaBreaking)]);

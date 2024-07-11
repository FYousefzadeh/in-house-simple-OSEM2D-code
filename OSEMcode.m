clc; clear all;

param.nx = 64;
param.ny = 64;
param.nz = 64;

param.sx = 64; % mm
param.sy = 64; % mm
param.sz = 64; % mm

% Detector setting, according to Varian Trilogy OBI (real size)
param.sv = 64;	% mm
param.su = 64;	% mm

%The real detector panel pixel density (number of pixels)
param.nv = 64;
param.nu = 64;

% X-ray source and detector setting
%param.DSD = 900000;    %  Distance source to detector
%param.DSO = 900000;	%  X-ray source to object axis distance

% angle setting
Dir = -1;   % gantry rotating direction
param.deg = -45:180/32:134;
%param.deg = StartAngle:AngularStep:(NumberOfFramesInRotation-1)*AngularStep+StartAngle;
param.deg = param.deg*Dir;
param.nProj = length(param.deg);

% filter='ram-lak','cosine', 'hamming', 'hann'
param.filter='none'; % high pass for sintetic images

param.dx = param.sx/param.nx;
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

param.off_u = 0; param.off_v = 0; % detector rotation shift (real size)

% % % For fast CPU calculation % % %
param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx;
param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz;

param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v;


%% load Image and Normalizing Factors
One=ones(param.nx,param.ny,param.nz);
for i=1:param.nz
    PrNorm(:,:,i)=radon(One(:,:,i),param.deg);
end
PrNorm=permute(PrNorm,[1,3,2]);%/max(proj(:));
S1=ceil((size(PrNorm,1)-param.nx)/2);
PrNorm=PrNorm(S1+1:S1+param.nx,:,:);


proj=single(permute(squeeze(dicomread('00005 10sec.dcm')),[2,1,3]));

Iter_no=8;
Subiter=2;
tic
Image=OSEM(proj,Iter_no,Subiter,param);
toc
%%


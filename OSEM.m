function [ImageArchive]=OSEM(proj,Iteration,Subiter,param)

img = ones(param.nx, param.ny, param.nz);

clear unitproj

ImageArchive=zeros([size(img),Iteration]);

%%
Rot=zeros(Subiter,param.nProj/Subiter);
RotAngle=param.deg;
for j=1:Subiter
    for n=1:param.nProj/Subiter
        Rot(j,n)=RotAngle((Subiter*(n-1)+j));
    end
end


%% MLEM
for iter = 1:Iteration
    for Sub=1:Subiter
        
        
        
        %     tic
        %     fprintf('------------------------------\n');
        %     fprintf('       Please Wait       \n');
        %     fprintf('MLEM Algorithm under Runing\n');
        %     fprintf('  Press Ctrl+c to Abort \n');
        %     fprintf('______________________________\n');
        %     pause(0.01);
        
        for i=1:size(proj,2)
            Project(:,:,i)=radon(img(:,:,i),Rot(Sub,:));
        end
        Project=permute(Project,[1,3,2]);
        S1=ceil((size(Project,1)-param.nu)/2);
        Project=Project(S1+1:S1+param.nu,:,:);
        
        
        
        
        One=ones(param.nx,param.ny,param.nz);
        for i=1:param.nz
            PrNorm(:,:,i)=radon(One(:,:,i),Rot(Sub,:));
        end
        PrNorm=permute(PrNorm,[1,3,2]);%/max(proj(:));
        S1=ceil((size(PrNorm,1)-param.nx)/2);
        PrNorm=PrNorm(S1+1:S1+param.nx,:,:);
        
        Project=Project./PrNorm;
        
        clear PrNorm
        
        
        
        
        proj_ratio = proj(:,:,Sub:Subiter:param.nProj)./Project;
        clear Project
        
        
        
        proj_ratio(isnan(proj_ratio)) = 0;
        proj_ratio(isinf(proj_ratio)) = 0;
        
        for j=1:size(img,3)
            img_ratio(:,:,j) = iradon(squeeze(proj_ratio(:,j,:)), Rot(Sub,:), 'pchip', param.filter, param.nx)/(pi/2);
        end
        
        
        unitproj=ones(param.nu, param.nv, param.nProj/Subiter);
        for i=1:size(proj,2)
            Norimg(:,:,i) = iradon(squeeze(unitproj(:,i,:)), Rot(Sub,:), 'pchip', param.filter ,param.nx)/(pi/2);
        end
        
        img_ratio(isnan(img_ratio)) = 0;
        img_ratio(isinf(img_ratio)) = 0;
        img = img.*img_ratio./Norimg;
        
        clear Norimg
        %     toc
        
        ImageArchive(:,:,:,iter)=img;
        
        %     figure(2); imagesc(max(img(:,:,45),0)); axis off; axis equal; colormap gray; colorbar;
        %     title(['Iteration - ',num2str(iter)]);
        %     pause(0.01);
    end
end
end
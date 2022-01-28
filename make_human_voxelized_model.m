run C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\startup.m
skipcurr=4;
skipvolt=4;
N=16;
finemodel=mk_common_model('j2T3',N);
boundarynodes=finemodel.fwd_model.boundary;
boundarycoordinates1=finemodel.fwd_model.nodes(boundarynodes(:,1),:);
boundarycoordinates2=finemodel.fwd_model.nodes(boundarynodes(:,2),:);
b1=1;
for b=1:length(boundarycoordinates1)
    boundarycoordinates(b1,:)=boundarycoordinates1(b,:);
    boundarycoordinates(b1+1,:)=boundarycoordinates2(b,:);
    b1=b1+2;
end

for el=1:N
electrode_nodes(el,:)=finemodel.fwd_model.electrode(el).nodes;
end
electrode_coords1=finemodel.fwd_model.nodes(electrode_nodes(:,1),:);
electrode_coords2=finemodel.fwd_model.nodes(electrode_nodes(:,2),:);
ell=1;
for el=1:N
    electrode_coords(ell,:)=electrode_coords1(el,:);
    electrode_coords(ell+1,:)=electrode_coords2(el,:);
    electrode_coordinates(el,:)=(electrode_coords(ell,:)+electrode_coords(ell+1,:))/2;
    ell=ell+2;
end
%%%%%make coarse mesh
coarseres=41;
meshcoarse=-150:300/coarseres:150;
[X1c,Y1c]=meshgrid(meshcoarse);
xlc=X1c(:);
ylc=Y1c(:);
inc = inpolygon(xlc,ylc,boundarycoordinates(:,1),boundarycoordinates(:,2));
xc=xlc(inc); yc=ylc(inc);

figure(1)
plot(xc,yc,'*')
hold on
plot(boundarycoordinates(:,1),boundarycoordinates(:,2),'r')
hold on
plot(electrode_coordinates(:,1),electrode_coordinates(:,2),'ko')

vtxc= [xc(:),yc(:)];
% vtx=[vtx; electrode_coordinates];
for el=1:el
    elec_nodes{el}=[electrode_coordinates(el,1) electrode_coordinates(el,2)];
end
coarse_fwd_mdl=mk_fmdl_from_nodes(vtxc, elec_nodes,1,'finetest1');
figure(2)
show_fem(coarse_fwd_mdl)
title('coarse model')

%%%%make fine mesh
fineres=6*41;
meshfine=-150:2*300/fineres:150;
[X1f,Y1f]=meshgrid(meshfine);
xlf=X1f(:);
ylf=Y1f(:);
infm = inpolygon(xlf,ylf,boundarycoordinates(:,1),boundarycoordinates(:,2));
xf=xlf(infm); yf=ylf(infm);
%%%%%introduce z axis (3D)
height=max(xc);
zax=0:max(diff(xf)):height;
zc=0:max(diff(xc)):height;


x3dc=zeros(length(xc),length(zc));
y3dc=zeros(length(yc),length(zc));
z3dc=repmat(zc,[length(xc) 1]);
for ilc=1:length(xc)
    x3dc(ilc,:)=xc(ilc);
    y3dc(ilc,:)=yc(ilc);
end
xc3d=x3dc(:); yc3d=y3dc(:); zc3d=z3dc(:); 

x3df=zeros(length(xf),length(zax));
y3df=zeros(length(yf),length(zax));
z3df=repmat(zax,[length(xf) 1]);
for ilf=1:length(xf)
    x3df(ilf,:)=xf(ilf);
    y3df(ilf,:)=yf(ilf);
end
%%%%%%%%square placement ?_?_?...
electrodeind=1;
for el=1:2:N
    elec_nodes3d{electrodeind}=[electrode_coordinates(el,1) electrode_coordinates(el,2) height*0.35];
    electrodeind=electrodeind+1;
    elec_nodes3d{electrodeind}=[electrode_coordinates(el,1) electrode_coordinates(el,2) height*0.65];
    electrodeind=electrodeind+1;
    elec_nodes3d{electrodeind}=[electrode_coordinates(el+1,1) electrode_coordinates(el+1,2) height*0.35];
    electrodeind=electrodeind+1;
    elec_nodes3d{electrodeind}=[electrode_coordinates(el+1,1) electrode_coordinates(el+1,2) height*0.65];
    electrodeind=electrodeind+1;
end

figure
plot3(x3df,y3df,z3df,'b*')
axis tight

vtxf=[x3df(:),y3df(:), z3df(:)];
fine_fwd_mdl=mk_fmdl_from_nodes(vtxf, elec_nodes3d,1,'finetest1');
figure
show_fem(fine_fwd_mdl,[0,1])
title('fine model')

fprintf('passed refinement\n')

%%%%LEAD MODEL
[stim_lead, meas_sel_lead] = mk_stim_patterns(2*N,1,[0 skipvolt+1],[0 skipvolt+1],{},0.001);
fineleadmodel=fine_fwd_mdl;
fineleadmodel.stimulation= stim_lead;
imgfinelead = mk_image(fineleadmodel,1);

imgfinelead.fwd_solve.get_all_meas = 1;

leadfwd_solution=fwd_solve(imgfinelead);

fprintf('passed fwd\n')

green_diffs_lead=leadfwd_solution.volt;


fine_fnodes=fine_fwd_mdl.nodes;
for el=1:N*2
elnodes(el)=fine_fwd_mdl.electrode(el).nodes;
end
fine_fnodes(elnodes,:)=[];
% if norm(fine_fnodes-vtxf)~=0
%     error('non-matching nodes!')
% end
green_diffs_lead2=green_diffs_lead;
green_diffs_lead2(elnodes,:)=[];

x3dnew=fine_fnodes(:,1); y3dnew=fine_fnodes(:,2); z3dnew=fine_fnodes(:,3);

figure
scatter3(x3dnew,y3dnew,z3dnew,25,green_diffs_lead2(:,2),'filled')
axis equal
view([-30 30])
colormap jet



%%%% INJECT MODEL

if skipcurr==skipvolt
    green_diffs_inject2=green_diffs_lead2;
else
[stim_inject, meas_sel_inject] = mk_stim_patterns(2*N,2,[0 skipcurr+1],[0 skipcurr+1],{},0.001);
fineinjectmodel=fine_fwd_mdl;
fineinjectmodel.stimulation= stim_inject;
imgfineinject = mk_image(fineinjectmodel,1);
imgfineinject.fwd_solve.get_all_meas = 1;

injectfwd_solution=fwd_solve(imgfineinject);
green_diffs_inject=injectfwd_solution.volt;
green_diffs_inject2=green_diffs_inject;
green_diffs_inject2(elnodes,:)=[];
end

LF=length(x3dnew);
%%%%%%%%%% compute gradient numerically
%%%%%direction x
xunique=unique(x3dnew);
yunique=unique(y3dnew);
zunique=unique(z3dnew);
DFdX=zeros(LF,N);
for el=1:N*2
    for ylevel=1:length(yunique)
        for zlevel=1:length(zunique)
            xis_of_yzlevel=find(y3dnew==yunique(ylevel)&z3dnew==zunique(zlevel));
            DFdX(xis_of_yzlevel,el)=[diff(green_diffs_inject2(xis_of_yzlevel,el))./...
                diff(x3dnew(xis_of_yzlevel)); 0];
        end
    end
end

%%%%%direction y
DFdY=zeros(LF,N);
for el=1:N*2
    for xlevel=1:length(xunique)
        for zlevel=1:length(zunique)
            yis_of_xzlevel=find(x3dnew==xunique(xlevel)&z3dnew==zunique(zlevel));
            DFdY(yis_of_xzlevel,el)=[diff(green_diffs_inject2(yis_of_xzlevel,el))./...
                diff(y3dnew(yis_of_xzlevel)); 0];
        end
    end
end

%%%%%direction z
DFdZ=zeros(LF,N);
for el=1:N*2
    for xlevel=1:length(xunique)
        for ylevel=1:length(yunique)
            zis_of_xylevel=find(x3dnew==xunique(xlevel)&y3dnew==yunique(ylevel));
            DFdZ(zis_of_xylevel,el)=[diff(green_diffs_inject2(zis_of_xylevel,el))./...
                diff(z3dnew(zis_of_xylevel)); 0];
        end
    end
end

figure
scatter3(x3dnew,y3dnew,z3dnew,25,DFdX(:,1),'filled')
axis equal 

figure
scatter3(x3dnew,y3dnew,z3dnew,25,DFdY(:,1),'filled')
axis equal 

figure
scatter3(x3dnew,y3dnew,z3dnew,25,DFdZ(:,1),'filled')
axis equal 

%%%%%%find coarse mesh pairs indices in fine mesh
%%%(to get differenceGreen and gradientVoltage samples from the fine model)

coarsefindices=zeros(length(xc3d),1);
pointsno=1;
eps=1e-05;
for point=1:length(xc3d)
    indice=find(xc3d(point)>x3dnew-eps&xc3d(point)<x3dnew+eps&...
        yc3d(point)>y3dnew-eps&yc3d(point)<y3dnew+eps&...
        zc3d(point)>z3dnew-eps&zc3d(point)<z3dnew+eps);
    coarsefindices(pointsno)=indice;
    pointsno=pointsno+1;
end

%%%%%% coarse model potentials (lead)
green_diffs_coarse=green_diffs_lead2(coarsefindices,:);
green_diffs_coarse=(green_diffs_coarse-mean(green_diffs_coarse,1))*1000;
figure
scatter3(xc3d,yc3d,zc3d,70,green_diffs_coarse(:,32),'filled')
axis equal

%%%%%% coarse model gradients (injection)
green_diffs_inject_coarsex=DFdX(coarsefindices,:);
figure
scatter3(xc3d,yc3d,zc3d,70,green_diffs_inject_coarsex(:,2),'filled')
axis equal

green_diffs_inject_coarsey=DFdY(coarsefindices,:);
figure
scatter3(xc3d,yc3d,zc3d,70,green_diffs_inject_coarsey(:,2),'filled')
axis equal

green_diffs_inject_coarsez=DFdZ(coarsefindices,:);
figure
scatter3(xc3d,yc3d,zc3d,70,green_diffs_inject_coarsez(:,2),'filled')
axis equal


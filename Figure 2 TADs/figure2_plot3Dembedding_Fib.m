%% Set view angles for different models
close all
view_points = [ 140.70   46.2; % TAD A 119-122-123
                 80.7    50.4; % TAD T1, 335-336-340
                 75.0    11.4; % TAD B 124-126-127
                186.3    72.6]; % TAD B 124-125-126

%% Construct edges
for i_plot = 1:4
clear edges edgeWt nodeIds

    % Define node IDs, edges/hyperedges, and their weights
    if i_plot == 1 % TAD T1, 334-337-338-340 hyperegde
        nodeIds=[334:341];
        edges{1}=[334 335];             edgeWt(1)=1;
        edges{end+1}=[335 336];         edgeWt(end+1)=1;
        edges{end+1}=[336 337];         edgeWt(end+1)=1;
        edges{end+1}=[338 339];         edgeWt(end+1)=1;
        edges{end+1}=[339 340];         edgeWt(end+1)=1;
        edges{end+1}=[340 341];         edgeWt(end+1)=1;
        edges{end+1}=[334 337 338 340]; edgeWt(end+1)=10000000;

    elseif i_plot == 2 % TAD T1, 335-336-340 hyperdege
        nodeIds=[334:341];
        edges{1}=[334 335];             edgeWt(1)=1;
        edges{end+1}=[335 336];         edgeWt(end+1)=1;
        edges{end+1}=[336 337];         edgeWt(end+1)=1;
        edges{end+1}=[338 339];         edgeWt(end+1)=1;
        edges{end+1}=[339 340];         edgeWt(end+1)=1;
        edges{end+1}=[340 341];         edgeWt(end+1)=1;
        edges{end+1}=[335 336 340];     edgeWt(end+1)=100000;

    elseif i_plot == 3 % TAD T2, 343-344-345
        nodeIds=[342:346];
        edges{1}=[342 343];             edgeWt(1)=1;
        edges{end+1}=[343 344];         edgeWt(end+1)=1;
        edges{end+1}=[344 345];         edgeWt(end+1)=1;
        edges{end+1}=[345 346];         edgeWt(end+1)=1;
        edges{end+1}=[343 344 345];     edgeWt(end+1)=100000;
        
    elseif i_plot == 4 % TAD T2, 343-344-346
        nodeIds=[342:346];
        edges{1}=[342 343];             edgeWt(1)=1;
        edges{end+1}=[343 344];         edgeWt(end+1)=1;
        edges{end+1}=[344 345];         edgeWt(end+1)=1;
        edges{end+1}=[345 346];         edgeWt(end+1)=1;
        edges{end+1}=[343 344 346];     edgeWt(end+1)=100000;

%     elseif i_plot == 5 % Inter TAD ==> NOT UPDATED FOR FIB (not in use)
%         nodeIds=[118:127];
%         edges{1}=[118 119];             edgeWt(1)=1;
%         edges{end+1}=[119 120];         edgeWt(end+1)=1;
%         edges{end+1}=[120 121];         edgeWt(end+1)=1;
%         edges{end+1}=[121 122];         edgeWt(end+1)=1;
%         edges{end+1}=[122 123];         edgeWt(end+1)=1;
%         edges{end+1}=[123 124];         edgeWt(end+1)=1;
%         edges{end+1}=[124 125];         edgeWt(end+1)=1;
%         edges{end+1}=[125 126];         edgeWt(end+1)=1;
%         edges{end+1}=[126 127];         edgeWt(end+1)=1;
%         edges{end+1}=[118 123 126];     edgeWt(end+1)=100000;
%         edges{end+1}=[121 122 124];     edgeWt(end+1)=100000;
    end

    %% Create Incidence, Adjacency, and Distance Matrices
    exponent=1/3;
    eps=10^-8;

    % Create incidence matrix
    H=zeros(length(nodeIds),length(edges));
    for i=1:length(edges)
       ed=edges{i}; 
       for j=1:length(ed)
            n1=find(nodeIds==ed(j));
            H(n1,i)=1;
       end
    end

    % Clique expansion for adjacency matrix
    A=H*sparse(1:length(edgeWt),1:length(edgeWt),edgeWt,length(edgeWt),length(edgeWt))*H';
    A=A-diag(diag(A));
    A(A==0)=eps;

    % Convert to distance matrix
    D=1./(A.^exponent);
    D=D-diag(diag(D));
    cord=mdscale(D,3);

    %% Plot 3-D curve containing points and hyperedges of interest
    % c = jet(length(x)+1);       % can be array of size Nx3 or Nx1
    h = fnplt(cscvn(cord'),'r',3);

    figure
    plot3(cord(:,1),cord(:,2),cord(:,3),'co','linewidth',5)
    hold on
    fnplt(cscvn(cord'),'r',3)
    shading interp
    text(cord(:,1),cord(:,2),cord(:,3),num2str(nodeIds'))
    axis off

    for i=1:length(edges)
        if length(edges{i})>=3
            ed=edges{i};
            for j=1:length(ed)
                n1=find(nodeIds==ed(j));
                for k=j+1:length(ed)
                    n2=find(nodeIds==ed(k));
                    plot3([cord(n1,1) cord(n2,1)],[cord(n1,2) cord(n2,2)],[cord(n1,3) cord(n2,3)],'k-.')
                end
            end
        end
    end
    set(gca,'View',view_points(i_plot,:))
end

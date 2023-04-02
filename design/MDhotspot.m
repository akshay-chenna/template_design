%% Peptide design from MD without MMPBSA-WSAS

clc 
clearvars
%% 1. Load relevant data
% Load the delta area of binding per residue for NSP12:

delarea=readmatrix("delsasares.txt","FileType","text","OutputType","double"); 
delarea(end-71+1:end,:)=[]; 
%% 
% Load the cluster size and member information:

clustrows=readmatrix("clustrows25bb.txt","FileType","text","OutputType","double"); 
clustrows(isnan(clustrows))=0;
%% 
% Load the templates for peptide design from gmx mdat for the central structure 
% of that cluster:

mdmat=flip(reshape(readmatrix("mdmatclust25.txt","FileType","text","OutputType","double"),890,890,[]));
mdmat(:,end,:)=[]; % CASE SPECIFIC ! SHOULD HAVE RECTIFIED IN gmx mdmat ITSELF. 
mdmat(end,:,:)=[];
%% 2. Find the hotspots from each cluster

for cl=1:1%size(clustrows,2)
    harea=abs(delarea(:,clustrows(:,cl)>0));
    harea=mean(harea,2); % harea corresponds to the hotspot area
    harea(harea<0.001)=0;
    harea=harea/sum(harea);
%% 3. Segments from generated hotspots

    residues=find(harea);
    residues(length(residues)+1,1)=max(residues)+4; % pseudo value to run the "for" loop easily
    c=1;
    r=1;
    seg=[];
    for i=1:length(residues)-1
        if (residues(i+1)-residues(i))<=3
            res(r,c)=residues(i);
            r=r+1;
        else     
            res(r,c)=residues(i);
            maxval(c)=max(res(:,c))-min(res(res(:,c)>0,c))+1;
            r=1;
            c=c+1;
        end
    end
    clear c r    
    
    seg=zeros(max(maxval),size(res,2));
    for i=1:size(res,2)
        range=(min(res(res(:,i)>0,i)):max(res(:,i)))';
        seg(:,i)=[range;zeros(max(maxval)-length(range),1)];
    end
    clear i
%% 4. Find the inter-residue distances
% 1. Classify all inter-residue distances
% Classify residues based on the interresidue distance:

    mdmatA=triu(mdmat(:,:,cl));
    [x{1,1}(:,1), x{1,1}(:,2)]=find(mdmatA<=0.7 & mdmatA>0);
    [x{2,1}(:,1), x{2,1}(:,2)]=find(mdmatA<=1 & mdmatA>0.7);
    [x{3,1}(:,1), x{3,1}(:,2)]=find(mdmatA<=1.3 & mdmatA>1);
    calpha=x;
    clear mdmatA x;
% 2. Remove those within the same segment or not at the interface

     for i=1:size(calpha,1)
         for j=1:length(calpha{i,1})
             a=ismember(calpha{i,1}(j,1),seg); 
             b=ismember(calpha{i,1}(j,2),seg);
             if a==0 || b==0 % if either are not part of the interface, remove the entire row
                calpha{i,1}(j,:)=0;
             end
             [~,c1]=find(seg==calpha{i,1}(j,1));
             [~,c2]=find(seg==calpha{i,1}(j,2));
             if c1==c2 % if they are part of the same segment, delete the whole row
                 calpha{i,1}(j,:)=0;
             end
         end
         calpha{i,1}=reshape(nonzeros(calpha{i,1}),[],2); % remove zeros
     end
     clear a b c1 c2 i j
%% 5. Design peptides
% 1.  Intersegment peptides

     for h=1:length(calpha)
         c=1; %just a flag
         for i=1:size(calpha{h,1},1) 
%% 
% Find the position of the residue in the segment along with the maximum number 
% of residues in that segment:

             [r1,c1]=find(seg==calpha{h,1}(i,1));
             [r2,c2]=find(seg==calpha{h,1}(i,2));
             l1=length(find(seg(:,c1)));
             l2=length(find(seg(:,c2)));
%% 
% Construct peptide sequences: Could have used *int2aa* function from bioinformatics 
% toolbox to directly write amino acids instead of numbers

            %Both moving upwards
            for l=1:r1 
                for m=1:r2
                    pep{h,1}(60,c)=10000;
                    for n=1:m
                        for o=1:l
                            pep{h,1}(60-o,c)=seg(r1-o+1,c1);
                        end
                        pep{h,1}(60+n,c)=seg(r2-n+1,c2);
                    end
                        c=c+1;
                end
            end
            
            %Both moving downwards
            for l=r1:l1 
                for m=r2:l2
                    pep{h,1}(60,c)=10000;
                    for n=1:m-r2+1
                        for o=1:l-r1+1
                            pep{h,1}(60-o,c)=seg(r1+o-1,c1);
                        end
                        pep{h,1}(60+n,c)=seg(r2+n-1,c2);
                    end
                        c=c+1;
                end
            end
           
            %1 moves upwards, 2 moves down
             for l=1:r1 
                for m=r2:l2
                    pep{h,1}(60,c)=10000;
                    for n=1:m-r2+1
                        for o=1:l
                            pep{h,1}(60-o,c)=seg(r1-o+1,c1);
                        end
                        pep{h,1}(60+n,c)=seg(r2+n-1,c2);
                    end
                        c=c+1;
                end
             end
             
            %1 moves down and 2 moves up
             for l=r1:l1 
                for m=1:r2
                    pep{h,1}(60,c)=10000;
                    for n=1:m
                        for o=1:l-r1+1
                            pep{h,1}(60-o,c)=seg(r1+o-1,c1);
                        end
                        pep{h,1}(60+n,c)=seg(r2-n+1,c2);
                    end
                        c=c+1;
                end
             end
         end
     end
% 2. Peptides from the same segment

     c=1;
     for q=3:length(seg)                        %Set the minimum size of the peptide to three.
        for r=1:size(seg,2)
            n=max(find(seg(:,r)));
            for i=1:n+1-q
                x=0;
                a=1;
                while x<q
                    pep{h+1,1}(a,c)=seg(i+x,r);
                    x=x+1;                      % add upto q residues in that peptide chain
                    a=a+1;
                end
                c=c+1;                          %store peptides in different columns
            end
        end
     end
     clear a c c1 c2 h i j l l1 l2 m n o q r r1 r2 x 
     
% 3. Sort the peptides
% Add Alanine bridges
% Adding single alanine:

     if ~isempty(pep{2,1})
         for i=1:size(pep{2,1},2)
             a=length(find(pep{2,1}(1:59,i)));
             b=length(find(pep{2,1}(61:end,i)));
            if a<b || a==b
                temp{1,1}(:,i)=[pep{2,1}(1:59,i); 10000; 20000; pep{2,1}(61:end,i)];
            else
                temp{1,1}(:,i)=[pep{2,1}(1:59,i); 20000; 10000; pep{2,1}(61:end,i)];
            end
         end
        pep{2,1}=temp{1,1};
        clear temp
     end
%% 
% Add two alanines:

    if ~isempty(pep{3,1})
        for i=1:size(pep{3,1},2)
             a=length(find(pep{3,1}(1:59,i)));
             b=length(find(pep{3,1}(61:end,i)));
            if a<b || a==b
                temp{1,1}(:,i)=[pep{3,1}(1:59,i); 10000; 20000; 20000; pep{3,1}(61:end,i)];
            else
                temp{1,1}(:,i)=[pep{3,1}(1:59,i); 20000; 20000; 10000; pep{3,1}(61:end,i)];
            end
         end
        pep{3,1}=temp{1,1};
        clear a b i temp
        end
% Categorize the peptides by their length

    for x=1:size(pep,1)
        sz(x)=size(pep{x,1},1);
    end
    for i=1:size(pep,1)
        for j=1:max(sz)
            seq{j,i}=zeros(j,1);
        end
    end
    for i=1:size(pep,1)
        for j=1:size(pep{i,1},2)
            a=length(find(pep{i,1}(:,j)));
            seq{a,i}=[seq{a,i} pep{i,1}(find(pep{i,1}(:,j)),j)];
        end
    end
    clear a i j x pep sz
    for i=1:size(seq,1)
            a=[seq{i,1} seq{i,2} seq{i,3} seq{i,4}]';
            pep{i,1}=unique(a,'rows','stable')';
            pep{i,1}(:,1)=[];
    end
    pep=pep(~cellfun('isempty',pep)); % to remove empty cells in a cell array
    clear seq i a
    
%% 6. Peptide analysis
% 1. Find the hotpsot area contibution from each peptide and sort them descendingly 
% Finding the areas

    weight=length(find(clustrows(:,cl)))/max(max(clustrows)); %weight of the cluster
    for i=1:size(pep,1)
        for j=1:size(pep{i,1},2)
            k=pep{i,1}(pep{i,1}(:,j)<10000,j);
            area{i,1}(1,j)=sum(harea(k,1)).*weight;
        end
    end
    clear i j k
% Sorting them

    for i=1:size(pep,1)
        [area{i,1},I]=sort(area{i,1},'descend');
        pep{i,1}=pep{i,1}(:,I);
    end
    clear I i
% 2. Find the length of the shortest segment for cum frequency plot

    for i=1:size(pep,1)
        for j=1:size(pep{i,1},2)
            ln=length(find(pep{i,1}(:,j)<10000));
            if (i+2-ln)==0
                shortseg{i,1}(1,j)=ln;
            elseif (i+2-ln)==1
                l1=find(pep{i,1}(:,j)>=10000,1)-1;
                l2=i+2-l1-1;
                shortseg{i,1}(1,j)=min(l1,l2);
            elseif (i+2-ln)==2
                l1=find(pep{i,1}(:,j)>=10000,1)-1;
                l2=i+2-l1-2;
                shortseg{i,1}(1,j)=min(l1,l2);
            else
                l1=find(pep{i,1}(:,j)>=10000,1)-1;
                l2=i+2-l1-3;
                shortseg{i,1}(1,j)=min(l1,l2);
            end
        end
    end
    clear ln l2 l1 j i 
% 3. Cluster to which it belongs

    for i=1:size(pep,1)
        for j=1:size(pep{i,1},2)
            cluster{i,1}(1,j)=cl;
        end
    end
%% 7. Save relevant data

    catall=[pep area shortseg cluster];
    catall=cell2table(catall,'VariableNames',{'Seq','area','shortseg','Cluster'});
    store(cl).lib=catall;
    %clearvars -except cl store delarea mdmat clustrows
end
%% 8. Sort the generated data
% Merge the data

catdata=[];
for i=1:length(store)
    catdata=[catdata;store(i).lib];
end
%clear store
for i=1:height(catdata)
    sz(i)=size(catdata.Seq{i,1},1);
end
sz=sz';
for i=min(sz):max(sz)
    t=find(sz==i);
    merge(i-2).seq=[catdata.Seq{t}];
    merge(i-2).area=[catdata.area{t}];
    merge(i-2).shortseg=[catdata.shortseg{t}];
    merge(i-2).cluster=[catdata.Cluster{t}];
end
clear catdat i sz t 
%% 
% Sort into its final form:

for i=1:size(merge,2)
        [merge(i).area,I]=sort(merge(i).area,'descend');
        merge(i).seq=merge(i).seq(:,I);
        merge(i).shortseg=merge(i).shortseg(:,I);
        merge(i).cluster=merge(i).cluster(:,I);
end
clear I i
%save merge.mat merge
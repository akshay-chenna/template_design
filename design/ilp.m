%% SET-2 peptide analysis
%% Find the global alignment

clc
clearvars
load top50.mat
set2=string(top50.dat.aa{2,1});
clear top50
for i=1:length(set2)
    for j=i+1:length(set2)
        [~, aln,~]=nwalign(set2(i),set2(j),'ScoringMatrix',eye(24),'gapopen',0);
        fastawrite([num2str(i) '_' num2str(j) '.fasta'],['Peptide_' num2str(i)],aln(1,:));
        fastawrite([num2str(i) '_' num2str(j) '.fasta'],['Peptide_' num2str(j)],aln(3,:));
    end
end
%% 
% Proceed to TM-score calcuations followed by ILP
%% 
%% Integer Linear Programming
% Derive the binary matrix for tm scores and convert to peptides-conformation format

current=pwd;
cd ~/Desktop/origindomain/pepfold_complex/tmscores/cut_35/
selected=readmatrix("../../selected_35.txt","FileType","text","OutputType","double");
tm=[];
for i=1:length(selected)
    temp=readmatrix([num2str(selected(i)) '.txt'],"FileType","text","OutputType","double");
    
    y=reshape(temp,[],5);
    x=[zeros((i-1)*5,5);ones(5,5);y];
    tm=[tm x];
end
tm=tm+tm';
tm(tm==2)=1;
tm(tm<0.5)=0;
tm(tm>=0.5)=1;
cd(current);
clearvars -except tm selected
%% 
% Peptide-conformation format:

peptm=[];
for i=1:length(tm)/5
  temp=any(tm(:,(i-1)*5+1:i*5),2);
  peptm=[peptm temp];
end
clearvars -except peptm selected
% Run ILP

var=peptm.*(-1);
cond(1:size(peptm,1),1)=-1;
lbvar(1:size(peptm,2),1)=0;
ubvar(1:size(peptm,2),1)=1;
bnodes(1:size(peptm,2))=1;
%% actual ilp FUNCTION
%f=bnodes;
f=[10 9 11 10 10 9 8 9 7 8 9 8 12]; %peptide lengths as weights
A=var;
b=cond;
Aeq=[];
beq=[];
intcon=1:size(peptm,2);
lb=lbvar;
ub=ubvar;
options=optimoptions('intlinprog','RootLPAlgorithm','dual-simplex','CutGeneration','intermediate','Heuristics','rss-diving');
[x,fval, exitflag,output]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
peptides=selected(find(x))
% Node graphs

clf

temp1=[zeros(length(peptm),length(peptm)) peptm];
sq_peptm=[temp1; zeros(size(peptm,2),length(temp1))];
o=graph(sq_peptm,'upper');

subnodes=[1:length(peptm) (find(x)+length(peptm))'];
os=subgraph(o,subnodes);
cs=(length(subnodes)-length((find(x)+length(peptm)))+1):length(subnodes);
c=(length(sq_peptm)-size(peptm,2)+1):length(sq_peptm);
xpos1=[1:length(peptm) linspace(10,length(peptm)-10,size(peptm,2))];
ypos1=[ones(1,length(peptm))*1 ones(1,size(peptm,2))*2];
xpos2=xpos1(subnodes);
ypos2=ypos1(subnodes);

hold on
p=plot(o,'XData',xpos1,'YData',ypos1);
p.NodeLabel=[];
for i=1:length(c)
    lbp(i)={strcat('P',num2str(i))};
end
labelnode(p,c,lbp)
p.NodeFontSize=24;
p.NodeColor='k';
p.LineWidth=0.8;
p.EdgeColor=[0.35 0.35 0.35];
p.LineStyle='-';
p.MarkerSize=4;
ps=plot(os,'XData',xpos2,'YData',ypos2);
ps.NodeLabel=[];
ps.NodeColor='k';
ps.LineWidth=2.5;
ps.EdgeColor=[0.2 0.5 0.8];
ps.LineStyle='-';
ps.MarkerSize=4;
set(gca,'visible','off')
hold off
%exportgraphics(p,'ilp.eps','ContentType','vector','BackgroundColor','none')

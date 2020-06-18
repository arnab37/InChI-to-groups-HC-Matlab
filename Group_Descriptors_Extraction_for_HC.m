clc;
clear all;
close all;

T = readcell('Name_Plus_InChI.xlsx');
%%%%%%%%%%%%%%%%%%%%%% extracting only bond connectivity for HC %%%%%%%
expr1='/';
for i=1:size(T,1)
    InChIHC{i}=T{i,2};
    ssC=regexp(InChIHC{i},'C');
    ssH=regexp(InChIHC{i},'H');
    ss1=regexp(InChIHC{i},expr1);
    C1=0;H1=0;
    C1=cell2mat(extractBetween(InChIHC{i},ssC(1)+1,ssH(1)-1));
    
    fid = fopen('C1.txt','w');
    fprintf(fid,'%s\n', C1);
    fclose(fid);
    
    HMtype=load('C1.txt');
    C(i)=C1(1);
    
    Ha1=cell2mat(extractBetween(InChIHC{i},ssH(1)+1,ss1(2)-1));
    
    fid = fopen('Ha1.txt','w');
    fprintf(fid,'%s\n', Ha1);
    fclose(fid);
    
    HMtype=load('Ha1.txt');
    Ha(i)=Ha1(1);
   
    
    pos1=round(ss1(2))+2;
    pos2=round(ss1(3))-1;
    ConnHC{i}=cell(extractBetween(InChIHC{i},pos1,pos2));
    
    pos1=round(ss1(3))+2;
    if(length(ss1)>3)
        pos2=round(ss1(4))-1;
    else
         pos2=numel(InChIHC{i});
    end
    Hattach{i}=cell(extractBetween(InChIHC{i},pos1,pos2));
    D='D';
    if(contains(InChIHC{i},D)==1)
        ss1=regexp(InChIHC{i},D);
        pos1=round(ss1(1))-1;
        if(numel(InChIHC{i})>ss1(end))
            pos2=round(ss1(end))+1;
        else
            pos2=round(ss1(end));
        end
        Dattach{i}=cell(extractBetween(InChIHC{i},pos1,pos2));
    else
        Dattach{i}=[];
    end
    Ebb='b';
    if(contains(InChIHC{i},Ebb)==1)
        sEb1=regexp(InChIHC{i},D);
        pos1=round(sEb1(1))+1;
        pos2=(numel(InChIHC{i}));
        bCis{i}=cell(extractBetween(InChIHC{i},pos1,pos2));
        StrCis{i}=[];
    else
        bCis{i}=[];
        StrCis{i}=[];
    end
    
end
ConnHC=ConnHC';
InChIHC=InChIHC';
Hattach=Hattach';
Dattach=Dattach';
bCis=bCis';
StrCis=StrCis';

CycMolDes = Des1(C,Ha,InChIHC,ConnHC,Hattach,Dattach,bCis,StrCis);


function MolDes = Des1(C,Ha,InChIHC,ConnHC,Hattach,Dattach,bCis,StrCis)

%%%%%%%%%%%% vector for connectivity  %%%%%%%%%%%%%%%%%%
Eh='-';Ep1='(';Ep2=')';Ec=',';m=1;n=1;
for i=1:length(C)
    i
Bc2c=[];DPnum=[];Du=[];CCBond=[];HMtype=[];Repp=[];RS=[];MPComOCP=[];MPair=[];MPair1=[];MPair2=[];Pairs=[];WPair=[];
Ph_num=[];Pp1_num=[];Pp2_num=[];Pc_num=[];Phy=[];PCom=[];Pc_num=[];POp=[];PCp=[];AAA=[];BBB=[];,CCC=[];DDD=[];DDD1=[];
DDD2=[];B1=[];X=[];A1=[];Hout=[];
    Ph=regexp(ConnHC{i},Eh);
    Pp1=regexp(ConnHC{i},Ep1);
    Pp2=regexp(ConnHC{i},Ep2);
    Pc=regexp(ConnHC{i},Ec);
    NN=cellfun('length',ConnHC{i,1});
    if (isempty(Ph)==0)
        Ph_num=cell2mat(Ph);
    else
        Ph_num=[];
    end
    if (isempty(Pp1)==0)
        Pp1_num=cell2mat(Pp1);
    else
        Pp1_num=[];
    end
    if (isempty(Pp2)==0)
        Pp2_num=cell2mat(Pp2);
    else
        Pp2_num=[];
    end
    if (isempty(Pc)==0)
        Pc_num=cell2mat(Pc);
    else
        Pc_num=[];
    end
    AAA=[1 NN];
    
    Phy=Ph_num;
    POp=Pp1_num;
    PCp=Pp2_num;
    PCom=Pc_num;
    
    IdStr=[Ph_num Pp1_num Pp2_num Pc_num AAA];
    
    IdStr=sort(IdStr);
    
    RP=(sum(Pp1_num(:))+sum(Pc_num(:)))/numel(ConnHC{i,1}{1,1});
    
    DPcell_{i}(1)=(extractBetween(ConnHC{i},(IdStr(1)),(IdStr(2)-1)));
    for j=2:length(IdStr)-2
        DPcell_{i}(j)=(extractBetween(ConnHC{i},(IdStr(j)+1),(IdStr(j+1)-1)));
    end
    DPcell_{i}(j+1)=(extractBetween(ConnHC{i},(IdStr(j+1)+1),(IdStr(length(IdStr)))));
    DPcell_{i}=DPcell_{i}'; 
    
    fid = fopen('Dcell.txt','w');
    CT = DPcell_{i};
    fprintf(fid,'%s\n', CT{:});
    fclose(fid);
    

    DPnum=load('Dcell.txt');
    
    %%%%%%%%%%%%%%% DPnum is connectivity vector %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    H='H'; D='D';Eh='-';Ec=',';
    
    sHy=cell2mat(regexp(Hattach{i},Eh)); %%%%% Positions of Hyphen at the H-layer
    sH=cell2mat(regexp(Hattach{i},H));   %%%%% Positions of H at the H-layer
    sC=cell2mat(regexp(Hattach{i},Ec));  %%%%% Positions of Comma at the H-layer
    hM=numel(Hattach{i,1}{1,1});
    
    Hstr=[sC hM+1];
    Hstr=sort((Hstr));
    
    for j=1:length(sH)
        if(j<length(sH))
            HMt{i}(j,1)=(extractBetween(Hattach{i,1}{1,1},(sH(j)+1),(sH(j)+1)));
            if(contains((HMt{i}(j,1)),Ec)==0)
                HMtype1{i}(j,1)=HMt{i}(j,1);
            end
            if(contains((HMt{i}(j,1)),Ec)==1)
                HMtype1{i}(j,1)={'1'};
            end
        end 
        if(j==length(sH))
            if(sH(length(sH))==hM)
                 HMtype1{i}(j,1)={'1'};
            end
            if(sH(length(sH))<hM)
                 HMtype1{i}(j,1)=(extractBetween(Hattach{i,1}{1,1},(sH(j)+1),(sH(j)+1)));
            end
        end    
    end
    
    fid = fopen('HMtype.txt','w');
    CT = HMtype1{i};
    fprintf(fid,'%s\n', CT{:});
    fclose(fid);
    HMtype=load('HMtype.txt');  %%%%% It means the types of carbons, i.e., carbon with 3 hydrogens, or carbon with 2 hydrogen, etc.
    
    Commm=[0 sC];
    
    for j=1:length(sH)   %%%%%% I am getting only one Spl{j} out of this loop
        if(j==1)
            pos1=1;
            Spl11{j}=cell(extractBetween(Hattach{i},pos1,sH(j)-1));
        end
        if(j==2)
            pos1=1;
            Spl11{j}=cell(extractBetween(Hattach{i},sH(1)+2,sH(j)-1));
        end
        if(j==3)
            pos1=1;
            Spl11{j}=cell(extractBetween(Hattach{i},sH(2)+3,sH(j)-1));
        end
    end
    
    
    pp=1;
    for j=1:length(sH)
        Pt=cell(extractBetween(Spl11{j},1,1));
        if(contains(Pt,Ec)==1)
            mm=numel(Spl11{1,j}{1,1});
            Spl{j}=cell(extractBetween(Spl11{j},2,mm));
        end
        if(contains(Pt,Ec)==0)
            mm=numel(Spl11{1,j}{1,1});
            Spl{j}=cell(extractBetween(Spl11{j},1,mm));
        end
        Spl{j};
        TTT=[];AAA=[];BBB=[];CCC=[];DDD=[];XX=[];
        if(contains(Spl{j},Ec)==0) && (contains(Spl{j},Eh)==0)
            fid = fopen('Test.txt','w');
            CT = (Spl{j});
            fprintf(fid,'%s\n', CT{:});
            fclose(fid);
            AAA=load('Test.txt');
        else
            AAA=[];
        end
        AAA;
        if(contains(Spl{j},Ec)==0) && (contains(Spl{j},Eh)==1)
            SSplf=cell2mat(regexp(Spl{j},Eh));
            nn=numel(Spl{1,j}{1,1});
            if(nn==3)
                SSpl{1}=cell(extractBetween(Spl{j},SSplf(1)-1,SSplf(1)-1));
                SSpl{2}=cell(extractBetween(Spl{j},SSplf(1)+1,SSplf(1)+1));
            end
            if(nn==4)
                SSpl{1}=cell(extractBetween(Spl{j},SSplf(1)-1,SSplf(1)-1));
                SSpl{2}=cell(extractBetween(Spl{j},SSplf(1)+1,SSplf(1)+2));
            end
            if(nn==5)
                SSpl{1}=cell(extractBetween(Spl{j},SSplf(1)-2,SSplf(1)-1));
                SSpl{2}=cell(extractBetween(Spl{j},SSplf(1)+1,SSplf(1)+2));
            end
            SSpl;
            fid = fopen('Test.txt','w');
            CT = (SSpl{1});
            fprintf(fid,'%s\n', CT{:});
            fclose(fid);
            B1(1)=load('Test.txt');
            
            fid = fopen('Test1.txt','w');
            CT = (SSpl{2});
            fprintf(fid,'%s\n', CT{:});
            fclose(fid);
            B1(2)=load('Test1.txt');
                l=0;
                for k=B1(1):B1(2)
                    BBB(k)=B1(1)+l;
                    l=l+1;
                end
                l=0;
        else
            BBB=[];
        end
        BBB;
        if(contains(Spl{j},Ec)==1) && (contains(Spl{j},Eh)==0)
            SSplfc1=cell2mat(regexp(Spl{j},Ec));
            mm=numel(Spl{1,j}{1,1});
            SSplfc=[1 SSplfc1 mm];
            SSplfc=sort(SSplfc);
            nn=length(SSplfc)-1;
            SSpl1{1}=cell(extractBetween(Spl{j},SSplfc(1),SSplfc(2)-1));
            
            fid = fopen('Test.txt','w');
            CT = (SSpl1{1});
            fprintf(fid,'%s\n', CT{:});
            fclose(fid);
            xx=load('Test.txt');
            CCC=[xx];
            for k=2:nn-1
                SSpl1{k}=cell(extractBetween(Spl{j},SSplfc(k)+1,SSplfc(k+1)-1));
                CT = (SSpl1{k});
                fid = fopen('Test.txt','w');
                fprintf(fid,'%s\n', CT{:});
                fclose(fid); 
                xx=load('Test.txt');
                CCC=[CCC xx];
            end
            SSpl1{nn}=cell(extractBetween(Spl{j},SSplfc(nn)+1,SSplfc(nn+1)));
            fid = fopen('Test.txt','w');
            CT = (SSpl1{nn});
            fprintf(fid,'%s\n', CT{:});
            fclose(fid);
            xx=load('Test.txt');
            CCC=[CCC xx];
        else
            CCC=[];
        end
        CCC;
        Spl{j};%%% These are the various segmanets in H-layer. These segments are the H-connectivities of C atoms.
        if(contains(Spl{j},Ec)==1) && (contains(Spl{j},Eh)==1)
            SSplfc2=cell2mat(regexp(Spl{j},Ec)); %%%% Position of comma in Spl{j} 
            SSplfh2=cell2mat(regexp(Spl{j},Eh)); %%%% Position of Hyphen in Spl {j}
            mm=numel(Spl{1,j}{1,1}); %%%% length of Spl{j}
            STot1=[1 mm SSplfc2 SSplfh2];
            STot1=sort(STot1);
            STot=STot1;
            
            for x3=1:length(STot1)
                for x4=1:length(SSplfh2)
                    if(STot1(x3)==SSplfh2(x4))
                        STot(x3)=0;
                    end
                end
            end
            
            for x5=1:length(SSplfh2)
                if(SSplfh2(x5)==STot1(2))
                    Sph1=cell(extractBetween(Spl{j},1,SSplfh2(x5)-1));
                    Sph2=cell(extractBetween(Spl{j},SSplfh2(x5)+1,SSplfh2(x5)+2));
                    if(contains(Sph2,Ec)==1)
                        Sph2=cell(extractBetween(Spl{j},SSplfh2(x5)+1,SSplfh2(x5)+1));
                    end
                    fid = fopen('T1.txt','w');
                    CT = (Sph1);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX1=load('T1.txt');
                    fid = fopen('T1.txt','w');
                    CT = (Sph2);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX2=load('T1.txt');
                    XX=[];
                    for x3=XX1:XX2
                        XX(x3)=x3;
                    end
                    DDD=[DDD XX];
                end
                if(SSplfh2(x5)==STot1(end-1))
                    Sph1=cell(extractBetween(Spl{j},SSplfh2(x5)-2,SSplfh2(x5)-1));
                    if(contains(Sph1,Ec)==1)
                        Sph1=cell(extractBetween(Spl{j},SSplfh2(x5)-1,SSplfh2(x5)-1));
                    end
                    Sph2=cell(extractBetween(Spl{j},SSplfh2(x5)+1,STot1(end)));
                    fid = fopen('T1.txt','w');
                    CT = (Sph1);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX1=load('T1.txt');
                    fid = fopen('T1.txt','w');
                    CT = (Sph2);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX2=load('T1.txt');
                    XX=[];
                    for x3=XX1:XX2
                        XX(x3)=x3;
                    end
                    DDD=[DDD XX];
                end    
                if(SSplfh2(x5)~=STot1(end-1))&&(SSplfh2(x5)~=STot1(2))
                    Sph1=cell(extractBetween(Spl{j},SSplfh2(x5)-2,SSplfh2(x5)-1));
                    if(contains(Sph1,Ec)==1)
                        Sph1=cell(extractBetween(Spl{j},SSplfh2(x5)-1,SSplfh2(x5)-1));
                    end
                    Sph2=cell(extractBetween(Spl{j},SSplfh2(x5)+1,SSplfh2(x5)+2));
                    if(contains(Sph2,Ec)==1)||(contains(Sph2,Eh)==1)
                        Sph2=cell(extractBetween(Spl{j},SSplfh2(x5)+1,SSplfh2(x5)+1));
                    end
                    fid = fopen('T1.txt','w');
                    CT = (Sph1);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX1=load('T1.txt');
                    fid = fopen('T1.txt','w');
                    CT = (Sph2);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX2=load('T1.txt');
                    XX=[];
                    for x3=XX1:XX2
                        XX(x3)=x3;
                    end
                    DDD=[DDD XX];
                end
            end
            
            
            
            for x5=1:length(SSplfc2)
                if(SSplfc2(x5)==STot1(2))
                    Spc1=cell(extractBetween(Spl{j},1,SSplfc2(x5)-1));
                    fid = fopen('T1.txt','w');
                    CT = (Spc1);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX1=load('T1.txt');
                    DDD=[DDD XX1];
                end
                if(SSplfc2(x5)==STot1(end-1))
                    Spc1=cell(extractBetween(Spl{j},SSplfc2(x5)+1,STot1(end)));
                    fid = fopen('T1.txt','w');
                    CT = (Spc1);
                    fprintf(fid,'%s\n', CT{:});
                    fclose(fid);
                    XX1=load('T1.txt');
                    DDD=[DDD XX1];
                end    
                if(SSplfc2(x5)~=STot1(2))
                    Spl{j};
                    SSplfc2(x5)-3;
                    SSplfc2(x5)-1;
                    STot1;
                    Spc2=cell(extractBetween(Spl{j},SSplfc2(x5)-3,SSplfc2(x5)-1));
                    if(contains(Spc2,Eh)==0)
                       Spc1=cell(extractBetween(Spl{j},SSplfc2(x5)-2,SSplfc2(x5)-1));
                       if(contains(Spc1,Ec)==1)
                          Spc1=cell(extractBetween(Spl{j},SSplfc2(x5)-1,SSplfc2(x5)-1));
                       end
                        fid = fopen('T1.txt','w');
                        CT = (Spc1);
                        fprintf(fid,'%s\n', CT{:});
                        fclose(fid);
                        XX1=load('T1.txt');
                        DDD=[DDD XX1];
                    end
                end
            end
            
        else
            
            DDD=[];
            
        end
       
        TTT=[AAA BBB CCC DDD];
        TTT=nonzeros(TTT);
        for k=1:length(TTT)
           Hout(pp,1)=TTT(k);
           Hout(pp,2)=HMtype(j);
           pp=pp+1;
        end
    end
   
    DPAll=zeros(length(DPnum),2);
    DPAll(:,1)=DPnum(:);

    for j=1:length(DPnum)
        for k=1:length(Hout)
            if((DPnum(j)==Hout(k,1)))
                DPAll(j,2)=Hout(k,2);
            end
        end
    end
   %%% DPAll 1st col is the sequential numbers in the connectivty, 2nd
   %%% column is the number of H or D attached with each number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Bond connectivity atom to atom %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Operation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% connectivty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% with
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% paranthesis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%
    MPair=[]; POp;
    i;

    if(isempty(POp)==0)
        
        ppp=1;
    for j=1:length(POp) %%%%%%% this fucntion is equivalent to Pairs =[POp PCp]
        for k=1:length(PCp)
            if (POp(j)<PCp(k))
                Pairs(ppp,1)=POp(j);
                Pairs(ppp,2)=PCp(k);
                ppp=ppp+1;
            end
        end 
    end
    Pairs;
    ppp=ppp-1; %%% ppp is the number of rows in the Pairs matrix
    %%%%%%%%%%% the Pairs is the matrix of all Open + close parenthesis
    %%%%%%%%%%% combination, not sure whether they are the correct
    %%%%%%%%%%% combination.. if the paranthesis in the connectivity layer
    %%%%%%%%%%% is ()(), i.e., not nested then Pairs is correct. But if the
    %%%%%%%%%%% paranthesis is nested (()) then the combination is not
    %%%%%%%%%%% correct..
    
   %%%%%%%%%%%%%% Start of the loop for extracting correct paranthesis
   %%%%%%%%%%%%%% pair's index
   %%%%%%%%%%%%%% 
   ppo=1;
    for j=1:ppp
        k1=0;
        for k=1:length(POp)  %%%% ppp = length of POp or PCp
            if (Pairs(j,1)<POp(k)) && (Pairs(j,2)>POp(k))
                k1=k1+1;
            end
        end
        
        
        %%%%%%%%%%%%% k1 measures the count of nested open parenthesis, such as the open paranthesis
        %%%%% at the middle of (()...
        
        k2=0;
        for k=1:length(PCp)
            if (Pairs(j,1)<PCp(k)) && (Pairs(j,2)>PCp(k))
                k2=k2+1;
            end
        end
        
        
        %%%%%%%% k2 measures the count of nested close paranthesis, such as the open paranthesis
        %%%%% at the middle of ())...
        %%%%%%
        
        %%%%%%%%%%%%%%% currently k1 and k2 will be eqaul to zero if there is any correct pair exists, such as
        %%%%%%%%%%%%%%%% (())().... if you see the last pair then for the
        %%%%%%%%%%%%%%%% last pair k1=0 and k2=0; thus it is a correct
        %%%%%%%%%%%%%%%% pair.
        if(k1==0)&&(k2==0) 
            MPair1(ppo,1)=Pairs(j,1);
            MPair1(ppo,2)=Pairs(j,2);
            ppo=ppo+1;
        end
    end
    ppo=ppo-1;  %%% ppo is the number of rows of Mpair1 right now.
    
    %%%%%%% MPair1 contains the those pairs from the Pairs maxtrix which has the zeros of ( and ) inside them %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% More refinement is required to remove large unncessary pairs
    %%%%%%%% %%%%%
    MPair1;
    %%%% Wrong pairs for now, WPair = Pairs - MPair1 %%%%%%%%%%
    if(length(POp)==size(MPair1,1))
        MPair=MPair1;
    end
    
    if(length(POp)>size(MPair1,1))
        
    POpm=[];POpm1=[];PCpm=[];PCpm1=[];
    POpm=POp;PCpm=PCp;
    for jk=1:length(POp)
        for lm=1:size(MPair1,1)
            if(POp(jk)==MPair1(lm,1))
                POpm(jk)=0;
            end
        end
    end
    
    for jk=1:length(PCp)
        for lm=1:size(MPair1,1)
            if(PCp(jk)==MPair1(lm,2))
                PCpm(jk)=0;
            end
        end
    end
    
    PCpm1=nonzeros(PCpm);POpm1=nonzeros(POpm);
    
    ppp=1;
    for j=1:length(POpm1) %%%%%%% this fucntion is equivalent to Pairs =[POp PCp]
        for k=1:length(PCpm1)
            if (POpm1(j)<PCpm1(k))
                WPair(ppp,1)=POpm1(j);
                WPair(ppp,2)=PCpm1(k);
                ppp=ppp+1;
            end
        end 
    end
    
    Pairs;
    MPair1;
    WPair;
    MPair2 = PairParanthesisCorrection(WPair);
    x2=1;
    for jj=1:size(MPair1,1)
        MPair(x2,:)=MPair1(jj,:);
        x2=x2+1;
    end
    for jj=1:size(MPair2,1)
        MPair(x2,:)=MPair2(jj,:);
        x2=x2+1;
    end
    end
    end
    
    MPair;
    
    %%%%%%%%%%%%%%%%% End of the actual paranthesis finding %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Insertion of comma starts from the next line %%%%%%%
     
    MPComOCP=zeros(length(PCom),3);  %% MPComOCP = Modified Pairs of Comma, Openned, Closed Paranthesis %%%
    
    if (isempty(MPair)==0)
    for j=1:length(PCom)
        mmm=1;
        for k=1:size(MPair,1)
            if(PCom(j)>MPair(k,1)) && (PCom(j)<MPair(k,2))
                Checking(mmm,1)=MPair(k,1);
                Checking(mmm,2)=MPair(k,2);
                Checking(mmm,3)=MPair(k,2)-MPair(k,1);
                mmm=mmm+1;
            end
        end
        Mini=min(Checking(:,3));
        for l=1:size(Checking,1)
            if(Checking(l,3)==Mini)
                MPComOCP(j,1)=Checking(l,1);
                MPComOCP(j,2)=Checking(l,2);
                MPComOCP(j,3)=PCom(j);
            end
        end
    end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of operation with parenthesis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in conncetivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%% Line 468 to Line 720 is to decode the bond
    %%%%%%%%%%%%%%%%% connectivity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DPnum;
    %%%%%%%%%%%%%%%%%%%%  Recognition of the starting and ending of each number %%%%%%%%%%%%%%%%%
     pk=1;
    Symbl=[Phy POp PCp PCom];
    Symbl=sort(Symbl);
    for s1s=1:length(DPnum) %%%% the end is in 1044 line %%%%%
        PpSs=num2str(DPnum(s1s)); %%% PpSs is the # of a carbon atom in string
        PS=[];
        PS11=cell2mat(regexp(ConnHC{i},PpSs));  %%% Positions of that number in a string
        mm=numel(ConnHC{i,1}{1,1});
        xx1=1;
        strlength(PpSs);
        if(strlength(PpSs)==1)  %%%%% correction of the position of 1 digit index %%%%
            for jj=1:length(PS11)
                if (PS11(jj)>1)&&(PS11(jj)<mm)
                    yT11=cell(extractBetween(ConnHC{i},PS11(jj)-1,PS11(jj)-1));
                    yT12=cell(extractBetween(ConnHC{i},PS11(jj)+1,PS11(jj)+1));
                    if((contains(yT11,Eh)==1)||(contains(yT11,Ep1)==1)||(contains(yT11,Ep2)==1)||(contains(yT11,Ec)==1))&&((contains(yT12,Eh)==1)||(contains(yT12,Ep1)==1)||(contains(yT12,Ep2)==1)||(contains(yT12,Ec)==1))
                        PS(xx1)=PS11(jj);
                        xx1=xx1+1;
                    end
                end
                if (PS11(jj)==1)
                    yT12=cell(extractBetween(ConnHC{i},PS11(jj)+1,PS11(jj)+1));
                    if((contains(yT12,Eh)==1)||(contains(yT12,Ep1)==1)||(contains(yT12,Ep2)==1)||(contains(yT12,Ec)==1))
                        PS(xx1)=PS11(jj);
                        xx1=xx1+1;
                    end
                end
                if (PS11(jj)==mm)
                    yT12=cell(extractBetween(ConnHC{i},PS11(jj)-1,PS11(jj)-1));
                    if((contains(yT12,Eh)==1)||(contains(yT12,Ep1)==1)||(contains(yT12,Ep2)==1)||(contains(yT12,Ec)==1))
                        PS(xx1)=PS11(jj);
                        xx1=xx1+1;
                    end
                end
                
            end
        end
        xx1=1;
      if(strlength(PpSs)==2)  %%%%% correction of the position of 2 digit index %%%%
          for jj=1:length(PS11)
              PS(xx1)=PS11(jj);
              xx1=xx1+1;
          end
      end
          
          
          
    
      for jj=1:length(PS)
            
            %%%%%%%%%%%%%%%%%%%%  BACKWARD PASS i.e., towards Left from the atom %%%%%%%%%%%%%%%%%
        if (PS(jj)~=1) 
            yT=cell(extractBetween(ConnHC{i},PS(jj)-1,PS(jj)-1)); 
                %%%%% Just previous symbol %%%
           for mnop=1:1 
            if(contains(yT,Eh)==1) && ((PS(jj)-1)==Symbl(1)) %%%% if it contains - then %%%%
                y1=char(extractBetween(ConnHC{i},1,PS(jj)-2));
                PpSs;
                fid = fopen('T1.txt','w');
                CT = (PpSs);
                fprintf(fid,'%s\n', CT);
                fclose(fid);
                X(pk,1)=load('T1.txt');
                y1;
                fid = fopen('T1.txt','w');
                CT = (y1);
                fprintf(fid,'%s\n', CT);
                fclose(fid);
                X(pk,2)=load('T1.txt');
                pk=pk+1;
            end
            if(contains(yT,Eh)==1) && ((PS(jj)-1)~=Symbl(1))   
               PrevS=cell(extractBetween(ConnHC{i},PS(jj)-3,PS(jj)-3)); %%%% other cases with single digits
               if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                  y1=char(extractBetween(ConnHC{i},PS(jj)-2,PS(jj)-2));
                  PpSs;
                  fid = fopen('T1.txt','w');
                  CT = (PpSs);
                  fprintf(fid,'%s\n', CT);
                  fclose(fid);
                  X(pk,1)=load('T1.txt');
                  y1;
                  fid = fopen('T1.txt','w');
                  CT = (y1);
                  fprintf(fid,'%s\n', CT);
                  fclose(fid);
                  X(pk,2)=load('T1.txt');
                  pk=pk+1;
               end
               PrevS=cell(extractBetween(ConnHC{i},PS(jj)-4,PS(jj)-4));
               if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                  y1=char(extractBetween(ConnHC{i},PS(jj)-3,PS(jj)-2));
                  PpSs;
                  fid = fopen('T1.txt','w');
                  CT = (PpSs);
                  fprintf(fid,'%s\n', CT);
                  fclose(fid);
                  X(pk,1)=load('T1.txt');
                  y1;
                  fid = fopen('T1.txt','w');
                  CT = (y1);
                  fprintf(fid,'%s\n', CT);
                  fclose(fid);
                  X(pk,2)=load('T1.txt');
                  pk=pk+1;
               end
            end
          end
        
          if(contains(yT,Ec)==1) %%%% if it contains , then %%%%
             for k=1:size(MPComOCP,1)        %%%%%%%%%%%%% if there is comma then MPComOCP and MPair wil not be empty
                 if(MPComOCP(k,3)==(PS(jj)-1))&&(MPComOCP(k,1)==Symbl(1))
                    y1=char(extractBetween(ConnHC{i},1,MPComOCP(k,1)-1));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');                       
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end
                 if(MPComOCP(k,3)==(PS(jj)-1))&&(MPComOCP(k,1)~=Symbl(1))
                 PrevS=cell(extractBetween(ConnHC{i},MPComOCP(k,1)-2,MPComOCP(k,1)-2));  
                    if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                       y1=char(extractBetween(ConnHC{i},MPComOCP(k,1)-1,MPComOCP(k,1)-1));
                       PpSs;
                       fid = fopen('T1.txt','w');
                       CT = (PpSs);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,1)=load('T1.txt'); 
                       y1;
                       fid = fopen('T1.txt','w');
                       CT = (y1);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,2)=load('T1.txt');
                       pk=pk+1;
                    end
                    PrevS=cell(extractBetween(ConnHC{i},MPComOCP(k,1)-3,MPComOCP(k,1)-3));
                    if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                       y1=char(extractBetween(ConnHC{i},MPComOCP(k,1)-2,MPComOCP(k,1)-1));
                       PpSs;
                       fid = fopen('T1.txt','w');
                       CT = (PpSs);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,1)=load('T1.txt');                    
                       y1;
                       fid = fopen('T1.txt','w');
                       CT = (y1);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,2)=load('T1.txt');
                       pk=pk+1;
                    end
                 end
             end
          end
                                             
          if(contains(yT,Ep2)==1) %%%% if it contains ) then %%%%
             for k=1:size(MPair,1)
                 if(MPair(k,2)==(PS(jj)-1))&&(MPair(k,1)==Symbl(1))
                    y1=char(extractBetween(ConnHC{i},1,MPair(k,1)-1));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');    
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end
                 if(MPair(k,2)==(PS(jj)-1))&&(MPair(k,1)~=Symbl(1))
                    PrevS=cell(extractBetween(ConnHC{i},MPair(k,1)-2,MPair(k,1)-2));            
                    if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                       y1=char(extractBetween(ConnHC{i},MPair(k,1)-1,MPair(k,1)-1));
                       PpSs;
                       fid = fopen('T1.txt','w');
                       CT = (PpSs);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,1)=load('T1.txt'); 
                       y1;
                       fid = fopen('T1.txt','w');
                       CT = (y1);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,2)=load('T1.txt');
                       pk=pk+1;
                    end
                    PrevS=cell(extractBetween(ConnHC{i},MPair(k,1)-3,MPair(k,1)-3));  
                    if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                       y1=char(extractBetween(ConnHC{i},MPair(k,1)-2,MPair(k,1)-1));
                       PpSs;
                       fid = fopen('T1.txt','w');
                       CT = (PpSs);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
                       X(pk,1)=load('T1.txt');  
                       y1;
                       fid = fopen('T1.txt','w');
                       CT = (y1);
                       fprintf(fid,'%s\n', CT);
                       fclose(fid);
%                        ConnHC{i}
%                        PpSs
%                        y1
%                        MPair(k,1)-1
                       X(pk,2)=load('T1.txt');
                       pk=pk+1;
                    end
                 end
             end
          end           
                
          if(contains(yT,Ep1)==1) %%%% if it contains( then %%%%
             for k=1:size(MPair,1)
                 if(MPair(k,1)==(PS(jj)-1))&&(MPair(k,1)==Symbl(1))
                    y1=char(extractBetween(ConnHC{i},1,PS(jj)-2));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end
                 if(MPair(k,1)==(PS(jj)-1))&&(MPair(k,1)~=Symbl(1))
                    PrevS=cell(extractBetween(ConnHC{i},MPair(k,1)-2,MPair(k,1)-2));                                
                    if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                      y1=char(extractBetween(ConnHC{i},MPair(k,1)-1,MPair(k,1)-1));
                      PpSs;
                      fid = fopen('T1.txt','w'); 
                      CT = (PpSs);
                      fprintf(fid,'%s\n', CT);
                      fclose(fid);
                      X(pk,1)=load('T1.txt');                      
                      y1;
                      fid = fopen('T1.txt','w');
                      CT = (y1);
                      fprintf(fid,'%s\n', CT);
                      fclose(fid);
                      X(pk,2)=load('T1.txt');
                      pk=pk+1;
                    end
                    PrevS=cell(extractBetween(ConnHC{i},MPair(k,1)-3,MPair(k,1)-3));     
                    if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                      y1=char(extractBetween(ConnHC{i},MPair(k,1)-2,MPair(k,1)-1));
                      PpSs;
                      fid = fopen('T1.txt','w');
                      CT = (PpSs);
                      fprintf(fid,'%s\n', CT);
                      fclose(fid);
                      X(pk,1)=load('T1.txt');
                      y1;
                      fid = fopen('T1.txt','w');
                      CT = (y1);
                      fprintf(fid,'%s\n', CT);
                      fclose(fid);
                      X(pk,2)=load('T1.txt');
                      pk=pk+1;
                    end
                 end
             end
          end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FORWAD PASS, i.e, Right from the atom %%%%%%%%%%
        
        if(strlength(PpSs)==2)
            aa=PS(jj)+2;
            %%%%%% aa reads the postion next to a number
        end
        if(strlength(PpSs)==1)
            aa=PS(jj)+1; 
            %%%%% aa reads the postion next to a number
        end  
        
        if (aa<=Symbl(end))
            yT=cell(extractBetween(ConnHC{i},aa,aa));
        
        
        for mnop=1:1    
        if(contains(yT,Eh)==1)&&(aa==Symbl(end))  %%%%%%%%%%%%%%%%%%% it is for - %%%
           mm=numel(ConnHC{i,1}{1,1});
           y1=char(extractBetween(ConnHC{i},aa+1,mm));
           PpSs;
           fid = fopen('T1.txt','w');
           CT = (PpSs);
           fprintf(fid,'%s\n', CT);
           fclose(fid);
           X(pk,1)=load('T1.txt');           
           y1;
           fid = fopen('T1.txt','w');
           CT = (y1);
           fprintf(fid,'%s\n', CT);
           fclose(fid);
           X(pk,2)=load('T1.txt');
           pk=pk+1;
        end
        if(contains(yT,Eh)==1)&&(aa<Symbl(end))
        PrevS=cell(extractBetween(ConnHC{i},aa+2,aa+2));
          if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
            y1=char(extractBetween(ConnHC{i},aa+1,aa+1));
            PpSs;
            fid = fopen('T1.txt','w');
            CT = (PpSs);
            fprintf(fid,'%s\n', CT);
            fclose(fid);
            X(pk,1)=load('T1.txt');
            y1;
            fid = fopen('T1.txt','w');
            CT = (y1);
            fprintf(fid,'%s\n', CT);
            fclose(fid);
            X(pk,2)=load('T1.txt');
            pk=pk+1;
          end
          PrevS=cell(extractBetween(ConnHC{i},aa+3,aa+3));
          if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
            y1=char(extractBetween(ConnHC{i},aa+1,aa+2));
            PpSs;
            fid = fopen('T1.txt','w');
            CT = (PpSs);
            fprintf(fid,'%s\n', CT);
            fclose(fid);
            X(pk,1)=load('T1.txt');
            y1;
            y1;
            fid = fopen('T1.txt','w');
            CT = (y1);
            fprintf(fid,'%s\n', CT);
            fclose(fid);
            X(pk,2)=load('T1.txt');
            pk=pk+1;
          end
         end
        end

                  
        if(contains(yT,Ep1)==1) %%%% Contains ( 
          for k=1:size(MPair,1)
              if(aa==MPair(k,1))&&(MPair(k,2)==Symbl(end))
                    y1=char(extractBetween(ConnHC{i},MPair(k,2)+1,numel(ConnHC{i,1}{1,1})));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');    
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
              end
              if(aa==MPair(k,1))&&(MPair(k,2)~=Symbl(end))
                 PrevS=cell(extractBetween(ConnHC{i},MPair(k,2)+2,MPair(k,2)+2));
                 if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                    y1=char(extractBetween(ConnHC{i},MPair(k,2)+1,MPair(k,2)+1));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');    
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end
                 PrevS=cell(extractBetween(ConnHC{i},MPair(k,2)+3,MPair(k,2)+3));
                 if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                    y1=char(extractBetween(ConnHC{i},MPair(k,2)+1,MPair(k,2)+2));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');  
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end    
                 PrevS=cell(extractBetween(ConnHC{i},aa+2,aa+2));
                 if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                    y1=char(extractBetween(ConnHC{i},aa+1,aa+1));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');  
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end           
                 PrevS=cell(extractBetween(ConnHC{i},aa+3,aa+3));
                 if(contains(PrevS,Eh)==1)||(contains(PrevS,Ep1)==1)||(contains(PrevS,Ep2)==1)||(contains(PrevS,Ec)==1)
                    y1=char(extractBetween(ConnHC{i},aa+1,aa+2));
                    PpSs;
                    fid = fopen('T1.txt','w');
                    CT = (PpSs);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,1)=load('T1.txt');   
                    y1;
                    fid = fopen('T1.txt','w');
                    CT = (y1);
                    fprintf(fid,'%s\n', CT);
                    fclose(fid);
                    X(pk,2)=load('T1.txt');
                    pk=pk+1;
                 end
              end
          end
        end
        
        if(contains(yT,Ec)==1)              
           for l=1:size(MPComOCP,1)
             if(MPComOCP(l,2)~=Symbl(end))
               PrevS1=cell(extractBetween(ConnHC{i},MPComOCP(l,3)+2,MPComOCP(l,3)+2));
               PrevS2=cell(extractBetween(ConnHC{i},MPComOCP(l,1)-2,MPComOCP(l,1)-2));
               if(contains(PrevS1,Eh)==1)||(contains(PrevS1,Ep1)==1)||(contains(PrevS1,Ep2)==1)||(contains(PrevS1,Ec)==1)
                 if(contains(PrevS2,Eh)==1)||(contains(PrevS2,Ep1)==1)||(contains(PrevS2,Ep2)==1)||(contains(PrevS2,Ec)==1)
                   y1=char(extractBetween(ConnHC{i},MPComOCP(l,3)+1,MPComOCP(l,3)+1));
                   y2=char(extractBetween(ConnHC{i},MPComOCP(l,1)-1,MPComOCP(l,1)-1));
                   y1;
                   fid = fopen('T1.txt','w');
                   CT = (y1);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,1)=load('T1.txt');     
                   y2;
                   fid = fopen('T1.txt','w');
                   CT = (y2);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,2)=load('T1.txt');
                   pk=pk+1;
                 end
               end
               PrevS1=cell(extractBetween(ConnHC{i},MPComOCP(l,3)+3,MPComOCP(l,3)+3));
               PrevS2=cell(extractBetween(ConnHC{i},MPComOCP(l,1)-2,MPComOCP(l,1)-2));
               if(contains(PrevS1,Eh)==1)||(contains(PrevS1,Ep1)==1)||(contains(PrevS1,Ep2)==1)||(contains(PrevS1,Ec)==1)
                 if(contains(PrevS2,Eh)==1)||(contains(PrevS2,Ep1)==1)||(contains(PrevS2,Ep2)==1)||(contains(PrevS2,Ec)==1)
                   y1=char(extractBetween(ConnHC{i},MPComOCP(l,3)+1,MPComOCP(l,3)+2));
                   y2=char(extractBetween(ConnHC{i},MPComOCP(l,1)-1,MPComOCP(l,1)-1));
                   y1;
                   fid = fopen('T1.txt','w');
                   CT = (y1);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,1)=load('T1.txt');     
                   y2;
                   fid = fopen('T1.txt','w');
                   CT = (y2);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,2)=load('T1.txt');
                   pk=pk+1;
                 end
               end
               PrevS1=cell(extractBetween(ConnHC{i},MPComOCP(l,3)+3,MPComOCP(l,3)+3));
               PrevS2=cell(extractBetween(ConnHC{i},MPComOCP(l,1)-3,MPComOCP(l,1)-3));
               if(contains(PrevS1,Eh)==1)||(contains(PrevS1,Ep1)==1)||(contains(PrevS1,Ep2)==1)||(contains(PrevS1,Ec)==1)
                 if(contains(PrevS2,Eh)==1)||(contains(PrevS2,Ep1)==1)||(contains(PrevS2,Ep2)==1)||(contains(PrevS2,Ec)==1)
                   y1=char(extractBetween(ConnHC{i},MPComOCP(l,3)+1,MPComOCP(l,3)+2));
                   y2=char(extractBetween(ConnHC{i},MPComOCP(l,1)-2,MPComOCP(l,1)-1));
                   y1;
                   fid = fopen('T1.txt','w');
                   CT = (y1);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,1)=load('T1.txt');     
                   y2;
                   fid = fopen('T1.txt','w');
                   CT = (y2);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,2)=load('T1.txt');
                   pk=pk+1;
                 end
               end                   
               PrevS1=cell(extractBetween(ConnHC{i},MPComOCP(l,3)+2,MPComOCP(l,3)+2));
               PrevS2=cell(extractBetween(ConnHC{i},MPComOCP(l,1)-3,MPComOCP(l,1)-3));
               if(contains(PrevS1,Eh)==1)||(contains(PrevS1,Ep1)==1)||(contains(PrevS1,Ep2)==1)||(contains(PrevS1,Ec)==1)
                 if(contains(PrevS2,Eh)==1)||(contains(PrevS2,Ep1)==1)||(contains(PrevS2,Ep2)==1)||(contains(PrevS2,Ec)==1)
                   y1=char(extractBetween(ConnHC{i},MPComOCP(l,3)+1,MPComOCP(l,3)+1));
                   y2=char(extractBetween(ConnHC{i},MPComOCP(l,1)-2,MPComOCP(l,1)-1));
                   y1;
                   fid = fopen('T1.txt','w');
                   CT = (y1);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,1)=load('T1.txt');
%                    ConnHC{i}
%                    PrevS1
%                    PrevS2
%                    y1
%                    y2
                   y2;
                   fid = fopen('T1.txt','w');
                   CT = (y2);
                   fprintf(fid,'%s\n', CT);
                   fclose(fid);
                   X(pk,2)=load('T1.txt');
                   pk=pk+1;
                 end
               end
           end
        end
        end
        
        end
      end
      
    end
    
    

 
    %%%%%%%%%%%%%%%%% Line 468 to Line 980 is to decode the bond
    %%%%%%%%%%%%%%%%% connectivity,  is the matrix where thefirst column is the atom position and second column is 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% is the atoms connected to the atom at the
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% first column %%%%
    
    

    BC2C=X;
    X=[];
    DPnum;
    BC2C;
    Bc2c=bugfix(BC2C);
    BC2C=[];
    %%%%%%%%%%%%%%%%%%%%%% Bc2c is the all bond connectivity %%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    Du=unique(DPnum);
    
    for k=1:length(Du)
        l=0;
        for jj=1:size(Bc2c,1)
            if(Du(k)==Bc2c(jj,1))
                l=l+1;
            end
        end
    CCBond(k,1)=Du(k);
    CCBond(k,2)=l;
    end
    CCBond;
    TBond=zeros(length(Du),2);
    TBond(:,:)=CCBond(:,:);
    Hout;
    for k=1:length(Du)
        for l=1:size(Hout,1)
            if(Du(k)==Hout(l,1))
                TBond(k,2)=CCBond(k,2)+Hout(l,2);
            end
        end
    end
    CCBond;
    TBond;
    A1=zeros(size(Bc2c,1),5);
    A1(:,1)=Bc2c(:,1);
    A1(:,4)=Bc2c(:,2);
    Bc2c;
    for k=1:size(A1,1)
        for l=1:length(Du)
            if(A1(k,1)==Du(l))
                A1(k,2)=TBond(l,2);
                A1(k,3)=CCBond(l,2);
            end
            if(A1(k,4)==Du(l))
                A1(k,5)=TBond(l,2);
            end
        end
    end
    TBond;
    CCBond;
    for jj=1:size(A1,1)
        for k=1:size(TBond,1)
            if(A1(jj,1)==TBond(k,1))
                A1(jj,2)=TBond(k,2);
            end
            if(A1(jj,4)==TBond(k,1))
                A1(jj,5)=TBond(k,2);
            end
        end
    end
    A2=[];
    A2=A1(:,1:3);
       for jj=1:size(A2,1)-1
       for kk=jj+1:size(A2,1)
           if(A2(jj,1)==A2(kk,1))&&(A2(jj,2)==A2(kk,2))&&(A2(jj,3)==A2(kk,3))&&(jj~=kk)&&((A2(jj,1)~=0))&&((A2(jj,2)~=0))&&((A2(jj,2)~=0))
                A2(jj,1)=0;
                A2(jj,2)=0;
                A2(jj,3)=0;
           end
       end
   end
   A2;
   A3=[];
   pp=1;
   for jj=1:1:size(A2,1)
       if(sum(A2(jj,:))~=0)
           A3(pp,1)=A2(jj,1);
           A3(pp,2)=A2(jj,2);
           A3(pp,3)=A2(jj,3);
           pp=pp+1;
       end
   end
   A3;
    %%%%%%%%%%%%%%%%%%%%%%%% Cycles finding %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Rings = phenyl and others %%%%%%%%%%%%%%%%%%%%
    SR5=nonzeros(DPnum);
    Ncount = histc(SR5, unique(SR5));
    Sr5u=unique(sort(SR5));Repp=[];
    C11=[Sr5u Ncount];
    for k=1:size(C11,1)
        if (C11(k,2)==1)
            C11(k,1)=0;
            C11(k,2)=0;
        end
    end
    C11;
    Repp=unique(nonzeros(C11(:,1)));SmallRingSize=0;LargeRingSize=0;C6H5=0;
    C6H4=0; C6H3=0; C6H2=0; C6H=0;C6=0;
    
    Repp;
    
    if (isempty(Repp)==0)
    Cyc4=fun(Bc2c, Repp);

        for iii=1:size(Cyc4,1)
            Sr4=Cyc4(iii,:);
            RS(iii)=length(unique(nonzeros(Sr4)));
        end
        SmallRingSize=min((RS));
        LargeRingSize=max((RS));
        
        for jj=1:length(RS)
            mp=0;mH=0;
            if (RS(jj)==6)
                Sr7=nonzeros(unique(Cyc4(jj,:)));
                xx11=0;
                for kk=1:length(Sr7)
                    for ll=1:size(A3,1)
                        if(Sr7(kk)==A3(ll,1))&&(A3(ll,2)==3)
                            xx11=xx11+1;
                        end
                    end
                end
                if(xx11==6)
                    for kk=1:length(Sr7)
                        for ll=1:size(A3,1)
                            if(Sr7(kk)==A3(ll,1))&&(A3(ll,2)==3)
                               A3(ll,2)=5;
                            end
                        end
                    end
                end
                
             end
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% End of the ring business %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Descriptors entry %%%%%%%%%%%%%%%%%%%%%%

   SB1=0;SB2=0;SB3=0;SB4=0;DB1=0;DB2=0;DB3=0;BB1=0;BB2=0;TB1=0;TB2=0;A4=[];xx=1;
   for jj=1:size(A3,1)
       if (A3(jj,2)==4)
           if(A3(jj,3)==1)
               SB1=SB1+1; %% C(C)(H)3
           end
           if(A3(jj,3)==2)
               SB2=SB2+1; %% C(C)2(H)2
           end
           if(A3(jj,3)==3)
               SB3=SB3+1; %% C(C)3(H)
           end
           if(A3(jj,3)==4)
               SB4=SB4+1; %% C(C)4
           end
       end
       
       if (A3(jj,2)==3)
           if(A3(jj,3)==1)
               DB1=DB1+1; %% Cd(C)(H)2
           end
           if(A3(jj,3)==2)
               DB2=DB2+1; %% Cd(C)2(H)
           end
           if(A3(jj,3)==3)
               DB3=DB3+1; %% Cd(C)2
           end
       end
       
       if (A3(jj,2)==5)
           if(A3(jj,3)==2)
               BB1=BB1+1; %% CB(CB)2(H)
           end
           if(A3(jj,3)==3)
               BB2=BB2+1; %% CB(CB)2(C)
           end
       end
       
       if (A3(jj,2)==2)
           if(A3(jj,3)==1)
               TB1=TB1+1; %% Ct(Ct)(H)
           end
           if(A3(jj,3)==2)
               TB2=TB2+1; %% Ct(Ct)(C)+ Ca(Cd)2
               A4(xx)=A3(jj,1);
               xx=xx+1;
           end
       end 
   end
       
           Ca=0;Ct=0;
   if (isempty(A4)==0)
       Sr8=[];xx=1;
       for jj=1:length(A4)
           for kk=1:size(Bc2c)
               if(A4(jj)==Bc2c(kk,1))
                   Sr8(xx)=Bc2c(kk,2); %%% connectd carbon 
                   xx=xx+1;
               end
           end
           Sr9=[];xx=1;
           for kk=1:length(Sr8)
               for ll=1:size(A3,1)
                   if (Sr8(kk)==A3(ll,1))
                       Sr9(xx,1)=A3(ll,2);  %%% Total bonds of connected carbon
                       Sr9(xx,2)=A3(ll,3);  %%% CC bond of connected carbon
                       xx=xx+1;
                   end
               end
           end
           
           if((Sr9(1,1)==3)&&((Sr9(1,2)==1)||(Sr9(1,2)==2)||(Sr9(1,2)==3)))&&((Sr9(2,1)==3)&&((Sr9(2,2)==1)||(Sr9(2,2)==2)||(Sr9(2,2)==3)))
               Ca=Ca+1;
           end
           if((Sr9(1,1)==3)&&((Sr9(1,2)==1)||(Sr9(1,2)==2)||(Sr9(1,2)==3)))&&((Sr9(2,1)==2)&&((Sr9(2,2)==2)))
               Ca=Ca+1;
           end
           if((Sr9(1,1)==2)&&((Sr9(1,2)==2)))&&((Sr9(2,1)==2)&&((Sr9(2,2)==2)))
               Ca=Ca+1;
           end
       end
   end
       
       Ct=TB2-Ca; %%% Ct(Ct)(C)
       Ca; %%% Ca(C)2
           
           
    %%%%%%%%%%%%%%%% positions of double and triple bond %%%%%%%
    pL=0;
    for k=1:length(DPnum)
        for l=1:size(TBond,1)
            if(DPnum(k)==TBond(l,1))
                if(TBond(l,2)==2)||(TBond(l,2)==3)
                    pL=pL+k;
                end
            end
        end
    end
    
    RelPosDTBond=pL/length(DPnum);
    
    %%%%%%%%%%%%%%%%%%%% finding cis count %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cis=0;
    minus='-';
    
    if(isempty(bCis{i})==0)
        mm=numel(bCis{i,1}{1,1});
        Ptt=cell(extractBetween(bCis{i},mm,mm));
        if(contains(Ptt,minus)==1)
            cis=cis+1;
        end
        if(contains(bCis{i},Ec)==1)
            sCC=cell2mat(regexp(bCis{i},Ec));
            for jj=1:length(sCC)
                Ptt=cell(extractBetween(bCis{i},sCC(jj)-1,sCC(jj)-1));
                if(contains(Ptt,minus)==1)
                    cis=cis+1;
                end
            end
        end
    end

                
    %%%%%%%%%%%%%%%%%%%%%%%%%% DPH contains number of H connected to each %
    %%%%%%%%%%%%%%%%%%%%%%%%%% C atom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dtotal=0;
    if(isempty(Dattach{i})==0)
        sD1=cell2mat(regexp(Dattach{i},D));
        sC1=cell2mat(regexp(Dattach{i},Ec));
        Dstr=[sC1  numel(Dattach{i})+1];
        Dstr=sort(Dstr);
        for j=1:length(sD1)
            if((sD1(j)+1)==(Dstr(j)-1))
                nD1{i}(j,1)=cell(extractBetween((Dattach{i}),(sD1(j)+1),(Dstr(j)-1)));
            end
            if((sD1(j)+1)>(Dstr(j)-1))
                nD1{i}(j,1)={'1'};
            end   
        end
    
    fid = fopen('nD.txt','w');
    CT = nD1{i};
    fprintf(fid,'%s\n', CT{:});
    fclose(fid);
    nD=load('nD.txt');
    
    Dtotal=sum(nD(:));
    end
    
    %%%%%%%%%%%%%%%%%%%%% Descriptors %%%%%%%%%%%%%%%%%%%
    
    MolDes(i,1)=SB1;SB1=0; %% Cs(C)(H)3
    MolDes(i,2)=SB2;SB2=0; %% Cs(C)2(H)2
    MolDes(i,3)=SB3;SB3=0; %% Cs(C)3(H)
    MolDes(i,4)=SB4;SB4=0; %% Cs(C)4
    MolDes(i,5)=DB1;DB1=0; %% Cd(C)(H)2
    MolDes(i,6)=DB2;DB2=0; %% Cd(C)2(H)
    MolDes(i,7)=DB3;DB3=0; %% Cd(C)2
    MolDes(i,8)=BB1;BB1=0; %% CB(CB)2(H)
    MolDes(i,9)=BB2;BB2=0; %% CB(CB)2(C)
    MolDes(i,10)=TB1;TB1=0; %% Ct(Ct)(H)
    MolDes(i,11)=Ct;Ct=0; %% Ct(Ct)(C)
    MolDes(i,12)=Ca;Ca=0; %% Ca(C)2
    MolDes(i,13)=RelPosDTBond;RelPosDTBond=0; % Relative positions of double and triple bonds %%%
    MolDes(i,14)=cis;cis=0; %% Cis count
    MolDes(i,15)=RP;RP=0; %% Realtive positions of branching 
    MolDes(i,16)=SmallRingSize;SmallRingSize=0; %% Smallest Ring know to create more ring strain except cyc-hexane
    MolDes(i,17)=LargeRingSize;LargeRingSize=0;
    MolDes(i,18)=C(i);  %% Total carbon atoms
    MolDes(i,19)=Ha(i)-Dtotal;  %% Total H atoms
    MolDes(i,20)=Dtotal;Dtotal=0;  %% Total D atoms
    
    A1=[];A2=[];A3=[];A4=[];Repp=[];Cyc4=[];

end
end

function ring= fun(Bc2c, Repp)
A=Bc2c;
A2=Bc2c;
R=Repp(1);

a=max(A(:,1));

Cyc=zeros(2*a,a+1);

for i=1:size(Cyc,1)
    Cyc(i,1)=R;
end

for i=1:size(Cyc,1)
    B(:)=A(:,1);
    C(:)=A(:,2);
    B1(:)=A2(:,2);
    C1=A2(:,1);
    for j=2:size(Cyc,2)
        p=0;
        for k=1:length(B)
            if(p==0)&&(Cyc(i,j-1)==B(k))&&(j==2)
                p=p+1;
                Cyc(i,j)=C(k);
                A(k,1)=0;
                A(k,2)=0;
            end
            if(p==0)&&(Cyc(i,j-1)==B(k))&&(j>2)&&(Cyc(i,j-2)~=C(k))
                p=p+1;
                Cyc(i,j)=C(k);
                A(k,1)=0;
                A(k,2)=0;
            end
        end
        for k=1:length(B1)
            if(p==0)&&(Cyc(i,j-1)==B1(k))&&(j==2)
                p=p+1;
                Cyc(i,j)=C1(k);
                A2(k,1)=0;
                A2(k,2)=0;
            end
            if(p==0)&&(Cyc(i,j-1)==B1(k))&&(j>2)&&(Cyc(i,j-2)~=C1(k))
                p=p+1;
                Cyc(i,j)=C1(k);
                A2(k,1)=0;
                A2(k,2)=0;
            end
        end
    end
end
Cyc;
Cyc1=zeros(size(Cyc,1),size(Cyc,2));
p=1;
for i=1:size(Cyc,1)
    Sr=nonzeros(Cyc(i,:));
    
    if(length(unique(nonzeros(Sr)))<length(nonzeros(Sr)))
        for j=1:length(Sr)
             Cyc1(p,j)=Sr(j);
        end
             p=p+1;
    end
end
p=1;
Cyc1;
for i=1:size(Cyc1,1)
    if (sum(Cyc1(i,:)~=0))
        Cyc21(p,:)=Cyc1(i,:);
        p=p+1;
    end
    
end
Cyc21;
for i=1:size(Cyc21,1)
    Sr2=Cyc21(i,:);
    for j=2:length(Sr2)
        p=0;
        if(p==0)&&(Sr2(1)==Sr2(j))
            for k=j+1:length(Sr2)
                Sr2(k)=0;
                Cyc21(i,k)=0;
            end
            p=p+1;
        end
    end
end

p=1;

for i=1:size(Cyc21,1)-1
    for j=i+1:size(Cyc21,1)
        if(isempty(setxor(Cyc21(i,:),Cyc21(j,:)))==1)
            Sr3=Cyc21(i,:);
            for k=1:length(Sr3)
                Cyc21(i,k)=0;
            end
        end
    end
end

p=1;
for i=1:size(Cyc21,1)
    if (sum(Cyc21(i,:)~=0))
        Cyc2(p,:)=Cyc21(i,:);
        p=p+1;
    end
    
end    
Cyc2;

p=1;
for i=1:size(Cyc2,1)
    Sr5=nonzeros(Cyc2(i,:));
    Ncount = histc(Sr5, unique(Sr5));
    Sr5u=unique(sort(Sr5));
    C1=[Sr5u Ncount];
    for k=1:size(C1,1)
        if (C1(k,2)==1)
            C1(k,1)=0;
            C1(k,2)=0;
        end
    end
    C=unique(nonzeros(C1(:,1)));
    
    if(length(C)==1)&&(C(1)==Cyc2(i,1))
        Cyc3(p,:)=Cyc2(i,:);
        p=p+1;
    end
    
    if((length(C)==1)&&(C(1)~=Cyc2(i,1)))||((length(C)>1)&&((C(1)~=Cyc2(i,1))||(C(2)~=Cyc2(i,1))))
        if(length(C)==1)
            a=C(1);
        end
        if(length(C)>1)
            if (C(1)~=Cyc2(i,1))
                a=C(1);
            else
                a=C(2);
            end
        end
        
        P=1;
        for j=1:size(Cyc2,2)
            Sr5d(j)=Cyc2(i,j)-a;
            if(P<3)&&(Sr5d(j)==0)
                m(P)=j;
                P=P+1;
            end
        end
        for j=1:m(1)-1
                Cyc2(i,j)=0;
        end
        for j=m(2)+1:size(Cyc2,2)
                Cyc2(i,j)=0;
        end
        Cyc3(p,:)=Cyc2(i,:);
        p=p+1;
    end
end

ring(:,:)=Cyc3(:,:);
end 

function Bc2c = bugfix(BC2C)
BC2C;
A=BC2C;
 
 p=1;
for j=1:size(A,1)
     for i=1:size(A,1)
         if(((A(j,1)==A(i,2))&&(A(j,2)==A(i,1)))||((A(j,1)==A(i,1))&&(A(j,2)==A(i,2)))&&(i~=j)&&((A(j,1)~=0)||(A(j,2)~=0)))
             A(j,1)=0;
             A(j,2)=0;
         end
     end
end

p=1;
for i=1:size(A,1)
    if(sum(A(i,:))~=0)
        A1(p,:)=A(i,:);
        p=p+1;
    end
end

A2=zeros(2*size(A1,1),2);

for i=1:size(A1,1)
    A2(i,:)=A1(i,:);
    A2(size(A1,1)+i,1)=A1(i,2);
    A2(size(A1,1)+i,2)=A1(i,1);
end

Bc2c(:,:)=A2(:,:);
end


function MPair2 = PairParanthesisCorrection(WPair)
WPair;

A=WPair(:,1);POp=A;
B=WPair(:,2);PCp=B;
L=size(WPair,1);

x=1;
for k=1:L
    if(isempty(A)==0)&&(isempty(B)==0)
        p=1;
        for i=1:length(A)
            for j=1:length(B)
                if(A(i)<B(j))
                    Diff(p,1)=A(i);
                    Diff(p,2)=B(j);
                    Diff(p,3)=B(j)-A(i);
                    p=p+1;
                end
            end
        end
        WPair;
        Diff;
        p1=0;
        for i=1:size(Diff,1)
            if(p1==0)&&(Diff(i,3)==min(Diff(:,3)))
                p1=p1+1;
                MPair2(x,1)=Diff(i,1);
                MPair2(x,2)=Diff(i,2);
                for j=1:length(POp)
                    if(POp(j)==MPair2(x,1))
                        POp(j)=0;
                    end
                end
                for j=1:length(PCp)
                    if(PCp(j)==MPair2(x,2))
                        PCp(j)=0;
                    end
                end
                x=x+1;
            end
        end
        A=sort(nonzeros(POp));
        B=sort(nonzeros(PCp));
        Diff=[];
    end
end
end

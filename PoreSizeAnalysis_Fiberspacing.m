clear all
% close all
% Analysis based on David Weitz method on distribution of fiber spacings
% load a binary image


[filename, filepath] = uigetfile('*.tif'); 
        fname = [filepath filename];
        
        
InfoImage=imfinfo(fname);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
Image=zeros(nImage,mImage,NumberImages,'uint16');

% Pixel size [um]
px_size=1/InfoImage(1).XResolution;


TifLink = Tiff(fname, 'r');
    for i=1:NumberImages
   TifLink.setDirectory(i);
   Image(:,:,i)=TifLink.read();
    end
    
TifLink.close();
  
slice_counter=1;
z=1;

FinalImage=Image(:,:,z);

figure
hold on
imagesc(FinalImage)
colormap('gray')
title({filename},'FontSize',12)
xlim([0 mImage])
ylim([0 nImage])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'DataAspectRatio',[1 1 1]);
axis('ij');
hold off



FinalImage(FinalImage==255)=1;

% RowToAnalyze=FinalImage;  % (y-direction)
RowToAnalyze=FinalImage'; % (x-direction)

% RowToAnalyze=RowToAnalyze';
% 

%

i=1;

%  This part is adapted from binary time trace analysis, thus the variable
%  names hav "time"

for w=1:numel(RowToAnalyze(1,:)) 
    
On_diff=diff([0; RowToAnalyze(:,w); 0] == 1);
Off_diff=diff([1; RowToAnalyze(:,w); 1] == 0);

% eval(['Onbins' num2str(th) '=find(On_diff == -1)- find(On_diff == 1);'])
% eval(['Offbins' num2str(th) '=find(Off_diff == -1)- find(Off_diff == 1);'])

OnbinsStart=find(On_diff == 1);
OnbinsEnd=find(On_diff == -1);

OffbinsStart=find(Off_diff == 1);
OffbinsEnd=find(Off_diff == -1);

Onbins=OnbinsEnd-OnbinsStart;
Offbins=OffbinsEnd-OffbinsStart;

% Onbins=find(On_diff == -1)- find(On_diff == 1);
% Offbins=find(Off_diff == -1)- find(Off_diff == 1);

% [Offbins_int,~,~]=find(RowToAnalyze(:,w)==0);
% [Onbins_int,~,~]=find(RowToAnalyze(:,w));

% Ontimes=Onbins.*px_size;
% Offtimes=Offbins.*px_size;


% Ontimes_int(numel(Ontimes),1)=0;
% for l=1:numel(Ontimes)
% Ontimes_int(l,:)=(sum(bt(OnbinsStart(l):1:OnbinsEnd(l)))-bt(OnbinsEnd(l)))/Onbins(l);
% end
% 
% Offtimes_int(numel(Offtimes),1)=0;
% for l=1:[numel(Offtimes)-1]
% Offtimes_int(l,:)=(sum(bt(OffbinsStart(l):1:OffbinsEnd(l)))-bt(OffbinsEnd(l)))/Offbins(l);
% end
% Offtimes_int(end,:)=(sum(bt(OffbinsStart(numel(Offtimes)):1:numel(bt)))-bt(end))/(Offbins(end)-1);

% bt=bt';
% Offbins= find(Off_diff == -1)- find(Off_diff == 1);

OnOffTh(i).Ontimes=Onbins;
OnOffTh(i).Offtimes=Offbins;

% Dwell-time intensities

% OnOffTh(i).Offtimes_int=Offtimes_int.*1000; %convert to counts per second
% % OnOffTh(i).Ontimes_int=Ontimes_int.*1000;

% Time trace length

mean_Off(i)=mean(Offbins);
mean_On(i)=mean(Onbins);

% mean_Off_int(i)=mean(Offtimes_int);
% mean_On_int(i)=mean(Ontimes_int);

median_Off(i)=median(Offbins);
median_On(i)=median(Onbins);

% k(i)=numel(Ontimes)/(sum(Offtimes)+sum(Ontimes)); % Mean turnover rate (1/s)
% k(i)=numel(Ontimes)/(sum(Offtimes)+sum(Ontimes)); % Mean turnover rate (1/s) mod. 14.3.2015

i=i+1;

clear Ontimes_int Offtimes_int

end

% Combine fiber and void distances from different slices

%
All_Fiberspace = cell2mat({OnOffTh.Ontimes}');
All_Voidspace = cell2mat({OnOffTh.Offtimes}');
% unit still in pixels! conversion to real units as last step

% Dwelltime histogram

binsize=5;

OffHist=hist(All_Voidspace,[binsize:binsize:max(All_Voidspace)]);
OffHist=OffHist';


% %normalized histograms
OffHist_norm=OffHist./numel(All_Voidspace);


xaxis=[binsize:binsize:max(All_Voidspace)]'; %x values for histogram
x_axis_units=xaxis.*px_size;
% xaxis=xaxis';

mean_Off_hist=mean(All_Voidspace);

mean_Off_hist_units=mean(All_Voidspace).*px_size
% 
% figure
% axes('YScale','log','YMinorTick','on',...
%     'FontSize',12);
% hold on



% % Create semilogy  
% semilogy(x_axis_units,OffHist_norm,'.k','MarkerSize',10);
% xlabel('\xi [\mum]','FontSize',14);
% ylabel('P(\xi)','FontSize',14);
% title({filename},'FontSize',12) 
% hold off


[fitresult_single, gof_single] = createFit_Single_exp(x_axis_units, OffHist_norm)
coeffvalues=coeffvalues(fitresult_single);
p_single=coeffvalues;

% poresize_single=(-1/p_single(2))*px_size

poresize_single=(-1/p_single(2))

clear coeffvalues

[fitresult_double, gof_double] = createFit_double_exp(x_axis_units, OffHist_norm)
coeffvalues=coeffvalues(fitresult_double);
p_double=coeffvalues;

% pore1_double=(-1/p_double(2))*px_size
% pore2_double=(-1/p_double(4))*px_size

pore1_double=(-1/p_double(2))
pore2_double=(-1/p_double(4))

PoreCombined_double=(pore1_double*pore2_double)/(pore1_double+pore2_double)

fracAmp1=abs(p_double(1))/(abs(p_double(1))+abs(p_double(3)))
fracAmp2=abs(p_double(3))/(abs(p_double(1))+abs(p_double(3)))


%% Fit Off-time histograms

clear cu
figure
axes1 = axes('YScale','log','YMinorTick','on',...
    'FontSize',12);
hold(axes1,'all');
colo='grcmbkkkkk';
semilogy(x_axis_units,OffHist_norm,'.k','MarkerSize',10);
st=1; %to be able to ignore part of the curve in fitting...

% st_end=6;
 st_end=numel(OffHist_norm);

err=zeros(5,1);
w=1;
for i=[5]
   
    [p,err(i),cu{w},R2]=fitDwelltimehist_poremod2022(xaxis(st:st_end),OffHist_norm(st:st_end),i);
    semilogy(x_axis_units(st:st_end),cu{w},colo(i),'LineWidth',2);
    disp([sprintf('Model %d: R2=%g; Parameters=[',i,R2), sprintf(' %g', p),' ]']) ;
    w=w+1;
end

% st=20; %to be able to ignore part of the curve in fitting...
% st_end=Offbins;
% err=zeros(5,1);
% for i=1:3
%     [p,err(i),cu,R2]=fitDwelltimehist(xaxis(st:st_end),OffHist_norm(st:st_end),i);
%     semilogy(xaxis(st:st_end),cu,colo(i),'LineWidth',2);
%     disp([sprintf('Model %d: R2=%g; Parameters=[',i,R2), sprintf(' %g', p),' ]']) ; 
% end
Poresize=p(2).*px_size
%      ylim([0 1.2]) % Petri
% title(sample,'FontSize',14) % Petri
xlabel('\xi [\mum]','FontSize',14);
ylabel('P(\xi)','FontSize',14);


 legend('data',sprintf('single exponent fit: Mean poresize = %.3f um',Poresize),'double exp')
 
%  xlim([0 0.4])
%  ylim([0.0001 0.1])title({filename},'FontSize',12) 
title({filename},'FontSize',12) 
hold off


%%





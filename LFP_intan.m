% input 0.1-10kHz 
% output spikes


%% Loading the mat file
close all;clear all;clc;  
read_Intan_RHD2000_file


% [time,act_ch,data,aux] = read_intan_data;
%%
num=3; %channel 번호(9~12)
% amplifier_channels(num) %임피던스 확인, 숫자넣기
%% 부분추출 그래프그리기
f=frequency_parameters.amplifier_sample_rate;
% second_start= 0.1%시작하는 초
% second_end= 60 %끝나는 초
% plot(t_amplifier(1,f*second_start:f*second_end), amplifier_data(num,f*second_start:f*second_end)); %신호 확인, 숫자넣기
% 
% %% rms 구하기
% y = rms( amplifier_data(num,:));
%     timeline = 0:(1/Fs):Trec-1/Fs;


%% 이건 data전체일때
Trec = round(length(amplifier_data(num,:))/f); % total recording time
    set(0,'defaultfigurepos',[100 200 2000 400]);

figure, plot(t_amplifier, amplifier_data(num,:)); %신호 확인, 숫자넣기
%  xlim([130 131]);
% ylim([-200 200]);

%axis([1 2 -800 800]);xlabel('second(s)');ylabel('uV');title('Channel 10')
%%
% wavclus_handles_AM;

 data = double(amplifier_data(num,:))';

%%
% [m,n]=size(data);
% data=reshape(data,[1,m]);


    % Filtering (low pass filter for LFP)
%     Fc = 300; %300
%     [b,a] = butter(4,Fc/(f/2),'low');
%     lpf_data= filtfilt(b,a,data);
%     clear a b
        
    Fc = [1000 4000];
    [b,a] = butter(4,Fc/(f/2),'stop');
    lpf_data= filtfilt(b,a,data);
%     clear a b
    
    
    % Desampling 1000 Hz
    dFs = 3000;
    lpf_data_de= resample(lpf_data,dFs,f);
    x_de= 1/dFs:1/dFs:length(lpf_data_de)/dFs;

    % STFT 생성
    R= 500; %500                 % welch window size=  1초
    noverlap = R/2;          % 50 % overlap (1.25로 나누면 좀더 촘촘해짐)
    window1= hanning(R);
    NFFT =2048*4;              %2^nextpow2(L);
    [S, F, T]= spectrogram(lpf_data_de, window1, noverlap, NFFT, dFs);
%     clear R window1 noverlap NFFT
%                     figure(), clf
%                     colorbar = [0 2000];    % color map을 조절하며 원하는 scale 찾음
%                     imagesc(T,F,abs(S),colorbar);
%                     axis xy
%                     xlabel('time(s)'), ylabel('frequency(Hz)'), title('');
%                     ylim([0 Fc]), xlim([1 Trec])
%     
    % Normalized STFT
    before_fin = floor(double(length(T)*(20/Trec)));
    before_stft = abs(S(:,1:before_fin));
    during_after_stft= abs(S(:,1:end));
    
    baseline_value= mean(before_stft,2);
    mul= ones(1,length(T));
    baseline_matrix= baseline_value*mul;
    
    % Normalizaed: (X-baseline)/baseline
    normalized_stft= (during_after_stft-baseline_matrix)./baseline_matrix;
    
%     figure(), \%clf
    set(0,'defaultfigurepos',[100 100 2000 400])
    figure
    set(gca, 'FontSize',12)
    colorbar= [0 1];
    imagesc((T(1:end)),F,normalized_stft,colorbar);%
    axis xy
    xlabel('time(s)'), ylabel('frequency(Hz)'),
%     xlim([60 120]);
    ylim([800 1100]);
    colorbar;

    %% Normalized STFT smooting
    
    PSF = fspecial('average',[8 16]);
    Blurred_n = imfilter(normalized_stft,PSF,'conv');
    
    figure
    
    imagesc((T(1:end)),F,Blurred_n,[-0.5 0.5]);
    
    axis xy
    xlabel('time(s)','FontSize',10), ylabel('frequency(Hz)','FontSize',10)
    
    set(gca, 'FontSize',10),xlim([0 Trec]),ylim([0 150])
    hold on;
    
    set(gcf, 'Position', [100 300 2000 400]); % graph size 고정
    yL = get(gca,'YLim');
%     line([60 60],yL,'Color','r');
%     line([80 80],yL,'Color','r');

                
                %%
%                 
Gamma_power_s1 = zeros(1,1679);
Theta_power_s1 = zeros(1,1679);
          %% power-time plot
            % total 1-150 Hz
                F_index=find(F>1&F<150);
                F_index2= find(F>1&F<4);
                F_index3= find(F>5&F<8); 
                F_index4= find(F>9&F<12); 
                F_index5= find(F>13&F<25); 
                F_index6= find(F>26&F<150); 
%                 Wn=F_index;
%                 Wn_2=F_index2;
%                 Wn_3=F_index3;
%                 Wn_4=F_index4;
%                 Wn_5=F_index5;
%                 Wn_6=F_index6;
%                 
                toal_power=mean(Blurred_n(F_index,:));

                figure,             
                subplot(6,1,1), plot(T,toal_power,'k'),
                set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); 
                xlim([0 Trec]),ylim([-1 2.5]), title('Total power (1-150 Hz)');

                delta_power=mean(Blurred_n(F_index2,:));
                subplot(6,1,2), plot(T,delta_power,'k'),
                set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); 
                xlim([0 Trec]),ylim([-1 2.5]), title('Delta power(1-4 Hz)');

                theta_power=mean(Blurred_n(F_index3,:));
                subplot(6,1,3), plot(T,theta_power,'k'),
                set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); 
                xlim([0 Trec]),ylim([-1 2.5]), title('Theta power(5-8 Hz)');
                
                alpha_power=mean(Blurred_n(F_index4,:));
                subplot(6,1,4), plot(T,alpha_power,'k'),
                set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); 
                xlim([0 Trec]),ylim([-1 2.5]), title('Alpha power(9-12 Hz)');
                
                beta_power=mean(Blurred_n(F_index5,:));
                subplot(6,1,5), plot(T,beta_power,'k'),
                set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); 
                xlim([0 Trec]),ylim([-1 2.5]), title('Beta power(13-25 Hz)');
                
                gamma_power=mean(Blurred_n(F_index6,:));
                subplot(6,1,6), plot(T,gamma_power,'k'),
                set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); 
                xlim([0 Trec]),ylim([-1 2.5]), title('Gamma power(26-150 Hz)');
                
%                 baseline_mean= mean(Gamma_power(:,1:78)'); % 1439*(20/360)
%                 stdz=std(Gamma_power(1,1:78));
                %high_standard= baseline_mean+ 1*stdz;
                %low_standard= baseline_mean- 1*stdz;
%                 zz= 0:0.1:Trec; 
%                 %HS= high_standard*ones(size(zz));
%                 %LS= low_standard*ones(size(zz));
%                 hold on; plot(zz, baseline_mean,'b')
%                 %hold on; plot(zz, HS,'b:')
%                 %hold on; plot(zz, LS,'b:')
%                 set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); xlabel('time (sec)');
%                 xlim([0 Trec]),ylim([-.5 1]), title('Gamma power');
%                 yL = get(gca,'YLim');
%                 line([20 20],yL,'Color','r');
%                 line([40 40],yL,'Color','r');                
            
% %%
%             % Theta 4-9 Hz
%                 F_index4= find(F>1&F<11);
%                 Wn2=F_index4;
%                 Theta_power=mean(Blurred_n(Wn2,:));
% 
% %                 figure
%                  subplot(2,1,2),plot(T,Theta_power,'k')
%                 
%                 baseline_mean2= mean(Theta_power(:,1:78)'); % 1439*(20/360)
%                 %stdz=std(Theta_power(1,1:78));
%                 %high_standard= baseline_mean+ 1*stdz;
%                 %low_standard= baseline_mean- 1*stdz;
%                 zz= 0:0.1:Trec; 
%                 %HS= high_standard*ones(size(zz));
%                 %LS= low_standard*ones(size(zz));
%                 hold on; plot(zz, baseline_mean2,'b')
%                 %hold on; plot(zz, HS,'b:')
%                 %hold on; plot(zz, LS,'b:')
%                 set(gca, 'FontSize',10);
%                 ylabel('Spectral power'); xlabel('time (sec)'), title('Theta power');
%                 xlim([0 Trec]),ylim([-.5 1])
%                 yL = get(gca,'YLim');
% %                 line([20 20],yL,'Color','r');
% %                 line([40 40],yL,'Color','r');                
     
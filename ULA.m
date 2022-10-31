clear all;
close all;

load("spcom_10sep.mat");
doa_est_plot(X,Delta,M,N);
load("spcom_50sep.mat");
doa_est_plot(X,Delta,M,N);

% Function to plot DOA estimate using
function [] = doa_est_plot(X,Delta,M,N)
    theta = -90:90;
    theta_rad = theta*pi/180;
    Rx = zeros([M M]); 
    
    %% Calculating sample covariance matrix
    for i = 1:M
        for j = 1:M
            Rx(i,j) = (1/N)*X(i,:)*X(j,:)';
        end
    end
    
    
    %% Calculating the spatial power spectrum using classical beamform
    spect_est_bartlett = zeros(length(theta_rad),1);
    for i = 1:length(theta_rad)
        a = sv(theta_rad(i),M,Delta);
        spect_est_bartlett(i) = (a'*Rx*a)/(a'*a);
    end
    
    Pbartlett = abs(spect_est_bartlett);
    Pbartlett_db = 10*log10(abs(spect_est_bartlett));
    figure;
    plot(theta,Pbartlett_db);
    hold on;
    
    %% Calculating the spatial power spectrum using MVDR
    spect_est_mvdr = zeros(length(theta_rad),1);
    for i = 1:length(theta_rad)
        a = sv(theta_rad(i),M,Delta);
        spect_est_mvdr(i) = 1/(a'/Rx*a);
    end
    
    Pmvdr = abs(spect_est_mvdr);
    Pmvdr_db = 10*log10(abs(spect_est_mvdr));
    plot(theta,Pmvdr_db);
    
    
    
    %% Calculating the spatial power spectrum using MUSIC
    
    [U,Lamda] = eig(Rx);
    Lamda = abs(diag(Lamda));
    Un = [];
    for i = 1:length(Lamda)
        if Lamda(i) < 0.01
            Un = [Un, U(:,i)];
        end
    end
    
    
    spect_est_music = zeros(length(theta_rad),1);
    for i = 1:length(theta_rad)
        a = sv(theta_rad(i),M,Delta);
        spect_est_music(i) = 1/(a'*(Un*Un')*a);
    end
    
    Pmusic = abs(spect_est_music);
    Pmusic_db = 10*log10(abs(spect_est_music));
    plot(theta,Pmusic_db);
    

    %% Estimation of DOA using ESPRIT

   [U,Lamda] = eig(Rx);
   Lamda = diag(Lamda);
   nsources = 0;
   for i = 1:length(Rx)
       if (Lamda(i) > 0.01)
           nsources = nsources + 1;
       end
   end

   Us = U(:,1:nsources);
   Phi = Us(1:end-1,:)\Us(2:end,:);
   Lamda = eig(Phi);
   DoA = -asind(angle(Lamda)/(2*pi*Delta));



   plot(DoA,zeros(length(DoA),1),'o')
   title('Direction of Arrival Estimation')
   xlabel("Angle (degrees)")
   ylabel("Power (dB)")
   legend("Bartlett","MVDR","MUSIC","ESPRIT")

end

%function to calculate the steering vector
function a = sv(theta, M, Delta)
    a = exp(2j*pi*Delta*(0:M-1)*sin(theta))';
end


clear all;
close all;

load("spcom_10sep.mat");
doa_est_plot(X,Delta,M,N);
load("spcom_50sep.mat");
doa_est_plot(X,Delta,M,N);

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
    
    Pbartlett = 10*log10(abs(spect_est_bartlett));
    figure;
    plot(theta,Pbartlett);
    hold on;
    
    %% Calculating the spatial power spectrum using MVDR
    spect_est_mvdr = zeros(length(theta_rad),1);
    for i = 1:length(theta_rad)
        a = sv(theta_rad(i),M,Delta);
        spect_est_mvdr(i) = 1/(a'/Rx*a);
    end
    
    Pmvdr = 10*log10(abs(spect_est_mvdr));
    plot(theta,Pmvdr);
    
    
    
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
        spect_est_music(i) = (a'*a)/(a'*(Un*Un')*a);
    end
    
    Pmusic = 10*log10(abs(spect_est_music));
    plot(theta,Pmusic);
    

    %% Estimation of DOA using ESPRIT

   [U,~] = eig(Rx);
   Us = U(:,1:2);
   Phi = Us(1:end-1,:)\Us(2:end,:);
   Lamda = eig(Phi);
   DoA = -asind(angle(Lamda)/(2*pi*Delta));



   plot(DoA,zeros(length(DoA),1),'o')
   title('Direction of Arrival Estimation')
   xlabel("Angle (degrees)")
   ylabel("Power (dB)")
   legend("Bartlett","MVDR","MUSIC","ESPRIT")

end

%% Functions
%function to calculate the steering vector
function a = sv(theta, M, Delta)
    a = exp(2j*pi*Delta*(0:M-1)*sin(theta))';
end


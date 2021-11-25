%%% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo da sequência positiva

% ENTRADA COM DADOS DE ÂNGULO

% phix - ângulos 
% f0 - frequência nominal
% Fs - Frequência de amostragem


function [freq_final] = estimador_freq_deriv_ang_Ain(phi0, phi1, phi2, f0, h, p0, p1, p2,N,m,Fs)

num_wins = length(phi0);

% desvio_fase = zeros(3,num_wins-1);
% freq_meas = ones(3,num_wins-1);
% freq_final = ones(1,num_wins)*ref_seg(1,1);% ***** cria vetor c/ dados de freq de ref *******


dt = m/Fs;                % passo da janela (s)
ajuste = pi/90;

% % Cria ângulo constante de referência - para janelamento 'm' < N 
% % (desvio angular dif de zero mesmo com freq nominal)
num_cic = f0*N/12800;
Ang_ref = (m/N)*2*pi*num_cic;
%Ang_ref = 0;

for i=1:num_wins
    
    desvio_fase(1,i) = phi0(i) - p0(i) - Ang_ref;
    desvio_fase(2,i) = phi1(i) - p1(i) - Ang_ref;
    desvio_fase(3,i) = phi2(i) - p2(i) - Ang_ref;

    % Corrige desvio ângular maior que 180º
    for j=1:3
        % Enquanto houver desvios >ou< que +/-180º executa correção
        while desvio_fase(j,i)>pi || desvio_fase(j,i)<-pi
            if desvio_fase(j,i)>pi
                desvio_fase(j,i) = desvio_fase(j,i)-2*pi;
            elseif desvio_fase(j,i)<-pi
                desvio_fase(j,i) = desvio_fase(j,i)+2*pi;
            end
        end
%         % Enquanto houver desvios >ou< que 6º executa correção
%         while desvio_fase(j,i)>ajuste || desvio_fase(j,i)<-ajuste
%             if desvio_fase(j,i)>ajuste
%                 desvio_fase(j,i) = ajuste;%desvio_fase(j,i)-2*pi;
%             elseif desvio_fase(j,i)<-ajuste
%                 desvio_fase(j,i) = -ajuste;
%             end
%         end
        % Finalmente, calcula freq. com desvios corrigidos
        freq_meas(j,i) = h*f0 + desvio_fase(j,i)/(2*pi*dt); 
    end         
    freq_final(i) = mean(freq_meas); % media das tres fases
end

end
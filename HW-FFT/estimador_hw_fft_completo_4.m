


function [Amp_a,Amp_b,Amp_c,ang_a,ang_b,ang_c,freq_final] = estimador_hw_fft_completo_4(Va,Vb,Vc,t,refs,N,m,Fs,f0,noise,param,type)



%-------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)

[~, ~, ~, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);
num_cic = length(ref_seg(1,:));




num_h = 1;
Tm = 0.02*(m/N);
Tin = 0:Tm:(num_cic*Tm);


%Fe = ones(num_cic,num_h);
Fe(1,:) = ref_seg(1,:)'*num_h;
Fes(1,:) = ref_seg(1,:)'*num_h;

%Fem = ones(num_cic,num_h);
Fem(:,1) = ref_seg(1,:)*num_h';

T = 0.02;     % tempo total (em segundos)
f = ref_seg(1,1);
%-------------------------------------------------------------------------
for w=1:num_cic
    
    if w==1
        Fs = round(ref_seg(1,1)*N);
    else
        Fs = round(Fes(w-1)*N);   % ITERA VETOR DE FREQUÊNCIAS
    end
        
    for h=1:num_h  % itera componentes harmonicas


        freqs(1,h) = f*h;        % Vetor de frequências

        phases(1,h) = deg2rad(0);
        phases(2,h) = deg2rad(-120);
        phases(3,h) = deg2rad(+120);

        amps(1,h) = 0.1;
        amps(2,h) = 0.1;
        amps(3,h) = 0.1;
    end

    amps(1,1) = 1;
    amps(2,1) = 1;
    amps(3,1) = 1;


    % ------------------------------------------------------------------------------------------------------------------
    % GERAÇÃO DO SINAL
    [Va, Vb, Vc, t, f_ref, a_ref, p_ref, mod] = gerador_sinais_Tin(freqs, phases, amps, Fs, noise, type, param, N,Tin(w));
    
%     figure, hold on, grid on, stem(Va)
    %figure, hold on, grid on, stem(a)
    
    % %-------------------------------------------------------------------------------------------------------------------
    % % SEGMENTAÇÃO (JANELAMENTO)
    % % Cria ângulo constante de referência - para janelamento 'm' < N 
    % % (desvio angular dif de zero mesmo com freq nominal)
    num_cic = f0*N/12800;
    Ang_ref = (m/N)*2*pi*num_cic;
    %Ang_ref = 0;
    
    % REFERÊNCIAS
    for i=1:(T/0.02)

        for fase=1:3
            amp_ref(:,:,fase) = a_ref(:,(i-1)*N+1,fase)';
            phi_ref(:,:,fase) = p_ref(:,(i-1)*N+1,fase)';  

        end
        freq_ref = f_ref(:,(i-1)*N+1)';

        freq_ref2(w,:) = f_ref(:,(i-1)*N+1)';

        t_seg(i,1) = t(:,(i-1)*N+1);

    end

    if w==1 % ATUALIZA VALORES DA PRIMEIRA ITERAÇÃO
        
        p0(:) = phi_ref(1,:,1) - Ang_ref; 
        p1(:) = phi_ref(1,:,2) - Ang_ref;
        p2(:) = phi_ref(1,:,3) - Ang_ref;
        % ESTIMADOR FASORIAL DFT
%         [~, ~, ~, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va, Vb, Vc, Fs, f0);
%         [freq_final_fft] = estimador_freq_deriv_ang_Ain(phia(1,h+1), phib(1,h+1), phic(1,h+1), f0, h, p0(1,h), p1(1,h), p2(1,h),N,m,Fs);
%         %[freq_final_fft] = estimador_freq_deriv_ang_Ain(phi0, phi1, phi2, f0, h, p0(1,h), p1(1,h), p2(1,h),N,m,Fs);
%         
%         % atualiza angulos de entrada (ang anterior)
%         p0(1,h) = phia(1,h+1) - Ang_ref; 
%         p1(1,h) = phib(1,h+1) - Ang_ref; 
%         p2(1,h) = phic(1,h+1) - Ang_ref; 



        Fe =  freq_ref;
        %Fs = round(Fe(w)*N);   % ITERA VETOR DE FREQUÊNCIAS
    end



    %--------------------------------------------------------------------------------------------------------------------
    % ESTIMADOR FASORIAL DFT
    [~, ~, ~, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va, Vb, Vc, Fs, f0);

    granN = (Fs/N);
    gradeN = 0:granN:granN*(N-1);
    %bin = round(granN/f0)+1;       % local da componente fundamental no espectro

    %--------------------------------------------------------------------------------------------------------------------
%     % CÁLCULO DE COMPONENTES SIMÉTRICAS
%     bin = round(granN/f0)+1;       % local da componente fundamental no espectro
%     
%     pa(w,:) = phia;
%     pb(w,:) = phib;
%     pc(w,:) = phic;
%     % COM DADOS FFT
%     [A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa(:,bin)', Ab(:,bin)', Ac(:,bin)', pa(:,bin)', pb(:,bin)', pc(:,bin)');
% 
%     % ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva
%     [freq_final2] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva
    dt = m/12800;
    % COM DADOS FFT
    for h=1:length(freqs)

        [freq_final_fft] = estimador_freq_deriv_ang_Ain(phia(1,h+1), phib(1,h+1), phic(1,h+1), f0, h, p0(1,h), p1(1,h), p2(1,h),N,m,Fs);
        %[freq_final_fft] = estimador_freq_deriv_ang_Ain(phi0, phi1, phi2, f0, h, p0(1,h), p1(1,h), p2(1,h),N,m,Fs);
        %Ang_ref = 0;
        desvio_fase(1,i) = phia(1,h+1) - p0(1,h) - Ang_ref;
        desvio_fase(2,i) = phib(1,h+1) - p1(1,h) - Ang_ref;
        desvio_fase(3,i) = phic(1,h+1) - p2(1,h) - Ang_ref;
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
            % Finalmente, calcula freq. com desvios corrigidos
            freq_meas(j,i) = f0 + desvio_fase(j,i)/(2*pi*dt); 
        end
        freq_final3(i) = mean(freq_meas); % media das tres fases
        % atualiza angulos de entrada (ang anterior)
        p0(1,h) = phia(1,h+1);
        p1(1,h) = phib(1,h+1);
        p2(1,h) = phic(1,h+1);

        % Frequencia estimada
        if w==1
            Fe(w+1,h) = ref_seg(1,w+1);
        else
            %Fe(w+1,h) = freq_final2(1,2);
            Fe(w+1,h) = freq_final3;
        end
        
        
        ROCOF(w,h) = ((Fe(w+1,h)) - (Fe(w,h)))/0.04;

        Fem(w+1,h) = Fe(w+1,h);% - ROCOF(w,h)/(Fe(w+1,h));
        Fe(w+1,h) = (Fe(w+1,h) + Fe(w,h))/2;% - ROCOF(w,h)/(Fe(w+1,h));
        %Fe(w+1,h) = mean(freq_final2(1,1:3))
    end


    Amp_a(w,:) = Aa;
    Amp_b(w,:) = Ab;
    Amp_c(w,:) = Ac;

    ang_a(w,:) = phia;
    ang_b(w,:) = phib;
    ang_c(w,:) = phic;

    freq_final(1,w) = Fem(w,1);

end

end
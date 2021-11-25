


function [Aa, Ab, Ac, phia, phib, phic,freq_final] = estimador_hwfft_completo_5(Va, Vb, Vc, t, refs, N, m, Fs,f0)

%---------------------------------------------------------------------------------------%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);

% % Plot amostragem (ciclo a ciclo)
% w = 5;   % numero de janelas no plot
% figure
% for i=1:w
%     subplot(w,1,i), hold on, grid on
%     plot(t_seg(i,:),Va_seg(i,:),'o-')
% end
% %Plot referências
% figure, subplot(3,1,1), hold on, grid on
% plot(t_seg(:,1), ref_seg(1,:), 'o-')
% xlabel('tempo (s)'), ylabel('freq (Hz)'), title('Referências (segmentada)')
% subplot(3,1,2), hold on, grid on
% plot(t_seg(:,1),ref_seg(2,:), 'o-'), plot(t_seg(:,1),ref_seg(3,:), 'o-'),plot(t_seg(:,1),ref_seg(4,:), 'o-')
% xlabel('tempo (s)'), ylabel('amp (pu)')
% subplot(3,1,3), hold on, grid on, plot(t_seg(:,1),ref_seg(5,:), 'o-'),plot(t_seg(:,1),ref_seg(6,:), 'o-'),plot(t_seg(:,1),ref_seg(7,:), 'o-')
% xlabel('tempo (s)'), ylabel('fase (rad)')

% Va_seg = flip(Va_seg,2);
% Vb_seg = flip(Vb_seg,2);
% Vc_seg = flip(Vc_seg,2);

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL DFT
[~, ~, ~, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va_seg, Vb_seg, Vc_seg, Fs, f0);

granN = (Fs/N);
gradeN = 0:granN:granN*(N-1);
bin = round(granN/f0)+1;
% % Plot componente fundamental
% figure, subplot(2,1,1), hold on, grid on
% plot(Aa(:,2), 'o-'),plot(Ab(:,2), 'o-'),plot(Ac(:,2), 'o-')
% title('Amplitude (pu) - DFT'),ylabel('pu'), xlabel('Ciclo')
% subplot(2,1,2), hold on, grid on, plot(phia(:,2), 'o-'),plot(phib(:,2), 'o-'),plot(phic(:,2), 'o-')
% title('Fase (rad) - DFT'),ylabel('rad'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS

bin = round(granN/f0)+1;       % local da componente fundamental no espectro

% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa(:,bin)', Ab(:,bin)', Ac(:,bin)', phia(:,bin)', phib(:,bin)', phic(:,bin)');

% % Plot fasores da componente fundamental
% figure, subplot(2,3,1), hold on, grid on, plot(A0r, 's-'),plot(A0, 'o-')
% title('Amplitudes (pu)'),ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on, plot(rad2deg(phi0r), 's-'), plot(rad2deg(phi0), 'o-')
% title('Fase (º)'),ylabel('º'), xlabel('Ciclo')
% subplot(2,3,2), hold on, grid on, plot(A1r, 's-'), plot(A1, 'o-')
% title('Amplitudes (pu)'),ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on, plot(rad2deg(phi1r), 's-'),plot(rad2deg(phi1), 'o-')
% title('Fase (º)'),ylabel('º'), xlabel('Ciclo')
% subplot(2,3,3), hold on, grid on, plot(A2r, 's-'), plot(A2, 'o-')
% title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on, plot(rad2deg(phi2r), 's-'),plot(rad2deg(phi2), 'o-')
% title('Fase (º)'), ylabel('º'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva

% COM DADOS FFT
[freq_final] = estimador_freq_deriv_ang_saturador(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

% figure, hold on, grid on, plot(ref_seg(1,1:end), 'o-'), plot(freq_final, 'o-')
% title('Frequência (Hz)'), legend('Referência', 'Estimativa (DFT)'),ylabel('Hz'), xlabel('Ciclo')
% phia(:,bin) = phia(:,bin) + ((f0 - freq_final)*pi/f0)';
% phib(:,bin) = phib(:,bin) + ((f0 - freq_final)*pi/f0)';
% phic(:,bin) = phic(:,bin) + ((f0 - freq_final)*pi/f0)';

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL DFT
[~, ~, ~, Aa, Ab, Ac, phia, phib, phic] = estimador_hwfft(Va_seg, Vb_seg, Vc_seg, Fs, f0, freq_final);

granN = (Fs/N);
gradeN = 0:granN:granN*(N-1);
bin = round(granN/f0)+1;
% % Plot componente fundamental
% figure, subplot(2,1,1), hold on, grid on
% plot(Aa(:,2), 'o-'),plot(Ab(:,2), 'o-'),plot(Ac(:,2), 'o-')
% title('Amplitude (pu) - DFT'),ylabel('pu'), xlabel('Ciclo')
% subplot(2,1,2), hold on, grid on, plot(phia(:,2), 'o-'),plot(phib(:,2), 'o-'),plot(phic(:,2), 'o-')
% title('Fase (rad) - DFT'),ylabel('rad'), xlabel('Ciclo')


end
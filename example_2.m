clear ; clc
% computes an impedance spectrum sequentially in the presence of
% degradation and plots the output as a Nyquist plot

flag = 'A'; % use to pick recombination parameter set A or B

mu = 1e-2; % degradation rate (s-1)

f_list = logspace(6,-1,71); % list of frequencies to measure

[Z,sol] = generate_spectrum(f_list,mu,flag);


%%

figure(1)
plot(Z,'-o','MarkerFaceColor',lines(1),'MarkerSize',3)
hold on
ind = 21;
hold off
set(gca,'YDir','reverse','DataAspectRatio',[1,1,1])
xlabel('Re(z) [\Omega cm^2]')
ylabel('Im(z) [\Omega cm^2]')
set(gcf,'Position',[50,50,600,400])

%% === functions used by the above ===
function [Z, mastersol] = generate_spectrum(f_list,degrate,flag)
    tic
    rec_fac = @(t) degrate*t + 1; % construct degradation factor function

    exptime = 0; % begin at experiment time t=0
    for i = 1:length(f_list)
        tic
        params = paper_params(@(t) rec_fac(t+exptime),f_list(i),'IS',flag);
        duration = params.tstar2t(params.time(end)); % time taken to measure frequency
        
        try 
            sol{i} = numericalsolver(params); % attempt to call solver
            Z(i) = single_impedance(sol{i}); % fit impedance to current output
            fprintf('frequency %s solved in %s s \n',num2str(i), num2str(toc))
        catch me
            sol{i} = struct('time',params.tstar2t(params.time),...
                'J',nan(length(params.time),1),...
                'V',nan(length(params.time),1));
            rethrow(me)
        end
    
        exptime = exptime + duration; % progress experiment time for continuity of degradation factor
    
    end
    
    % compile master solution file by stitching together individual sol
    % files
    mastersol = sol{1};
    for i = 2:length(sol)
        mastersol.time = [mastersol.time, mastersol.time(end) + sol{i}.time(2:end)];
        mastersol.V = [mastersol.V; sol{i}.V(2:end)];
        mastersol.J = [mastersol.J; sol{i}.J(2:end)];
    end
end

function Z = single_impedance(sol)

    params = sol.params;
    
    FREQ = params.FREQ;
    nwaves = params.num_cycles;
    analysis = 1; % how many waves to analyse
    Z = [];
    ind = 1:(101+nwaves*100);
    
    % Get current
    J = sol.J(ind)*1e-3;
    t = sol.time(ind);
    
    J = J(end-analysis*100:end);
    t = t(end-analysis*100:end);
    
    try
        fit = FourierFit(t-t(1),J,FREQ);
        Vp = params.Vp;
        Jp = fit.Sp;
        theta = fit.theta+pi;
    
        Z = Vp/Jp*exp(-1i*theta); % output impedance in units of Ohm cm2
    catch me
        ind([1,end])
        rethrow(me)
    end
end

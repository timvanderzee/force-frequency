clear all; close all; clc
load('quad_parms.mat','parms','fl','fv');

% angle settings
parms.ce.amin = 1e-3; % minimal excitation
parms.exp.phi = 90;
parms = cfxc.calc_x0(parms);

% stimulation settings
parms.exp.stim_type = 'u_func';
parms.exp.A = 1;
parms.type = 'CaFaXC';
tmax = 0.3;
parms.exp.tstop = tmax;
parms.exp.pw = .005; % [s]
parms.exp.u_func = @(t, parms) parms.exp.A .* (.5 + .5 * square(2*pi*parms.exp.freq*t, parms.exp.pw*100 * parms.exp.freq)) .* (t < parms.exp.tstop);

if ishandle(1), close(1); end; figure(1)
if ishandle(2), close(2); end; figure(2)

freqs = logspace(0,3,20); 
Fpeak_freqs = nan(size(freqs));  
colors = parula(length(freqs));

% pulse widths
pws = linspace(1,10,length(freqs))/1000; % [s]

for j = 1:length(pws)
parms.exp.pw = pws(j);
    
for f = 1:length(freqs)
    
    % set frequency
    disp(['Freq = ', num2str(freqs(f)), ' Hz'])
    parms.exp.freq = freqs(f);
    
    % simulate
    [t,x] = ode113(@cfxc.sim_muscle, [0 tmax], [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)], parms.set.odeopt, parms);
    [y,X] = cfxc.get_model_output(t, x, parms);

    W = cumtrapz(-y.lce,y.Fce);
    
    figure(j)
    
    subplot(231);
    plot(t, parms.exp.u_func(t,parms),'color',colors(f,:)); hold on;
    
    subplot(232); 
    plot(t,y.Ca,'color',colors(f,:)); hold on;
    
    subplot(233); 
    plot(t,y.Fa,'color',colors(f,:)); hold on; 
    
    subplot(234)
    plot(t,y.Fce,'color',colors(f,:)); hold on
    plot(t,y.Fse,':','color',colors(f,:));
    plot(t,y.Fce+y.Fpe,'--','color',colors(f,:));

    subplot(235);
    plot(t,y.lce-y.lce(1),'color',colors(f,:)); hold on
    
    subplot(236);
    plot(t,W,'color',colors(f,:)); hold on
        
    Fpeak_freqs(f) = max(y.Fse) - min(y.Fse);
end


% make nice
titles = {'Excitation', 'Calcium activation', 'Force facilitation', 'Force','Length','Work'};
ylabels = {'u (0-1)', 'a_{Ca} (0-1)', 'a_F (0-1)', 'Force (N)', 'Length (m)', 'Work (J)'};

for i = 1:6
    subplot(2,3,i);
    box off
    title(titles{i})
    xlabel('Time (s)')
    ylabel(ylabels{i})
end

set(gcf,'units','normalized','position', [0 .3 .6 .4])

% summary figure
figure(100)
plot(freqs, Fpeak_freqs,'linewidth',2,'color',colors(j,:))
axis([0 100 0 max(Fpeak_freqs)])
box off
hold on
xlabel('Stimulation rate (Hz)'); ylabel('Peak force (N)')
title('Force - frequency')

set(gcf,'units','normalized','position', [.6 .3 .2 .4])
end


%% R2.8: frequency convergence along MatCont limit-cycle branches
% Independent reference: the PERIOD of a numerically continued limit cycle.
% Local predictor: the critical eigenvalue at its originating Hopf point.
%
% Requires the files in Honza_codes/diagram_DeltaMd:
%   EP_Md_Delta_0.5_a.mat, EP_Md_Delta_0.5_b.mat,
%   LC_Md_Delta_0.5_left_dense.mat, LC_Md_Delta_0.5_right.mat.
%
% The main validation runs from saved MatCont continuations, without MatCont
% on the MATLAB path. Set rerunFromHopf=true to demonstrate init_H_LC and a
% NEW continuation from the right-hand Hopf point (MatCont 7.6 required).
%
% All frequencies are in Hz, periods in ms, and parameter distances are
% absolute changes of M_D1 along the Delta_{v_theta}=0.5 slice.

clearvars; close all; clc
Path = 'C:\prace\MatCont7p6\Systems';
%% USER SETTINGS
scriptDir = fileparts(mfilename('fullpath'));
candidateDirs = {
    fullfile(scriptDir,'Honza_codes','diagram_DeltaMd'), ...
    fullfile(scriptDir,'tmp','honza','Honza_codes','diagram_DeltaMd'), ...
    fullfile(scriptDir,'diagram_DeltaMd')};
dataDir = '';
for k = 1:numel(candidateDirs)
    if isfile(fullfile(candidateDirs{k},'EP_Md_Delta_0.5_b.mat'))
        dataDir = candidateDirs{k};
        break
    end
end
assert(~isempty(dataDir), 'Set dataDir to your MatCont diagram_DeltaMd folder.');

maxDistance = 0.06;                 % M_D1 window around each Hopf point
targetDistances = [1e-3 1e-2 5e-2];
saveFigure = false;                 % export PDF beside this script if true
rerunFromHopf = false;              % optional new MatCont continuation

%% LOAD EQUILIBRIUM AND PERIODIC-ORBIT CONTINUATIONS
epLeft = load(fullfile(dataDir,'EP_Md_Delta_0.5_b.mat'));
epRight = load(fullfile(dataDir,'EP_Md_Delta_0.5_a.mat'));
lcLeft = load(fullfile(dataDir,'LC_Md_Delta_0.5_left_dense.mat'));
lcRight = load(fullfile(dataDir,'LC_Md_Delta_0.5_right.mat'));

% Find the H events by their recorded MatCont labels. The approximate
% coordinates only disambiguate among several H events in each EP file.
H1 = getHopfEvent(epLeft, 0.143, 'H1 (left, supercritical)');
H2 = getHopfEvent(epLeft, 0.645, 'H2 (left, subcritical)');
H3 = getHopfEvent(epRight,4.831, 'H3 (right, supercritical)');

% The left dense LC continuation starts from H2, passes through its main
% fold of cycles (LPC), and eventually reaches H1. Keep the two segments
% SEPARATE, otherwise cycles with different amplitudes could be conflated
% merely because their M_D1 coordinates are close.
folds = [lcLeft.s(strcmp(strtrim({lcLeft.s.label}),'LPC')).index];
folds = folds(folds > 10 & folds < size(lcLeft.x,2)-10);
assert(~isempty(folds),'Could not identify the interior LPC of the left branch.');
[~, imax] = max(lcLeft.x(end,folds));
foldIndex = folds(imax);

cases = struct([]);
cases(1) = buildCase(H1,lcLeft,foldIndex:size(lcLeft.x,2),maxDistance);
cases(2) = buildCase(H2,lcLeft,3:foldIndex,maxDistance);
cases(3) = buildCase(H3,lcRight,3:size(lcRight.x,2),maxDistance);

%% NUMERICAL RESULTS AND COMPACT DISTANCE TABLE
fprintf('Hopf-to-LC frequency comparison (Delta_vtheta = 0.5)\n');
fprintf('---------------------------------------------------------------------\n');
for j=1:numel(cases)
    c = cases(j);
    fprintf('%s: Md_H=%.9f, f_H=%.9f Hz, samples=%d\n', ...
        c.name,c.MdH,c.fH,numel(c.Md));
    fprintf('  requested dMd    actual |dMd|    f_LC [Hz]   |fLC-fH|/fLC [%%]\n');
    for t = targetDistances
        [~,ii] = min(abs(c.distance-t));
        fprintf('     %10.5g       %10.5g    %11.6f     %11.6g\n',...
            t,c.distance(ii),c.fLC(ii),c.errorPercent(ii));
    end
end

%% FIGURE: ACTUAL LC FREQUENCIES VS. HOPF LIMIT AND RELATIVE DIFFERENCE
fig = figure('Color','w','Position',[80 60 1160 900]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for j = 1:numel(cases)
    c = cases(j);
    [mdSort,ii] = sort(c.Md);
    nexttile(2*j-1)
    plot(mdSort,c.fLC(ii),'b.-','MarkerSize',8); hold on
    yline(c.fH,'k--','LineWidth',1.1);
    xline(c.MdH,':','Color',[.5 .5 .5]);
    xlabel('$M_{D1}$','Interpreter','latex');
    ylabel('$f_{\mathrm{LC}}\;[\mathrm{Hz}]$','Interpreter','latex');
    title(c.name,'Interpreter','none'); grid on
    if j == 1
        legend({'LC continuation','Hopf eigenfrequency'},'Location','best');
    end

    nexttile(2*j)
    good = c.distance > 0 & c.errorPercent > 0;
    loglog(c.distance(good),c.errorPercent(good),'.','MarkerSize',9,...
        'Color',[.1 .3 .7]); hold on
    xlabel('$|M_{D1}-M_{D1}^{H}|$','Interpreter','latex');
    ylabel('$100|f_{\mathrm{LC}}-f_H|/f_{\mathrm{LC}}\;[\%]$',...
        'Interpreter','latex');
    grid on
end
title(tl,'MatCont LC frequencies approaching the Hopf-onset eigenfrequency');
if saveFigure
    exportgraphics(fig,fullfile(scriptDir,'R2_8_Hopf_LC_validation.pdf'),...
        'ContentType','vector');
end

%% OPTIONAL: INITIALIZE AND CONTINUE A NEW LC BRANCH FROM H3
% This section explicitly demonstrates MatCont's init_H_LC machinery.
% It is NOT used to create the validation plot above, which comes from the
% published/archived continuation data. It may be run as an independent check.
if rerunFromHopf
    assert(exist('init_H_LC','file')==2 && exist('cont','file')==2,...
        'Add MatCont to the MATLAB path before setting rerunFromHopf=true.');
    assert(exist('Dopamine_NMM_fast4','file')==2,...
        'Add the directory containing Dopamine_NMM_fast4.m to the path.');

    % Parameter order in Dopamine_NMM_fast4:
    % [C,Iext,ujump,Md,B,gamma,Vr,Vt,Delta,cexc,Ta,K].
    par = [50,0,45,H3.MdH,0.2,1.1,-60,-43,0.5,0,6,10];
    activeParam = 4;           % M_D1
    amplitude = 1e-6;
    ntst = 150; ncol = 4;
    [x0,v0] = init_H_LC(@Dopamine_NMM_fast4,H3.state,par,...
                         activeParam,amplitude,ntst,ncol);
    opt = contset;
    opt = contset(opt,'MaxNumPoints',120);
    opt = contset(opt,'InitStepsize',1e-3);
    opt = contset(opt,'MinStepsize',1e-5);
    opt = contset(opt,'MaxStepsize',1e-2);
    opt = contset(opt,'Singularities',0);
    opt = contset(opt,'Multipliers',1);
    [xNew,~,sNew,~,~] = cont(@limitcycle,x0,v0,opt); %#ok<ASGLU>
    fprintf('New MatCont LC continuation: %d points, Md range %.6f to %.6f\n',...
        size(xNew,2),min(xNew(end,:)),max(xNew(end,:)));
end

%% LOCAL FUNCTIONS
function H = getHopfEvent(ep, approximateMd, label)
    assert(isfield(ep,'s') && isfield(ep,'f'),...
        'EP file must contain MatCont s (events) and f (eigenvalues).');
    events = ep.s(strcmp(strtrim({ep.s.label}),'H'));
    assert(~isempty(events),'No MatCont H events found.');
    mdAtEvents = arrayfun(@(e)ep.x(end,e.index),events);
    [delta,k] = min(abs(mdAtEvents-approximateMd));
    assert(delta < 0.1,'No Hopf event near requested Md=%.4f.',approximateMd);
    idx = events(k).index;
    eigvals = ep.f(:,idx);
    choices = find(imag(eigvals)>1e-10);
    assert(~isempty(choices),'No complex eigenvalue pair at MatCont H event.');
    [~,p] = min(abs(real(eigvals(choices))));
    lambda = eigvals(choices(p));
    assert(abs(real(lambda))<1e-5,'H critical pair not sufficiently close to imaginary axis.');
    H = struct('name',label,'MdH',ep.x(end,idx),...
        'fH',1000*abs(imag(lambda))/(2*pi),...
        'state',ep.x(1:4,idx));
end

function c = buildCase(H,lc,indices,maxDist)
    Md = double(lc.x(end,indices));
    T = double(lc.x(end-1,indices));
    freq = 1000 ./ T;
    d = abs(Md-H.MdH);
    valid = isfinite(Md) & isfinite(T) & T>0 & ...
        isfinite(freq) & d<=maxDist & d>1e-10;
    assert(nnz(valid)>5,'Not enough LC points near %s.',H.name);
    c = struct('name',H.name,'MdH',H.MdH,'fH',H.fH,...
        'Md',Md(valid),'fLC',freq(valid),'distance',d(valid),...
        'errorPercent',100*abs(freq(valid)-H.fH)./freq(valid));
end

%{ 
Copyright (C) 2023  N. Suhas Jagannathan
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
%}

% This script computes the oxygen consumption rate (OCR) of cells in the
% microbioreactor using the analytical solution to the ODE (6) in the main text.

function [k,t, O2_sat] = computeOCRbyAnalyticSol(data, kla)

O2_sat = 0.201; %Saturating O2 concentration at room temerature in mM O2.
O_m = O2_sat; %DO in fresh media
do_drive_frac = data(1:end,7)/100;
scale_sat = (1-do_drive_frac)+(do_drive_frac*4.762); %Computing scale factor for saturating O2 in media as a function of total O2 in the headspace

do_conc = (data(:,3)*O2_sat)/100;
scaled_do_sat = scale_sat(2:end)*O2_sat;

perf_ave = (data(1:end-1,6) + data(2:end,6))/2;
vcc_ave =  (data(1:end-1,8) + data(2:end,8))/2;
diff_t = diff(data(:,1));

comb_rate = kla+perf_ave;
e_term = exp(-comb_rate.*diff_t);

%Zt = kla*Osat + kp*Om - (kla+kp)*DO(t)
Z0 = (kla* scaled_do_sat) + (O_m*perf_ave) - (comb_rate.*do_conc(1:end-1));
Z1 = (kla* scaled_do_sat) + (O_m*perf_ave) - (comb_rate.*do_conc(2:end));

k = ((Z0.*e_term)-Z1)./(vcc_ave.*(e_term-1));
t = data(2:end,1);

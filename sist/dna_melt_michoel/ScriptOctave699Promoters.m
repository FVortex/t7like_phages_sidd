% Matlab script to reproduce the figures in the paper:
%
%   Tom Michoel and Yves van de Peer, "A helicoidal transfer matrix model
%      for inhomogeneous dna melting", Phys. Rev. E 73, 011908 (2006)
%      arXiv:q-bio/0507036
%
% Warning: it takes some time to do all these computations!

pkg load all
global dna_melt=dna_melt_preferences;


temp = 310;
sigma = -0.03;
bc = 'closed';
#seq=cmyc;
Gamma = (-0.1:0.001:0.0);
omega = (-0.04:0.0005:0.04)';

dir_name = '/home/hp/Promoters_dataset_pro_numerical_form/Promoters_dataset_pro_numerical_form/';
d=dir(dir_name)
dir_output_name='/home/hp/Probs/';
for i=1:length(d)
	if d(i).name(1)=="P"
		file = fopen(cstrcat(dir_name, d(i).name));
		cmyc = fgetl(file);
		fclose(file);
		cmyc = str2num(cmyc(:))';
		[prob,u,un] = melt_prob_lk(temp,sigma,cmyc,bc,Gamma,omega);
		matfile = fullfile(dir_output_name, cstrcat(d(i).name, '.mat'));
		save(matfile, 'prob');
	end
end



function [P, alpha, beta] = calc_tilt_correction(u, v, w, avg_interval)

	% implements the planar fit method by Steve Stage as described by Wilczak et al. (2001)
	% u, v, and w are the wind velocity components, avg_interval is the averaging interval in samples
	% if u, v and w are already averages, set avg_interval to 1
	
	if (avg_interval > 1)
		% calc mean over avg_interval points
		len = floor(length(u)/avg_interval);
		mu = zeros(len, 1);
		mv = zeros(len, 1);
		mw = zeros(len, 1);
		if (len > 0)
			for (i = 1:len)
				mu(i) = nanmean(u(((i-1)*avg_interval+1):i*avg_interval));
				mv(i) = nanmean(v(((i-1)*avg_interval+1):i*avg_interval));
				mw(i) = nanmean(w(((i-1)*avg_interval+1):i*avg_interval));
			end
		else
			P = eye(3);
			return;
		end
	else
		len = length(u);
		mu = u;
		mv = v;
		mw = w;
	end
	
	% sum of velocities
	su = nansum(mu);
	sv = nansum(mv);
	sw = nansum(mw);
	% sum of velocity products
	suv = nansum(mu.*mv);
	suw = nansum(mu.*mw);
	svw = nansum(mv.*mw);
	suu = nansum(mu.*mu);
	svv = nansum(mv.*mv);
	% matrix in equation b = H^(-1)*g
	H = [len su sv; su suu suv; sv suv svv];
	g = [sw; suw; svw];
	b = H\g;
	b1 = b(2);
	b2 = b(3);
	% rotation matrices C and D
	p31 = -b1/sqrt(b1^2 + b2^2 + 1);
	p32 = -b2/sqrt(b1^2 + b2^2 + 1);
	p33 = 1/sqrt(b1^2 + b2^2 + 1);
	sinalpha = p31;
	cosalpha = sqrt(p32^2 + p33^2);
	sinbeta = -p32/sqrt(p32^2 + p33^2);
	cosbeta = p33/sqrt(p32^2 + p33^2);
	C = [1 0 0; 0 cosbeta -sinbeta; 0 sinbeta cosbeta];
	D = [cosalpha 0 sinalpha; 0 1 0; -sinalpha 0 cosalpha];
	% final transform matrix P
	P = D'*C';
	alpha = asin(sinalpha);
	beta = asin(sinbeta);
	if (~isempty(find(isnan(P))))
		P = eye(3);
		alpha = 0;
		beta = 0;
	end
end

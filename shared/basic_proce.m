function d_ = basic_proce(d)
	
	% d is of size (nt x nr)
   
    d_ = d; % all operations in Matlab are in columns -
    
    % detrend
    %
    d_ = detrend(d_);
    
    % demean
    %
    d__mean = mean(d_);
    [nt,nr] = size(d);
    
    % taper
    %
    %tuk = tukeywin(nt , 0.1);

    for i=1:nr
    d_(:,i) = ( d_(:,i) - d__mean(i) ); % .* tuk;
    end
    
end
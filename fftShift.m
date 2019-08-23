function shifted_data = fftShift(data, samplingFreq, tdiff)
% Time shift function using the translation properties of the fft.
% This function assumes a real input signal.  
%
% Written by Ryan Schultz.
  
  % Check to make sure signal is real.
  if( ~isreal(data) )
      shifted_data=0;
      return;
  end;
  
  % Check if we're dealing with row or column vectors. 
  if( isrow(data) )
      rowflag=1;
  else
      rowflag=0;
      data=data';
  end;
  
  % Initialize.
  n=length(data);
  z=fft(data, n);
  
  % Apply shift theorem.
  wtau=-2i*pi*(samplingFreq/2)*tdiff*linspace(0,1,floor(n/2+1));
  z(1:floor(n/2+1))=exp(wtau).*z(1:floor(n/2+1));
  
  % Use the symmetry of fft on real data to 
  % assume values for the negative frequencies.
  shifted_data = ifft(z, n, 'Symmetric');
  
  % Output similarily organized vector as input.
  if(rowflag==0)
      shifted_data=shifted_data';
  end;
  
return;
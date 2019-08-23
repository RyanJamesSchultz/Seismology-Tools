function new_data = fftDeriv(data, samplingFreq, Nderiv)
  % Applies numerical differentiation/integration operators using the 
  % derivative property of the fft.  This function assumes a real input 
  % signal.  
  % The differentiation/integration order is consistent with the property;
  % Postive Nderiv are differentiation, negative are integration.
  %
  % Written by Ryan Schultz.
  
  %NOTE: Amplitudes should be fine for differentiation, maybe not
  %integration though.  Need to verify.
  
  % Check to make sure signal is real.
  if( ~isreal(data) )
      new_data=0;
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
  NFFT=2^(nextpow2(length(data))+1);
  z=fft(data, NFFT);
  
  % Apply differentiation/integration theorem.
  wtau=2i*pi*(samplingFreq/2)*linspace(0,1,floor(NFFT/2+1));
  z(1:floor(NFFT/2+1))=((wtau).^Nderiv).*z(1:floor(NFFT/2+1));
  
  % Add DC component if integrating.
  if(Nderiv<0)
      if(isinf(z(1)))
          z(1)=mean(data)*NFFT*abs(Nderiv);
      else
          z(1)=z(1)+mean(data)*NFFT*abs(Nderiv);
      end;
  end;
  
  % Use the symmetry of fft on real data to 
  % assume values for the negative frequencies.
  new_data = ifft(z, NFFT, 'Symmetric');
  new_data = new_data(1:n);
  
  % Output similarily organized vector as input.
  if(rowflag==0)
      new_data=new_data';
  end;
  
return;
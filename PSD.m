function [x,f,err] = PSD(data, dt, kinematic_flag, psd_flag, smooth_flag)
  % Function that will take input velocity time series in units of nm/s and
  % output Power Spectral Density (PSD), based on a variety of methods.  
  % Optional data processing  depends on the flag values.
  % 
  % References:
  % McNamara, D. E., & Buland, R. P. (2004). Ambient noise levels in the continental United States. Bulletin of the seismological society of America, 94(4), 1517-1527, doi: 10.1785/012003001.
  % Thomson, D. J. (1982). Spectrum estimation and harmonic analysis. Proceedings of the IEEE, 70(9), 1055-1096, doi: 10.1109/PROC.1982.12433.
  %
  % Written by Ryan Schultz.
  
  % Flag for either displacement, velocity, or acceleration PSD.
  if(strcmpi(kinematic_flag,'D'))
      x=fftDeriv(data,1/dt,-1);
  elseif(strcmpi(kinematic_flag,'V'))
      x=data;
  elseif(strcmpi(kinematic_flag,'A'))
      x=fftDeriv(data,1/dt,+1);
  else
      fprintf('Improper input for "kinematic_flag."  Aborting.\n');
      return;
  end;
  
  % Define vector lengths.
  N=length(data);
  n=floor(N/2+1);
  
  % Demean and convert from nm to m.
  x=(x-mean(x))/1e9;
  
  % Frequency axis.
  Fs=1/dt;
  f=(Fs/2)*linspace(0,1,n);
  
  % Compute PSD, and correct for taper.
  if(strcmpi(psd_flag,'FFT')) % Use simple fft method.
      taper=tukeywin(N,0.1);
      x=x.*taper;
      x=(2*dt/(taper'*taper))*abs(fft(x)).^2;
      x=x(1:n);
  elseif(strcmpi(psd_flag,'MT')) % Uses multi-taper method to estimate spectrum (Thomson, 1982).
      K=10;
      [x,f]=pmtm(x, (K+1)/2, 'onesided');
      x=x*pi/(Fs/2);
      f=f*(Fs/2)/pi;
  else % McNamara's (2004) original prescrition for PSD estimation.
      window=round(length(x)/13);
      taper=tukeywin(window,0.1);
      overlap=round(0.75*length(window));
      [x,f]=pwelch(x, taper, overlap, 'onesided');
      x=x*pi/(Fs/2);
      f=f*(Fs/2)/pi;
  end;
  
  % Units:
  % Displacement - (m/s^0)^2/Hz
  % Velocity     - (m/s^1)^2/Hz
  % Acceleration - (m/s^2)^2/Hz
  
  % Ignore zero frequency term.
  f=f(2:end);
  x=x(2:end);
  
  % Depending on input parameters, process
  if(strcmp(smooth_flag,'linear_raw'))
      return;
      
  elseif(strcmp(smooth_flag,'linear_smooth_lowess'))
      span=round(Fs);
      x=smooth(f,x,span,'lowess');
      x(x<0)=0;
      
  elseif(strcmp(smooth_flag,'linear_smooth_spline'))
      xc=fit(f,x,'smoothingspline');
      x=xc(f);
      
  elseif(strcmp(smooth_flag,'log_raw'))
      x=10*log10(x);
      
  elseif(strcmp(smooth_flag,'log_smooth')) % McNamara & Buland (2004)'s smoothing in log space.
      x=10*log10(x);
      
      % Resample PSD to be linear in log-space.
      p2end=floor( (log((Fs/2)/f(1))/log(2))*8 )/8;
      fs=(Fs/2)*2.^(-p2end:0.125:0);
      xs=zeros(size(fs));
      err=xs;
      
      % Smooth the PSD by taking full-octave averages in 1/8 octave intervals.
      for i=1:length(fs)
          xs(i)=mean(  x((f<(fs(i)*sqrt(2)))&(f>(fs(i)/sqrt(2))))  );
          err(i)=std(  x((f<(fs(i)*sqrt(2)))&(f>(fs(i)/sqrt(2))))  );
      end;
      f=fs;
      x=xs;
      
  elseif(strcmp(smooth_flag,'log_smooth_lowess'))
      x=10*log10(x);
      
      % Resample PSD to be linear in log-space.
      p2end=(ceil( (log((Fs/2)/f(1))/log(2))*8 )/8)-2;
      fs=(Fs/2)*2.^(-p2end:0.125:0);
      xs=zeros(size(fs));
      
      % Smooth the PSD in full-octave linear fits, by 1/8 octave intervals.
      for i=1:length(fs)
          ftemp=f((f<(fs(i)*sqrt(2)))&(f>(fs(i)/sqrt(2))));
          xtemp=x((f<(fs(i)*sqrt(2)))&(f>(fs(i)/sqrt(2))));
          xtemp=smooth( ftemp, xtemp, length(xtemp), 'lowess'  );
          xs(i)=interp1(ftemp,xtemp,fs(i),'linear');
      end;
      f=fs;
      x=xs;
      
  elseif(strcmp(smooth_flag,'log_smooth_spline'))
      x=10*log10(x);
      
      % Resample PSD to be linear in log-space.
      p2end=(ceil( (log((Fs/2)/f(1))/log(2))*8 )/8)-2;
      fs=(Fs/2)*2.^(-p2end:0.125:0);
      xs=zeros(size(fs));
      
      % Smooth the PSD in full-octave linear fits, by 1/8 octave intervals.
      for i=1:length(fs)
          ftemp=f((f<(fs(i)*sqrt(2)))&(f>(fs(i)/sqrt(2))));
          xtemp=x((f<(fs(i)*sqrt(2)))&(f>(fs(i)/sqrt(2))));
          xc=fit(ftemp',xtemp,'smoothingspline');
          xs(i)=xc(fs(i));
      end;
      f=fs;
      x=xs;
      
  else
      fprintf('Improper input for "process_flag."  Aborting.\n');
      return;
  end;
  
return;